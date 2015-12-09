# TODO:
# * special handling of Proline?
# * implement trajectory extension
# * implement collection of best structures
# * implement SAXS

from __future__ import print_function

import sys, time, os, glob, random, string, subprocess
import numpy as np
import scipy.constants as const
import mdtraj as md
import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit


class MCtrj(md.Trajectory):
    """
    Extends mdtraj.Trajectory by functionality to:
        * perform dihedral MC steps
        * minimize energy with OpenMM
        * compute Epot    with OpenMM
        * compute SAXS data fit (chi^2)
        * metroplolis criterion
    """

    def __init__(self, xyz=None, topology=None, time=None, unitcell_lengths=None, unitcell_angles=None, other=None, pdbfile="", logfile="mcmin.log", saxsfilename=None):

        # if created like the parent class
        if xyz is not None and topology is not None:
            try:
                super(MCtrj, self).__init__(xyz, topology, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)
            except TypeError as e:
                print(MCtrj)
                print(type(self))
                raise e

        # if created from parent class instance
        elif other is not None:
            super(MCtrj, self).__init__(other.xyz, other.top, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

        self.pdbfile                   = pdbfile
        self.logfilename               = logfile
        self.saxsfilename              = saxsfilename
        self.openmm_platform_is_set_up = False     # keep track of what has been initialised
        self.openmm_is_set_up          = False
        self.MC_is_set_up              = False     
        self.MC_attempt_counter        = 0         # counts the number of attempted MC moves
        self.MC_step_counter           = 0         # counts the number of accepted MC
        self.MC_moves                  = []        # holding all the MC moves
        self.resids                    = np.array([r.resSeq for r in self.top.residues], dtype=np.int)
        self.reporters                 = []

        # data storage
        self.Epot      = np.zeros(self.n_frames, dtype=np.float)
        self.SAXS_chi2 = np.zeros(self.n_frames, dtype=np.float)
        self.E_SAXS    = np.zeros(self.n_frames, dtype=np.float)

        # list of attributes added elsewhere
        self.angle_std  = 0.0  # maximum rotation angle with unit
        self.nrot       = 0    # number of angles to rotate each step

        # constant with units
        self.R = const.R /1e3 * unit.kilojoule/unit.kelvin/unit.mole  # gas constant

        # timer
        self.timer_run        = 0.0
        self.timer_MC_moves   = 0.0
        self.timer_minimize   = 0.0
        self.timer_metropolis = 0.0 
 

    def allocate_memory(self, n=None):
        """
        Allocate more memory for n new frames (double it by default)
        """
        if n is None:
            n = self.n_frames
        oldframes = self.n_frames
        newframes = self.n_frames + n

        xyz                 = np.zeros([newframes, self.n_atoms, 3], dtype=self.xyz.dtype)
        xyz[:oldframes,:,:] = self.xyz
        self.xyz            = xyz

        if len(self.time) > 1:
            dt              = self.time[1] - self.time[0]
        else:
            dt              = 1.0
        time                = np.zeros([newframes], dtype=self.time.dtype)
        time[:oldframes]    = self.time
        time[oldframes:]    = self.time[-1] + dt + dt*np.arange(newframes-oldframes)
        self.time           = time

        Epot                = np.zeros([newframes], dtype=self.Epot.dtype)
        Epot[:oldframes]    = self.Epot
        self.Epot           = Epot

        SAXS_chi2             = np.zeros([newframes], dtype=self.SAXS_chi2.dtype)
        SAXS_chi2[:oldframes] = self.SAXS_chi2
        self.SAXS_chi2        = SAXS_chi2 

        E_SAXS                = np.zeros([newframes], dtype=self.E_SAXS.dtype)
        E_SAXS[:oldframes]    = self.E_SAXS
        self.E_SAXS           = E_SAXS

        if self.unitcell_lengths is not None:
            lengths               = np.zeros([newframes,3], dtype=self.unitcell_lengths.dtype)
            lengths[:oldframes,:] = self.unitcell_lengths
            lengths[oldframes:,:] = self.unitcell_lengths[-1,:]
            self.unitcell_lengths = lengths

        if self.unitcell_angles is not None:
            angles               = np.zeros([newframes,3], dtype=self.unitcell_angles.dtype)
            angles[:oldframes,:] = self.unitcell_angles
            angles[oldframes:,:] = self.unitcell_angles[-1,:]
            self.unitcell_angles = angles


    def openmm_setup(self, ff='amber99sbildn.xml', water='amber99_obc.xml',
            T=300.0*unit.kelvin, friction=1.0/unit.picoseconds, dt=1e-8*unit.femtoseconds, nbCutoff=1.0*unit.nanometer, emtol=1.0*unit.kilojoule_per_mole,
            platform='OpenCL', properties={'OpenCLPrecision': 'mixed', 'OpenCLDeviceIndex': '0'}):
#            platform='CUDA', properties = {'CudaPrecision': 'mixed'}):
        """
        Setup up OpenMM for energy minimization and potential energy evaluation

        Parameters
        ----------
        ff (string):   OpenMM force field XML filename
        water(string): OpenMM water model XML filename
        T:             Temperature in unit.kelvin
        friction:      Langevin friction coefficient in unit.picoseconds**-1
        dt:            Integration timestep in unit.femtoseconds
        nbCutoff       non bonded cutoff in unit.nanometer
        platform:      OpenMM platform to use
        properties:    Properties (options) of the OpenMM platform
        """

        self.openmm_T        = T
        self.openmm_friction = friction
        self.openmm_dt       = dt
        self.openmm_nbCutoff = nbCutoff
        self.openmm_emtol    = emtol

        # create openmm topology
        self.openmm_topology = self.top.to_openmm() 

        # OpenMM initialization
        self.openmm_platform   = mm.Platform.getPlatformByName(platform)
        self.openmm_forcefield = app.ForceField(ff, water)
        self.openmm_system     = self.openmm_forcefield.createSystem(self.openmm_topology,
                bondedMethod=app.CutoffPeriodic, nonbondedCutoff=nbCutoff, constraints=app.HBonds)
        self.openmm_integrator = mm.LangevinIntegrator(T, friction, dt) 
        self.openmm_simulation = app.Simulation(self.openmm_topology, self.openmm_system, self.openmm_integrator, self.openmm_platform, properties)

        self.openmm_is_set_up = True



    def MC_setup(self, seed=None, logfile="mcmin.log", saxsfilename=None, **kwargs):
        """
        Set up Monte Carlo parameters

        Parameters:
        -----------
        seed (int):       random number generator seed
        logfile (string): filename for loging information
        kwargs:           additional keyword arguments are passed on to openmm_setup()
        """

        if type(seed) == type(int()):
            np.random.seed(seed)

        self.MC_attempt_counter = 0
        self.MC_step_counter    = 0
        self.MC_is_set_up       = True

        # also make sure OpenMM is set up
        if not self.openmm_is_set_up:
            self.openmm_setup(**kwargs)

        # minimize the energy of the initial structure
        self.minimize_energy()

        # set SAXS data file
        if saxsfilename is not None:
            self.saxsfilename = saxsfilename


            

    def add_ramachandran_moves(self, philist=None, psilist=None, angle_std=120*unit.degree, nrot=1, verbose=True):
        """
        Add dihedral rotation moves for the backbone phi/psi angles
        By default, ramachandran moves for all rotatable residues are added.

        Parameters:
        -----------
        philist:    residue IDs of phi angles to rotate (zero based residue numbering)
        psilist:    residue IDs of psi angles to rotate (zero based residue numbering)
        angle_std:  maximum rotation angle (with unit)
        nrot (int): number of rotations per step 
        """
        if philist is None:
            nphidihedrals = md.compute_phi(self[:1])[0].shape[0]
            philist    = range(1, nphidihedrals)
        if psilist is None:
            npsidihedrals = md.compute_psi(self[:1])[0].shape[0]
            psilist    = range(0, npsidihedrals)

        # add MC moves
        ndihedrals = len(philist) + len(psilist)
        frequency  = max(min(1.0, 1.0 * nrot / ndihedrals), 0.0)
        for residue in philist:
            self.MC_moves.append(MC_ramachandran_move(self, frequency=frequency, residue=residue, kind='phi', angle_std=angle_std, verbose=verbose))
        for residue in psilist:
            self.MC_moves.append(MC_ramachandran_move(self, frequency=frequency, residue=residue, kind='psi', angle_std=angle_std, verbose=verbose))
 



    def MC_run(self, nsteps=None, nattempts=None, verbose=False, **kwargs):
        """
        Run a monte carlo simulation

        Parameters
        ----------
        nsteps (int):                           number of accepted monte carlo steps to perform
        nattempts (int):                        number of monte calro attempts to try (overrides nsteps)
        kwargs:                                 keyword arguments passed on to MC_step()
        """

        if not self.MC_is_set_up:
            self.MC_setup(**kwargs)

        # reset timer
        self.timer_run        = 0.0
        self.timer_MC_moves   = 0.0
        self.timer_minimize   = 0.0
        self.timer_saxs       = 0.0
        self.timer_metropolis = 0.0
        starttime             = time.time()

        properties  = ["step", "attempt", "Etot (kJ/mol)", "E_SAXS (kJ/mol)", "Epot (kJ/mol)", "dEtot (kJ/mol)", "dE_SAXS (kJ/mol)", "dEpot (kJ/mol)", "ap", "accepted"]
        fieldwidths = [     8,         8,              14,                14,              14,               14,                 14,               14,   14,          8]
        self.reporters = []
        filereporter = StateReporter(properties=properties, fieldwidths=fieldwidths, outfile=self.logfilename)
        self.reporters.append(filereporter)
        if verbose:
            stdoutreporter = StateReporter(properties=properties, fieldwidths=fieldwidths, )
            self.reporters.append(stdoutreporter)
        
        # try this many MC attempts
        if type(nattempts) == type(int()):
            for n in range(nattempts):
                success = self.MC_step(verbose=verbose)
                # if successfull, allocate memory for one more frame
                if success and n < nattempts-1:
                    self.allocate_memory(1)

        # or perform this many successfull steps
        elif type(nsteps) == type(int()):
            # allocate enough memory
            framesNeeded = (nsteps + self.MC_step_counter+1) - self.n_frames
            if framesNeeded > 0:
                self.allocate_memory(framesNeeded)
 
            while self.MC_step_counter < nsteps:
                self.MC_step(verbose=verbose)

        # else abort
        else:
            print("Either nsteps or nattempts needs to be an integer!")
            sys.exit(1)


        self.timer_run = time.time() - starttime


        # final report
        print("")
        print("================================================================================")
        print("Statistics:")
        print("-----------")
        print("runtime:               {:12.2f} sec.".format(self.timer_run))
        print("# MC steps:            {:12d}".format(self.MC_step_counter))
        print("# MC attempts:         {:12d}".format(self.MC_attempt_counter)) 
        if self.MC_step_counter > 0:
            print("mean attempts/step:    {:12.2f}".format(1.0 * self.MC_attempt_counter / self.MC_step_counter))
            print("mean time/step:        {:12.2f} sec.".format(self.timer_run / self.MC_step_counter))
        if self.MC_attempt_counter > 0:
            print("mean time/attempt:     {:12.2f} sec.".format(self.timer_run / self.MC_attempt_counter))
        if self.timer_run > 0:
            print("time for MC moves:     {:12.2f} sec.  ({:5.2f}%)".format(self.timer_MC_moves,   100*self.timer_MC_moves  /self.timer_run))
            print("time for minimization: {:12.2f} sec.  ({:5.2f}%)".format(self.timer_minimize,   100*self.timer_minimize  /self.timer_run))
            if self.saxsfilename is not None:
                print("time for SAXS fit:     {:12.2f} sec.  ({:5.2f}%)".format(self.timer_saxs,   100*self.timer_saxs/self.timer_run))
            print("time for metropolis:   {:12.2f} sec.  ({:5.2f}%)".format(self.timer_metropolis, 100*self.timer_metropolis/self.timer_run))
        print("================================================================================")






    def MC_step(self, verbose=False, **kwargs):
        """
        Perform a single monte carlo step

        Parameters
        ----------
        verbose (bool): print status messages
        kwargs:         keyword arguments passed on to MC_setup()

        Returns
        -------
        success (bool): True if the new position has been accepted
        """
        if not self.MC_is_set_up:
            print("Call MC_setup() before MC_step().", file=sys.stderr)
            sys.exit(1)

        self.MC_attempt_counter += 1
        self.MC_step_counter    += 1

        # if we don't have any space left, allocate more memory (double it)
        if self.MC_step_counter >= self.n_frames:
            self.allocate_memory()

        # perform Monte Carlo moves
        st = time.time()
        self.xyz[self.MC_step_counter,:,:] = self.xyz[self.MC_step_counter-1,:,:]
        for move in self.MC_moves:
            move.move(self.xyz[self.MC_step_counter,:,:])
        self.superpose_chains()
        self.timer_MC_moves += time.time() - st

        # minimize energy
        st = time.time()
        self.minimize_energy()
        Etot      = self.Epot[self.MC_step_counter]
        Etot_prev = self.Epot[self.MC_step_counter-1]
        self.timer_minimize += time.time() - st


        # compute SAXS fit
        if self.saxsfilename is not None:
            st = time.time()
            self.SAXS_energy()
            Etot      += self.E_SAXS[self.MC_step_counter]
            Etot_prev += self.E_SAXS[self.MC_step_counter-1]
            self.timer_saxs += time.time() - st

        # metropolis criterion
        st = time.time()
        deltaE = (Etot - Etot_prev) * unit.kilojoule / unit.mole
        beta   = self.openmm_T * self.R
        acceptance_probability = 1.0
        if deltaE > 0.0 * unit.kilojoule / unit.mole:
            acceptance_probability = np.exp(- deltaE / beta)

        if np.random.rand() > acceptance_probability:
            # reject
            accepted = "rej"
            success  = False
        else:
            accepted = "->ACC"
            success = True
        self.timer_metropolis += time.time() - st

    
        # report status
        data = {"step":             self.MC_step_counter,
                "attempt":          self.MC_attempt_counter,
                "Etot (kJ/mol)":    Etot,
                "E_SAXS (kJ/mol)":  self.E_SAXS[self.MC_step_counter], 
                "Epot (kJ/mol)":    self.Epot[self.MC_step_counter],
                "dEtot (kJ/mol)":   Etot - Etot_prev,
                "dE_SAXS (kJ/mol)": self.E_SAXS[self.MC_step_counter] - self.E_SAXS[self.MC_step_counter-1],
                "dEpot (kJ/mol)":   deltaE.value_in_unit(deltaE.unit),
                "ap":               acceptance_probability,
                "accepted":         accepted}
        for reporter in self.reporters:
            reporter.report(data)

        if not success:
            self.MC_step_counter -= 1

        return success



    def superpose_chains(self):
        """
        Superimpose each chain on the previous frame
        """
        frame_indices = np.arange(2) + self.MC_step_counter - 1
        try:
            last2frames = self.slice(frame_indices, copy=False)
        except UnboundLocalError as e:
            print("There is a bug in trajectory.slice(). A pull request as already been submitted.") 
            sys.exit(1)

        for chain in self.top.chains:
            atom_indices = self.top.select("chainid {}".format(chain.index))
            thischain    = last2frames.atom_slice(atom_indices, inplace=False)
            thischain.superpose(thischain, frame=0, parallel=False)
            self.xyz[self.MC_step_counter, atom_indices, :] = thischain.xyz[-1,:,:]


#        for chain in self.top.chains:
#            atom_indices = self.top.select("chainid {}".format(chain.index))
#            try:
#                thisframe    = self.slice(self.MC_step_counter, copy=False)
#            except UnboundLocalError as e:
#                print("There is a bug in trajectory.slice(). A pull request as already been submitted.")
#            thisframe.superpose(self, frame=self.MC_step_counter-1, atom_indices=atom_indices, parallel=False)
#            thischain.
#            break



    def minimize_energy(self, frame=None):
        """
        Energy minimize the given frame or the current frame if not specified
        (as specified by self.MC_step_counter)
        """
        # set frame
        if frame is None:
            frame = self.MC_step_counter

        # set coordinates
        self.openmm_simulation.context.setPositions(self.xyz[frame])

        # evaluate energy
        state = self.openmm_simulation.context.getState(getEnergy=True)
        Epot_before_em = state.getPotentialEnergy()

        # save frame
        self.openmm_simulation.saveCheckpoint('openmm.cp')
        self.openmm_simulation.reporters = []
        self.openmm_simulation.reporters.append(app.PDBReporter('em/em_{:05}_b.pdb'.format(self.MC_attempt_counter), 1))
        self.openmm_simulation.step(1)
        self.openmm_simulation.reporters = []
        self.openmm_simulation.loadCheckpoint('openmm.cp')
 
        # minimize energy
        self.openmm_simulation.minimizeEnergy(tolerance=self.openmm_emtol) 
        state = self.openmm_simulation.context.getState(getPositions=True, getEnergy=True)

        # save frame
        self.openmm_simulation.saveCheckpoint('openmm.cp')
        self.openmm_simulation.reporters = []
        self.openmm_simulation.reporters.append(app.PDBReporter('em/em_{:05}_a.pdb'.format(self.MC_attempt_counter), 1))
        self.openmm_simulation.step(1)
        self.openmm_simulation.reporters = []
        self.openmm_simulation.loadCheckpoint('openmm.cp')
 
        # get potential energy
        Epot             = state.getPotentialEnergy()
        self.Epot[frame] = Epot.value_in_unit(unit.kilojoule/unit.mole)

        # get minimized coordinates
        pos = state.getPositions()
        self.xyz[frame,:,:] = np.array(pos.value_in_unit(unit.nanometer), dtype=self.xyz.dtype)



    def compute_Epot(self, frame=None):
        """
        Evaluate the potential energy of the given frame (self.MC_step_counter if unspecified)
        """

        # set frame
        if frame is None:
            frame = self.MC_step_counter

        # set coordinates
        self.openmm_simulation.context.setPositions(self.xyz[frame])

        # evaluate energy
        state = self.openmm_simulation.context.getState(getEnergy=True)
        Epot  = state.getPotentialEnergy()

        self.Epot[frame] = Epot.value_in_unit(unit.kilojoule/unit.mole)



    def get_rotation_angles(self):
        """
        Generate random rotation angles for phi and psi angles
        """
        # select angles
        phi_indices = self.philist[np.random.randint(len(self.philist), size=np.random.randint(self.nrot+1))]
        psi_indices = self.psilist[np.random.randint(len(self.psilist), size=self.nrot-len(phi_indices)  )]

        # draw random rotation angles
        phi_angles = self.angle_std * np.random.randn(len(phi_indices))
        psi_angles = self.angle_std * np.random.randn(len(psi_indices))

        # add up angles that have been selected multiple times
        if len(phi_indices) > 0:
            phi_indices_dict = {p: 0.0*self.angle_std.unit for p in set(phi_indices)}
            for ndx, phi_ndx in enumerate(phi_indices):
                phi_indices_dict[phi_ndx] += phi_angles[ndx]
            phi_indices, phi_angles = zip(*phi_indices_dict.iteritems())
        if len(psi_indices) > 0:
            psi_indices_dict = {p: 0.0*self.angle_std.unit for p in set(psi_indices)}
            for ndx, psi_ndx in enumerate(psi_indices):
                psi_indices_dict[psi_ndx] += psi_angles[ndx] 
            psi_indices, psi_angles = zip(*psi_indices_dict.iteritems())

        return phi_indices, phi_angles, psi_indices, psi_angles


    def rotate(self, phi_indices, phi_angles, psi_indices, psi_angles):
        """
        Rotate dihedral angles
        """

        # copy previous coordinates
        self.xyz[self.MC_step_counter,:,:] = self.xyz[self.MC_step_counter-1,:,:] 

        for phi_index, angle in zip(phi_indices, phi_angles):
            residue     = self.top.residue(phi_index)
            resid       = residue.resSeq
            N_atom      = residue.atoms_by_name('N').next()
            N_index     = N_atom.index
            CA_atom     = residue.atoms_by_name('CA').next()
            CA_index    = CA_atom.index
            point       = self.xyz[self.MC_step_counter, N_index , :]
            vector      = self.xyz[self.MC_step_counter, CA_index, :] - point
            vector     /= np.linalg.norm(vector)
            M           = rotation_matrix_arbitrary_axis(point, vector, angle)
            rotIndices  = list(self.top.select("residue {} and (sidechain or name C or name O) and not name H".format(resid)))
            rotIndices += list(self.top.select("residue > {}".format(resid))) 

            rotvec = np.ones([4,1], dtype=self.xyz.dtype)
            for rotIdx in rotIndices:
                rotvec[:3] = self.xyz[self.MC_step_counter, rotIdx, :].reshape([3,1])
                self.xyz[self.MC_step_counter, rotIdx, :] = np.dot(M, rotvec)[:3,0]

        for psi_index, angle in zip(psi_indices, psi_angles):
            residue     = self.top.residue(psi_index)
            resid       = residue.resSeq
            CA_atom     = residue.atoms_by_name('CA').next()
            CA_index    = CA_atom.index
            C_atom      = residue.atoms_by_name('C').next()
            C_index     = C_atom.index
            point       = self.xyz[self.MC_step_counter, CA_index, :]
            vector      = self.xyz[self.MC_step_counter, C_index , :] - point
            vector     /= np.linalg.norm(vector)
            M           = rotation_matrix_arbitrary_axis(point, vector, angle)
            rotIndices  = list(self.top.select("residue {} and name O".format(resid)))
            rotIndices += list(self.top.select("residue > {}".format(resid))) 

            rotvec = np.ones([4,1], dtype=self.xyz.dtype)
            for rotIdx in rotIndices:
                rotvec[:3] = self.xyz[self.MC_step_counter, rotIdx, :].reshape([3,1])
                self.xyz[self.MC_step_counter, rotIdx, :] = np.dot(M, rotvec)[:3,0]



    def center(self):
        """
        Center each frame of the trajectory
        """
        for f in range(self.n_frames):
            self.xyz[f,:,:] -= self.xyz[f,:,:].mean(0)


    def SAXS_energy_function(self, chi2, depth=10.0*unit.kilojoule_per_mole, width=5.0):
        """
        Morse potential like function to convert chi**2 to an energy in kJ/mol
        """
        energy = depth * (1.0 - np.exp(-chi2/width))**2 
        return energy


    def SAXS_energy(self, saxsfilename=None, frame=None, verbose=False):
        """
        Compute SAXS fit and convert chi**2 to potential Energy
        """

        if frame is None:
            frame = self.MC_step_counter 
        if saxsfilename is None:
            saxsfilename = self.saxsfilename 

        chi2 = self.compute_SAXS_chi2(saxsfilename=saxsfilename, frame=frame, verbose=verbose)

        # convert the chi**2 value of the SAXS fit to an energy
        energy = self.SAXS_energy_function(chi2)

        self.SAXS_chi2[frame] = chi2
        self.E_SAXS[frame]    = energy.value_in_unit(unit.kilojoule_per_mole)




    def compute_SAXS_chi2(self, saxsfilename=None, frame=None, verbose=False):
        """
        Compute the SAXS fit of the specified frame to the given experimental scattering data.
        Fitting is done with the program Crysol:
        Svergun D.I., Barberato C. & Koch M.H.J. (1995) J. Appl. Cryst. 28, 768-773.
        http://www.embl-hamburg.de/biosaxs/manuals/crysol.html

        Parameters:
        -----------
        saxsfilename (string): filename of the SAXS scattering data, defaults to self.saxs_filename
        frame (int):           which frame to use, defaults to MC_step_counter

        Returns
        -------
        chi2 (float):          the chi2 value after fitting the theoretical scattering data to the experimental curve
        """

        randomstring   = ''.join(random.choice(string.lowercase+string.uppercase+string.digits) for i in range(4))
        crysol_prefix  = "crysol_tmp_{}".format(randomstring)
        tmppdb         = "{}.pdb".format(crysol_prefix)
        crysol_outfile = "{}.out".format(crysol_prefix)

        if frame is None:
            frame = self.MC_step_counter
        if saxsfilename is None:
            saxsfilename = self.saxsfilename

        # write temporary pdb file
        trj = md.Trajectory(self.xyz[frame,:,:], self.topology)
        trj.save_pdb(tmppdb)

        # launch Crysol to compute SAXS fit
        crysol_cmd  = ["crysol", tmppdb, saxsfilename, "-cst", "-err", "-p {}".format(crysol_prefix)]
        if verbose:
            print(' '.join(crysol_cmd))
        with open(crysol_outfile, 'w') as of:
            subprocess.call(crysol_cmd, stdout=of, shell=False)

        # read chi**2
        with open(crysol_prefix + "-100.fit", 'r') as f:
            data = f.readline()
            try:
                chi2 = float(data.split(':')[-1])
            except ValueError:
                chi2 = float('inf')

        # remove temporary files
        for f in glob.glob(crysol_prefix + '*'):
            os.remove(f)
        if os.path.isfile('crysol_summary.txt'):
            os.remove('crysol_summary.txt')

        return chi2


    
def rotation_matrix_arbitrary_axis(r, v, t):
    """
    Rotation matrix around an arbitrary axis

    Parameters
    ----------
    r: array, (3,)
        point on the axis
    v: array, (3,)
        axis direction (normalized)
    t: float
        rotation angle (in rad or symtk.unit)

    Returns
    -------
    M : array(4,4)
        4 dimensionl rotation and translation matrix
    """

    M = np.zeros([4,4], dtype=np.float64)
    a, b, c = r
    u, v, w = v
    try:
        cost = np.cos(t.value_in_unit(unit.radian))
        sint = np.sin(t.value_in_unit(unit.radian))
    except AttributeError:
        cost = np.cos(t)
        sint = np.sin(t) 
    onemincost = 1 - cost

    M[0,0] = u**2 + (v**2 + w**2)*cost
    M[1,1] = v**2 + (u**2 + w**2)*cost
    M[2,2] = w**2 + (u**2 + v**2)*cost

    M[0,1] = u*v*onemincost - w*sint
    M[0,2] = u*w*onemincost + v*sint
    M[1,0] = u*v*onemincost + w*sint
    M[1,2] = v*w*onemincost - u*sint
    M[2,0] = u*w*onemincost - v*sint
    M[2,1] = v*w*onemincost + u*sint

    M[0,3] = (a*(v**2 + w**2) - u*(b*v + c*w))*onemincost + (b*w - c*v)*sint
    M[1,3] = (b*(u**2 + w**2) - v*(a*u + c*w))*onemincost + (c*u - a*w)*sint
    M[2,3] = (c*(u**2 + v**2) - w*(a*u + b*v))*onemincost + (a*v - b*u)*sint

    M[3,:] = np.array([0., 0., 0., 1.])

    return M
 











class MC_move(object):
    """
    Monte Carlo move base class.
    All possible monte carlo moves should derive from this class.
    Derived classes can all expose their own initializers,
    but must not change the signature of the move() method.

    Parameters:
    -----------
    freqency (float): frequency with which this move should be tried [0, 1]
    topology (mdtraj.topology): MDtraj topology
    """

    def __init__(self, topology, frequency=1.0):
        self.topology  = topology
        self.frequency = frequency


    def move(self, xyz, rand=None):
        """
        Perform the Monte Carlo move.

        Parameters:
        -----------
        xyz (numpy array of shape [n_atoms, 3]): coordinates from mdtraj trajectory. Change in place.
        rand (float):                            random number from interval [0,1) to decide if the move is to be performed
        """ 

        pass




class MC_ramachandran_move(MC_move):
    """
    Backbone dihedral (ramachandran) Monte Carlo move.

    Parameters:
    -----------
    trajectory (mdtraj.Trajectory):                MDtraj trajectory
    freqency (float):                       frequency with which this move should be tried [0, 1]
    angle_ndx (int, optional):              The angle needs to be specified with either the angle_ndx (0 is the first phi/psi angle),
                                            the residue ID (resid, needs to be unique if specified), or the residue sequence number (residue).
                                            Precedence order: angle_ndx > residue > resid
    residue (int, optional):                see above
    resid (int, optional):                  see above
    kind (string) :                         'phi' or 'psi'
    angle_std (unit.degree or unit.radian): standard deviation of the normally distributed random rotation angles
    """

    def __init__(self, trajectory, frequency=1.0, angle_ndx=None, residue=None, resid=None, kind='phi', angle_std=1.0*unit.degree, verbose=False):

        assert frequency >= 0.0 and frequency <= 1.0,                         "frequency must be in interval [0,1]"
        assert kind in ['phi', 'psi'],                                        "kind must be 'phi' or 'psi'"
        assert type(angle_std) == unit.quantity.Quantity,                     "angle_std must be of unit.degree or unit.radian"
        assert sum([i is not None for i in [angle_ndx, residue, resid]]) > 0, "One of [angle_ndx, residue, resid] must be an integer"

        # parent class initializer
        super(MC_ramachandran_move, self).__init__(trajectory.topology, frequency)

        self.verbose   = verbose
        self.kind      = kind
        self.angle_std = angle_std

        # get list of all ramachandran atom indices
        if self.kind == "phi":
            dihedral_atom_indices = md.compute_phi(trajectory[:1])[0]
        elif self.kind == "psi":
            dihedral_atom_indices = md.compute_psi(trajectory[:1])[0]
        ndihedrals = dihedral_atom_indices.shape[0] # number of dihedral angles

        # determine index of this dihedral angle
        if angle_ndx is not None:
            self.angle_ndx  = angle_ndx 
        elif residue is not None:
            residue_indices = [self.topology.atom(a).residue.index for a in dihedral_atom_indices[:,1]]
            self.angle_ndx  = residue_indices.index(residue)
        elif resid is not None:
            residue_ids = [self.topology.atom(a).residue.resSeq for a in dihedral_atom_indices[:,1]]
            if residue_ids.count(resid) != 1:
                self.angle_ndx = residue_ids.index(resid)
            else:
                print("{} angle not (unambiguoulsy) found for resid {}".format(self.kind, resid), file=sys.stderr)
                sys.exit(1)

        assert self.angle_ndx in range(ndihedrals), "dihedral angle index out of range"
        
        # determine indices of atoms spanning the rotation vector
        if self.kind == 'phi':
            atom1_name = 'N'
            atom2_name = 'CA'
        elif self.kind == 'psi':
            atom1_name = 'CA'
            atom2_name = 'C'
        self.atom1 = dihedral_atom_indices[self.angle_ndx,[self.topology.atom(a).name for a in dihedral_atom_indices[self.angle_ndx,:]].index(atom1_name)]
        self.atom2 = dihedral_atom_indices[self.angle_ndx,[self.topology.atom(a).name for a in dihedral_atom_indices[self.angle_ndx,:]].index(atom2_name)]

        self.residue = self.topology.atom(self.atom1).residue
        chainid      = self.residue.chain.index

        # determine atoms to rotate
        selection_texts = []
        if self.kind == 'phi':
            selection_texts.append("resid {} and (sidechain or name C or name O) and not name H and chainid {}".format(self.residue.index, chainid))
            selection_texts.append("resid > {} and chainid {}".format(self.residue.index, chainid))
        elif self.kind == 'psi':
            selection_texts.append("resid {} and name O and chainid {}".format(self.residue.index, chainid))
            selection_texts.append("resid > {} and chainid".format(self.residue.index))
 
        self.rotation_atoms = []
        for text in selection_texts:
            self.rotation_atoms += list(self.topology.select(text))

        self.symmetry_partners = []

        # report
        if self.verbose:
            print("MC move: {:3s} {:4d} {:3s} dihedral".format(self.residue.name, self.residue.resSeq, self.kind))
#            for atom_ndx in self.rotation_atoms:
#                print("\t{}".format(str(self.topology.atom(atom_ndx))))

        

    def add_symmetry_partner(self, other):
        """
        Add a symmetry partner to the list of symmetry partners.
        The partners will always be rotated by the same angle.
        This dihedrals should be added to the partner as well in case the parnter is also in the
        list of MC moves.

        Parameters
        ----------
        other (MC_ramachandran_move): The partner angle
        """
        self.symmetry_partners.append(other)


    def move(self, xyz, rand=None, angle=None):
        """
        Perform the Monte Carlo move.

        Parameters:
        -----------
        xyz (numpy array of shape [n_atoms, 3]): coordinates from mdtraj trajectory. Change in place.
        rand (float):                            random number from interval [0,1) to decide if the move is to be performed
        angle (unit.degree or unit.radian):      definite angle to rotate (choose randomly of None), float or int are also fine
        """

        if rand is None:
            rand = np.random.rand()

        if rand >= self.frequency:
            return

        if angle is None:
            angle = self.angle_std * np.random.randn()
        else:
            if type(angle) == unit.quantity.Quantity:
                pass
            elif type(angle) == type(float()) or type(angle) == type(int()):
                anlge = angle * unit.degree
            else:
                print("Angle must be unit.degree, unit.radian, float or int.", file=sys.stderr)
                sys.exit(1)

        point   = xyz[self.atom1,:]
        vector  = xyz[self.atom2,:] - point
        vector /= np.linalg.norm(vector)
        M       = rotation_matrix_arbitrary_axis(point, vector, angle)

        rotvec = np.ones([4,1], dtype=xyz.dtype)
        for atom_ndx in self.rotation_atoms:
            rotvec[:3] = xyz[atom_ndx,:].reshape([3,1])
            xyz[atom_ndx, :] = np.dot(M, rotvec)[:3,0]

        # rotate all symmetry partners
        for i, partner in enumerate(self.symmetry_partners):
            partner.move(xyz, rand=rand, angle=angle)

 

        




class StateReporter(object):
    """
    Generate Reports of the current simulation
    """

    def __init__(self, properties, fieldwidths, outfile=None, delimiter=' '):

        if outfile is not None:
            # backup old outfile if any
            if os.path.isfile(outfile):
                try:
                    maxn = max([int(s.strip('#').split('.')[-1]) for s in glob.glob("#" + outfile + ".*#")])
                except ValueError:
                    maxn = 0
                os.rename(outfile, "#" + outfile + ".{:d}#".format(maxn+1))
 
            # open outfile
            self.outstream = open(outfile, 'w')
        else:
            self.outstream = sys.stdout

        self.properties  = properties
        self.fieldwidths = fieldwidths
        self.delimiter   = delimiter

        # write header
        self.outstream.write('#')
        for i, prop in enumerate(self.properties):
            if i > 0:
                self.outstream.write(self.delimiter)
            self.outstream.write("{{:>{}s}}".format(fieldwidths[i]).format(prop))
        self.outstream.write('\n')


    def report(self, data):
        """
        Report data given as keyword arguments
        data = {property: value, ...} dictionary
        """

        for i, prop in enumerate(self.properties):
            if i > 0:
                self.outstream.write(self.delimiter)
            else:
                self.outstream.write(' ')

            value = data[prop]
            if type(value) == type(""):
                s = "{{:>{}s}}".format(self.fieldwidths[i]).format(value)

            elif type(value) == type(int()):
                s = "{{:>{}d}}".format(self.fieldwidths[i]).format(value)
            else:
                s = "{{:>{}.4e}}".format(self.fieldwidths[i]).format(value) 
            self.outstream.write(s)
        self.outstream.write('\n')
                


    def __del__(self):
        # if outstream is not stderr or stdout, close it
        if self.outstream.fileno() != sys.stdout.fileno() and self.outstream.fileno() != sys.stderr.fileno():
            self.outstream.close()








class MDrecorder(object):
    """
    Record and parse OpenMM output and make it available to other functions
    """
    def __init__(self, delimiter=','):
        self.categories = []
        self.data       = None
        self.length     = 0
        self.delimiter  = delimiter

    def __str__(self):
        stringrep  = ''
        stringrep += '#'
        for c, category in enumerate(self.categories):
            if c > 0:
                stringrep += ' ' + self.delimiter
            stringrep += category
        stringrep += '\n'
        for i in range(self.length):
            for c, category in enumerate(self.categories):
                if c > 0:
                    stringrep += ' ' + self.delimiter
                stringrep += str(self.data[c][i])
            stringrep += '\n'
        return stringrep

    def __repr__(self):
        return self.__str__()

    def write(self, s):
        """
        Insert data into this recorder
        """
        if s[0] == '#':
            self.categories = s.split('\t')
            self.data       = len(self.categories) * [[]]
            self.length     = 0
        elif s[0] != '\n':
            numbers = [float(n) for n in s.split('\t')]
            for i, n in enumerate(numbers):
                self.data[i].append(n)
            self.length += 1
    
    def save(self, filename):
        """
        Write contents of this recorder to a file
        """
        with open(filename, 'w') as of:
            of.write(self.__str__())
