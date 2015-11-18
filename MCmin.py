from __future__ import print_function

import sys, time
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

    def __init__(self, xyz=None, topology=None, time=None, unitcell_lengths=None, unitcell_angles=None, other=None, pdbfile=""):

        # if created like the parent class
        if xyz is not None and topology is not None:
            super(MCtrj, self).__init__(xyz, topology, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

        # if created from parent class instance
        elif other is not None:
            super(MCtrj, self).__init__(other.xyz, other.top, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

        self.pdbfile                   = pdbfile
        self.openmm_platform_is_set_up = False     # keep track of what has been initialised
        self.openmm_is_set_up          = False
        self.MC_is_set_up              = False     
        self.MC_attempt_counter        = 0         # counts the number of attempted MC moves
        self.MC_step_counter           = 0         # counts the number of accepted MC
        self.resids                    = np.array([r.resSeq for r in self.top.residues], dtype=np.int)

        # data storage
        self.Epot      = np.zeros(self.n_frames, dtype=np.float)
        self.SAXS_chi2 = np.zeros(self.n_frames, dtype=np.float)

        # list of attributes added elsewhere
        self.philist = []   # list of residue indices whose angles to rotate
        self.psilist = []
        self.maxrot  = 0.0  # maximum rotation angle with unit
        self.nrot    = 0    # number of angles to rotate each step

        # constant with units
        self.R = const.R /1e3 * unit.kilojoule/unit.kelvin/unit.mole  # gas constant
 

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

        self.openmm_emtol    = emtol
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
        self.openmm_Epot_recorder = MDrecorder()
        self.openmm_simulation.reporters.append(app.StateDataReporter(self.openmm_Epot_recorder, 1, potentialEnergy=True, separator='\t')) 

        self.openmm_is_set_up = True



    def MC_setup(self, philist=None, psilist=None, maxrot=120*unit.degree, nrot=1, seed=None, **kwargs):
        """
        Set up Monte Carlo parameters

        Parameters:
        -----------
        philist:    residue IDs of phi angles to rotate (numbering as in pdb file)
        psilist:    residue IDs of psi angles to rotate (numbering as in pdb file)
        maxrot:     maximum rotation angle (with unit)
        nrot (int): number of rotations per step
        seed (int): random number generator seed
        kwargs:     Additional keyword arguments are passed on to openmm_setup()
        """

        self.philist = philist
        self.psilist = psilist
        self.maxrot  = maxrot
        self.nrot    = nrot

        if type(seed) == type(int()):
            np.random.seed(seed)

        if self.philist is None:
            self.philist = range(1, self.n_residues-1)
        else:
            self.philist.sort()
        if self.psilist is None:
            self.psilist = range(0, self.n_residues-2)
        else:
            self.psilist.sort()
        self.philist = np.array(self.philist, dtype=np.int)
        self.psilist = np.array(self.psilist, dtype=np.int)

        self.phi_atom_indices = md.compute_phi(self)[0][self.philist-1,:]
        self.psi_atom_indices = md.compute_psi(self)[0][self.psilist  ,:]

        self.MC_attempt_counter = 0
        self.MC_step_counter    = 0
        self.MC_is_set_up       = True

        # also make sure OpenMM is set up
        if not self.openmm_is_set_up:
            self.openmm_setup(**kwargs)

        # minimize the energy of the initial structure
        self.minimize_energy()


    def MC_run(self, nsteps=None, nattempts=None, maxrot=120*unit.degree, nrot=1, **kwargs):
        """
        Run a monte carlo simulation

        Parameters
        ----------
        nsteps (int):                        number of accepted monte carlo steps to perform
        nattempts (int):                     number of monte calro attempts to try (overrides nsteps)
        maxrot (unit.degree or unit.radian): standard deviation of rotation angles
        nrot (int):                          how many angle rotation should be attempted per trial move
        kwargs:                              keyword arguments passed on to MC_step()
        """
        # allocate enough memory
        framesNeeded = (nsteps + self.MC_step_counter+1) - self.n_frames
        if framesNeeded > 0:
            self.allocate_memory(framesNeeded)
        
        # try this many MC attempts
        if type(nattempts) == type(int()):
            for n in range(nattempts):
                self.MC_step(maxrot=maxrot, nrot=nrot, **kwargs)
        # or perform this many successfull steps
        elif type(nsteps) == type(int()):
            while self.MC_step_counter < nsteps:
                self.MC_step(maxrot=maxrot, nrot=nrot, **kwargs)

        # abort
        else:
            print("Either nsteps or nattempts needs to be an integer!")
            sys.exit(1)




    def MC_step(self, maxrot=120*unit.degree, nrot=1, verbose=False, **kwargs):
        """
        Perform a single monte carlo step

        Parameters
        ----------
        maxrot (unit.degree or unit.radian): standard deviation of rotation angles
        nrot (int):                          how many angle rotation should be attempted per trial move 
        verbose (bool):                      print status messages
        kwargs:                              keyword arguments passed on to MC_setup()
        """
        if not self.MC_is_set_up:
            self.MC_setup(maxrot=maxrot, nrot=nrot, **kwargs)
        else:
            self.nrot   = nrot
            self.maxrot = maxrot

        self.MC_attempt_counter += 1
        self.MC_step_counter    += 1

        # if we don't have any space left, allocate more memory (double it)
        if self.MC_step_counter >= self.n_frames:
            self.allocate_memory()

        # get rotation angles
        phi_indices, phi_angles, psi_indices, psi_angles = self.get_rotation_angles()

        # rotate angles
        self.rotate(phi_indices, phi_angles, psi_indices, psi_angles)

        # minimize energy
        self.minimize_energy()

        # compute SAXS fit

        # metropolis criterion
        deltaE = (self.Epot[self.MC_step_counter] - self.Epot[self.MC_step_counter-1]) * unit.kilojoule / unit.mole
        beta   = self.openmm_T * self.R
#        print(deltaE)
#        print(self.R)
#        print(beta)
#        print(-deltaE/beta)
#        sys.exit()
        acceptance_probability = 1.0
        if deltaE > 0.0 * unit.kilojoule / unit.mole:
            acceptance_probability = np.exp(- deltaE / beta)
        if verbose:
            print("Step: {:5d}, attempt: {:5d}, dE: {:10.2f} {:s}, ap: {:4.2f}, ".format(self.MC_step_counter,
                self.MC_attempt_counter, deltaE.value_in_unit(deltaE.unit), deltaE.unit, acceptance_probability), end='') 
        if np.random.rand() > acceptance_probability:
            # reject
            if verbose:
                print("rejected")
            self.MC_step_counter -= 1
        else:
            if verbose:
                print("accepted")


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
 
        # minimize energy
        self.openmm_simulation.minimizeEnergy(tolerance=self.openmm_emtol) 
        state = self.openmm_simulation.context.getState(getPositions=True, getEnergy=True)
 
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
        phi_angles = self.maxrot * np.random.randn(len(phi_indices))
        psi_angles = self.maxrot * np.random.randn(len(psi_indices))

        # add up angles that have been selected multiple times
        if len(phi_indices) > 0:
            phi_indices_dict = {p: 0.0*self.maxrot.unit for p in set(phi_indices)}
            for ndx, phi_ndx in enumerate(phi_indices):
                phi_indices_dict[phi_ndx] += phi_angles[ndx]
            phi_indices, phi_angles = zip(*phi_indices_dict.iteritems())
        if len(psi_indices) > 0:
            psi_indices_dict = {p: 0.0*self.maxrot.unit for p in set(psi_indices)}
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
