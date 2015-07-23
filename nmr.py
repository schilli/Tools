# -*- coding: UTF-8 -*-
# Compute NMR observables from MD trajectories

from __future__ import print_function

import sys, os, time, string
import mdtraj
import numpy as np
import cPickle as pickle
import zipfile, bz2
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpi4py import MPI


# ============================================================================ #


class corrFunction(object):
    """
    correlation function data, additional information and manipulation methods

    Parameters
    ----------
    corr : (nvec, nframes) array
        Correlation functions

    std : (nvec, nframes) array
        Correlation function standard deviations

    error : (nvec, nframes) array
        Correlation function standard error of the mean

    info : dict
        Dictionary with information on correlation functions
    """

    def __init__(self, corr=None, std=None, error=None, info=None):
        self.corr        = corr
        self.std         = std
        self.error       = error
        self.resid       = info['bondvecinfo']['resid'     ]
        self.resindex    = info['bondvecinfo']['resindex'  ]
        self.resnames    = info['bondvecinfo']['resnames'  ]
        self.atomindex   = info['bondvecinfo']['atomindex' ]
        self.atomnames   = info['bondvecinfo']['atomnames' ]
        self.element     = info['bondvecinfo']['element'   ]
        self.chain       = info['bondvecinfo']['chain'     ]
        self.bondlength  = info['bondvecinfo']['bondlength']
        self.bondvec     = info['bondvecinfo']['bondvec'   ]
        self.fitgroup    = info['bondvecinfo']['fitgroup'  ]
        self.fit         = info['bondvecinfo']['fit'       ]
        self.dt          = info['bondvecinfo']['dt'        ]
        self.topfilename = info['topfilename']
        self.npzfilename = info['npzfilename']
        self.trjfilename = info['trjfilename'] 
        self.frames      = info['frames'     ]


# ============================================================================ #


class orderParameter(object):
    """
    Bond vector order parameters.

    Parameters
    ----------
    corrfilenames : list
        List with *.zip filenames containing correlation functions.

    converged : boolean, optional, default: True
        Use only converged correlation functions for averaging

    verbose : boolean, optional
        Report progress and other verbose output.

    **kwargs : optional keyword arguments
        All remaining keyword arguments are passed to the order parameter method.
    """

# ==================================== #

    def __init__(self, corrfilenames=None, converged=True, verbose=True, **kwargs):
        self.corrfilenames
        self.converged     = converged
        self.verbose       = verbose

        self.corrlist      = []

        self.method_kwargs = kwargs


        # for plotting
        self.figure  = None
        self.axs     = None
        self.cids    = None # connection ids for plot event handlers
        self.lines   = []   # line object handels
        self.corrset = 0    # which set of correlation functions from the corrlist should be plottet
        self.corridx = 0    # which correlatin function to plot next
        self.leftButtonPressed = False
 
# ==================================== #

    def load(self, corrfilenames=None):
        """
        Load sets of correlation functions.

        Parameters
        ----------
        corrfilenames : list
            List with *.zip filenames containing correlation functions.
 
        """

        if corrfilenames is not None:
            self.corrfilenames = corrfilenames

        starttime = time.time()

        # load correlation functions
        self.corrlist = []
        print("Loading {} set{} of correlation functions:".format(len(self.corrfilenames), *['s' if i>1 else '' for i in [len(self.corrfilenames)]]))
        for nf, filename in enumerate(self.corrfilenames):
            corr, corrstd, corrstdmean, info = load_corr(filename)
            self.corrlist.append(corrFunction(corr, corrstd, corrstdmean, info))
            if self.verbose:
                print("\rProgress: {:3.0f}%".format(100.0*nf/len(self.corrfilenames)), end="")
                sys.stdout.flush()

        # report runtime
        if self.verbose:
            print("\rLoading took: {:.2f} sec.".format(time.time() - starttime))

# ==================================== #

    def estimate(self, method="mean", converged=True, **kwargs):
        """
        Estimate bond vector order parameters from correlation functions.

        Parameters
        ----------
        method: string, optional
            The method to use for order parameter computation.
            Options are:
                "mean"          Use the mean of the final quarter as order parameter
                "single exp"    Fit correlation functions to single exponential
                "double exp"    Fit correlation functions to double exponential
                "extLS"         Use the extended least squares method (method 3) from:
                                JPCB 2008, 112, 6155-6158, pubs.acs.org/doi/abs/10.1021/jp077018h
                "iRED"          Use the iRED method from:
                                Brüschweiler, JACS 2002, 124, 4522-4534, pubs.acs.org/doi/abs/10.1021/ja012750u

        converged : boolean, optional, default: True
            Use only converged correlation functions for averaging
        """
 
        self.method    = method
        self.converged = converged

        # select order parameter estimation method
        if self.method == "mean":
            self.estimate_mean(self.corrlist, kwargs)
        else:
            print("Order parameter estimation method unknown: {}".format(self.method))
            sys.exit(1)
 
# ==================================== #

    def estimate_mean(self, converged=True, diffThreshold=0.02, stdThreshold=0.02):
        """
        Estimate bond vector order parameters as the mean of the last quarter of the correlation function.

        Parameters
        ----------
        converged : boolean, optional
            If True, only use converged correlation functions.
        diffThreshold : float, optional
            If the quarters 3 and 4 of the correlation function differ more than this threshold,
            they are not considered converged
        stdThreshold : float, optional
            If the mean standard deviation of quarter 3 and 4 is larger than this threshold,
            they are not considered converged.
        """

        for corrfun in self.corrlist:
            length = corrfun.corr.shape[1]

#        length            = corr.shape[1]
#        quarter           = length/4
#        thirdQuarter      = corr[:,2*quarter:3*quarter]
#        fourthQuarter     = corr[:,3*quarter:4*quarter]
#        fourthQuarterMean = fourthQuarter.mean(1)
#
#        difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
#        stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
#        convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)
#
#        return fourthQuarterMean, convergence 
 

# ==================================== #

    def plot_corr(self, event=None, corrset=0):
        """
        Plot correlation functions.
        
        Parameters
        ----------
        event : matplotlib mouse event
            Only there to update the plot on various matplotlib events, like mouse button presses
        corrset: integer, optional
            Which set of correlation functions to use for plotting
        """

        if self.figure is None:
            # create figure and axis
            self.figure = plt.figure()
            self.axs    = self.figure.add_subplot(111)  
            self.cids['button_press'  ] = self.fig.canvas.mpl_connect('button_press_event',   self._onclick)
            self.cids['button_release'] = self.fig.canvas.mpl_connect('button_release_event', self._onrelease)
            self.cids['motion'        ] = self.fig.canvas.mpl_connect('motion_notify_event',  self._onmove)
            self.cids['scroll'        ] = self.fig.canvas.mpl_connect('scroll_event',         self._onscroll)
            self.cids['close'         ] = self.fig.canvas.mpl_connect('close_event',          self._onclose)

            self.corrset = corrset
            self.corridx = 0
            self.lines   = []


#        self.corr        = corr
#        self.std         = std
#        self.error       = error
#        self.resid       = info['bondvecinfo']['resid'     ]
#        self.resindex    = info['bondvecinfo']['resindex'  ]
#        self.resnames    = info['bondvecinfo']['resnames'  ]
#        self.atomindex   = info['bondvecinfo']['atomindex' ]
#        self.atomnames   = info['bondvecinfo']['atomnames' ]
#        self.element     = info['bondvecinfo']['element'   ]
#        self.chain       = info['bondvecinfo']['chain'     ]
#        self.bondlength  = info['bondvecinfo']['bondlength']
#        self.bondvec     = info['bondvecinfo']['bondvec'   ]
#        self.fitgroup    = info['bondvecinfo']['fitgroup'  ]
#        self.fit         = info['bondvecinfo']['fit'       ]
#        self.dt          = info['bondvecinfo']['dt'        ]
#        self.topfilename = info['topfilename']
#        self.npzfilename = info['npzfilename']
#        self.trjfilename = info['trjfilename'] 
#        self.frames      = info['frames'     ]
 

        # remove old data
        while len(self.lines) > 0:
            line = self.lines.pop()
            line.remove() 

        # plot data
        corrFun     = self.corrlist[self.corrset]
        xdata       = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, corrFun.corr.shape[1])
        self.lines += self.axs.plot(self.xdata, corrFun.corr[self.corridx,:], 'b')

        # set axis limits
        self.axs.set_ylim(min(0, corrFun.corr.min()), max(1, corrFun.corr.max()))

        # annotate plot
        self.axs.set_title("{} {}".format(corrFun.resnames[self.corridx], corrFun.resids[self.corridx]))
        self.axs.set_ylabel("correlation")
        self.axs.set_xlabel("time [ps]")
        self.legend(loc="lower left")

        self.fig.canvas.draw()
        plt.show()
 
# ==================================== #

    def _onclick(self, event):
        if event.button == 1:
            self.leftButtonPressed = True
            try:
                xmin, xmax = self.axs.get_xlim()
                ncorr = self.corrlist[self.corrset].corr.shape[0]
                self.corridx = int(np.round((ncorr-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.corridx < 0:
                    self.corridx = 0
                elif self.corridx >= ncorr:
                    self.corridx = ncorr - 1
                self.plot()
            except:
                pass 

# ==================================== #
    
    def _onclose(self, event):
        self.figure  = None
        self.axs     = None
        self.cids    = {}
        self.lines   = []

# ============================================================================ #
 

def _fit_trj(trj, fitframe=0, ref=None, fitgroup="name CA", parallel=True):
    """
    Fit MD trajectory to reference based on the specified set of atoms.
    Fitting is done in place.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to fit.

    fitframe : int, optional
        Frame index of reference frame.

    ref : mdtraj trajectory
        Reference MDTRAJ trajectory.
        If none, the first trajectory is used as a reference.

    fitgroup : string or integer array
        mdtraj selection string or integer array of atom indices
        to select a set of atoms for RMSD fit.

    parallel : boolean, optional
        Toggle parallel fitting algorithm of mdtraj
    """
    if type(ref) == type(None):
        ref = trj

    try:
        fit_atomndx = trj.topology.select(fitgroup)
    except:
        fit_atomndx = fitgroup

    trj.superpose(ref, fitframe, atom_indices=fit_atomndx, parallel=parallel)
 

# ============================================================================ #

def _vec_atoms(trj, bondvec=None):
    """
    Compute bond atom indices.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory for topology information.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    Returns
    -------
    bondvec_ndx : (nvec, 2) array
        Integer array of all atom indices for bond vector atoms
    """

    # default option
    if type(bondvec) == type(None):
        selectiontext1   = "name N and not resname PRO and resid >= 1"
        selectiontext2   = "name H and not resname PRO and resid >= 1"
        atom1_ndx        = trj.topology.select(selectiontext1)
        atom2_ndx        = trj.topology.select(selectiontext2) 

    else:
        sel1 = bondvec[0]
        sel2 = bondvec[1]
        # if bond vectors are already specified with integer arrays
        try:
            bondvec_ndx = np.zeros([len(sel1), 2], dtype=np.int)
            bondvec_ndx[:,0] = np.array(sel1)
            bondvec_ndx[:,1] = np.array(sel2)
            return bondvec_ndx
        # if bond vectors are specified with selection texts
        except:
            atom1_ndx        = trj.topology.select(sel1)
            atom2_ndx        = trj.topology.select(sel2)  

    bondvec_ndx      = np.zeros([atom1_ndx.shape[0], 2], dtype=np.int)
    bondvec_ndx[:,0] = atom1_ndx
    bondvec_ndx[:,1] = atom2_ndx

    return bondvec_ndx


# ============================================================================ #


def _bond_vec(trj, bondvec_ndx):
    """
    Compute bond vectors from trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory

    bondvec_ndx : (nvec, 2) array
        Bond vector indices

    Returns
    -------
    bondvec : (nframes, nvec, 3) array
        Array of normalized bond vectors
    info : dict
        Dictionary with information on the bondvectors.
        Keys are:
        'dt'            time step in pico seconds
        'bondlength'    one bondlength per bond vector
        'resnames'      two sequences of resnames, one for each atom in bondvector
        'resids'                         resids
        'resindex'                       resindices
        'atomnames'                      atomnames
        'atomindex'                      atomindices
        'element'                        elements
        'chain'                          chains ids
    """

    # extract bond vector subtrajectories
    atom1_trj = trj.atom_slice(bondvec_ndx[:,0])
    atom2_trj = trj.atom_slice(bondvec_ndx[:,1])

    # compute bond vectors
    bondvec = atom1_trj.xyz - atom2_trj.xyz

    # normalize bond vectors
    bondlength = (((bondvec**2).sum(2))**0.5)
    bondvec /= bondlength.reshape([bondvec.shape[0], bondvec.shape[1], 1])

    info = {}
    info['dt'        ] = trj.timestep
    info['bondlength'] = bondlength
    info['resnames'  ] = [[], []]
    info['resid'     ] = [[], []]
    info['resindex'  ] = [[], []]
    info['atomnames' ] = [[], []]
    info['atomindex' ] = [[], []]
    info['element'   ] = [[], []]
    info['chain'     ] = [[], []]

    # get info on atoms and residues
    for atom1, atom2 in zip(atom1_trj.top.atoms, atom2_trj.top.atoms):
        info['resnames'  ][0].append(atom1.residue.name)
        info['resnames'  ][1].append(atom2.residue.name)
        info['resid'     ][0].append(atom1.residue.resSeq)
        info['resid'     ][1].append(atom2.residue.resSeq)
        info['resindex'  ][0].append(atom1.residue.index)
        info['resindex'  ][1].append(atom2.residue.index)
        info['atomnames' ][0].append(atom1.name)
        info['atomnames' ][1].append(atom2.name)
        info['atomindex' ][0].append(atom1.index) 
        info['atomindex' ][1].append(atom2.index) 
        info['element'   ][0].append(atom1.element) 
        info['element'   ][1].append(atom2.element)  
        info['chain'     ][0].append(atom1.residue.chain.index) 
        info['chain'     ][1].append(atom2.residue.chain.index)   

    return bondvec, info
 

# ============================================================================ #

 
def _angular_correlation_function(bondvec, verbose=True):
    """
    Compute 2nd order Legandre polynomial correlation functions of bond vectors.

    Parameters
    ----------
    bondvec : (nframes, nvectors, ndim) array
        normalized bond vectors
    verbose : boolean, optional
        Print progress messages or not

    Returns
    -------
    corr : (nvectors, nframes/2) array
        correlation functions
    corr_std : (nvectors, nframes/2) array
        standard deviations of correlation values
    corr_stdmean : (nvectors, nframes/2) array
        standard error of the mean
    """

    # compute bond vector correlation functions
    nframes, nvectors, ndim  = bondvec.shape
    corrnframes   = int(round(0.5*nframes))
    corr          = np.zeros([nvectors, corrnframes]) # correlation functions
    corr_std      = np.zeros([nvectors, corrnframes]) # standard deviation
    corr_stdmean  = np.zeros([nvectors, corrnframes]) # standard error of the mean
 
    for vector in range(nvectors):

        if verbose:
            print("\rProgress: {:3.0f}%".format(100.0*(vector+1)/nvectors), end="")
            sys.stdout.flush()

        # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
        P2 = np.polynomial.legendre.legval(np.dot(bondvec[:,vector,:], bondvec[:,vector,:].T), [0,0,1])

        # compute the correlation function for each lag time
        for frame in range(corrnframes):
            d = np.diagonal(P2, frame)
            corr        [vector, frame] = d.mean()
            corr_std    [vector, frame] = d.std()
            corr_stdmean[vector, frame] = corr_std[vector, frame] / d.shape[0]**0.5 

    if verbose:
        print()

    return corr, corr_std, corr_stdmean

 
# ============================================================================ #


def save_corr(savefilename, corr, corrstd, corrstdmean, bondvecinfo, topfilename, trjfilename, frames):
    """
    Save correlation functions to disk including all additional information.

    Parameters
    ----------
    savefilename : string
        Filename of the resulting zipfile containing all the data

    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    bondvecinfo : dict
        Dictionary with information on bondlength, resids and resnames

    trjfilename : string
        Filename of the trajectory the original MD data was taken from.
        Stored as is, so better provide a complete path.

    topfilename : string
        Filename of the topology the original MD data was taken from.
        Stored as is, so better provide a complete path.

    frames : (startframe, endframe) tuple
        Range of frames used for correlation function computation
    """

    tmpdirnamelen     = 8
    tmpdirname        = "tmpdir_" + ''.join(np.random.choice([i for i in string.letters + string.digits], size=tmpdirnamelen))
    os.mkdir(tmpdirname)
    savefilepath      = os.path.dirname(savefilename)
    npzfilename       = "corr.npz"
    infofilename      = "info.dat"

    if not os.path.isdir(savefilepath):
        os.mkdir(savefilepath)

    info = {}
    info['npzfilename'] = npzfilename
    info['bondvecinfo'] = bondvecinfo
    info['trjfilename'] = trjfilename
    info['topfilename'] = topfilename
    info['frames']      = frames
 
    # save arrays
    with open(tmpdirname + '/' + npzfilename, 'w') as outfile:
        #np.savez(outfile, corr, corrstd, corrstdmean)
        np.savez_compressed(outfile, corr=corr, corrstd=corrstd, corrstdmean=corrstdmean)

    # save info
    with open(tmpdirname + '/' + infofilename, 'w') as outfile:
        outfile.write(bz2.compress(pickle.dumps(info)))

    # pack both into zipfile
    with zipfile.ZipFile(savefilename, 'w') as outfile:
        outfile.write(tmpdirname + '/' + infofilename, arcname=infofilename)
        outfile.write(tmpdirname + '/' +  npzfilename, arcname=npzfilename )

    # remove npz and pickle files
    os.remove(tmpdirname + '/' + infofilename)
    os.remove(tmpdirname + '/' +  npzfilename)
    os.rmdir (tmpdirname)


# ============================================================================ #


def load_corr(filename):
    """
    Load correlation functions including all additional information.

    Parameters
    ----------
    filename : string
        Filename of the zipfile containing all the data

    Returns
    -------
    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    info
        Dictionary with information on the loaded correlation functions
    """

    # extract files
    with zipfile.ZipFile(filename, 'r') as infile:
        zipfilenames = infile.namelist()
        for zipfilename in zipfilenames:
            infile.extract(zipfilename)

    for zipfilename in zipfilenames:
        with open(zipfilename, 'r') as infile:
            # load info dictionary
            if zipfilename.find('.dat') >= 0:
                info = pickle.loads(bz2.decompress(infile.read()))

            # load corr arrays
            if zipfilename.find('.npz') >= 0:
                arrays      = np.load(infile)
                corr        = arrays['corr']
                corrstd     = arrays['corrstd']
                corrstdmean = arrays['corrstdmean']

        # remove extracted files
        os.remove(zipfilename)

    return corr, corrstd, corrstdmean, info


# ============================================================================ #
 

def bondvec_corr(trj, bondvec=None, fitgroup=None, parallel=True, saveinfo=None, verbose=True):
    """
    Compute bond vector correlation functions from mdtraj trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to compute the order parameters of.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.

    parallel : boolean, optional
        Do the fitting in parallel on as many processors as available on the machine.

    saveinfo : dict, optional
        Information needed to save the correlation functions.
        Required keys are:
            'topfilename'   the topology   file the data comes from
            'trjfilename'   the trajectory file the data comes from
            'frames'        tuple with initial and final frame
            'zipfilename'   name of the zipfile to store correlation functions in

    verbose : boolean, optional
        Toggle verbosity level.


    Returns
    -------
    corr : (nvec, t) array
        Bond vector correlation functions.

    bondvecinfo :
        Dictionary with information on the bond vectors.
    """

    # fit trajectory
    if type(fitgroup) != type(None):
        _fit_trj(trj, fitgroup=fitgroup, parallel=parallel)

    bondvec_ndx                = _vec_atoms(trj, bondvec)
    bondvectors, bondvecinfo   = _bond_vec(trj, bondvec_ndx)
    corr, corrstd, corrstdmean = _angular_correlation_function(bondvectors, verbose=verbose)

    # store additional bond vector information
    if type(fitgroup) != type(None):
        bondvecinfo['fit' ] = True
    bondvecinfo['fitgroup'] = fitgroup
    bondvecinfo['bondvec' ] = bondvec

    if type(saveinfo) == type(dict()):
        save_corr(saveinfo['zipfilename'], corr, corrstd, corrstdmean, bondvecinfo, saveinfo['topfilename'], saveinfo['trjfilename'], saveinfo['frames'])

    return corr, bondvecinfo

    
# ============================================================================ #


def bondvec_corr_batch_mpi(topfilename, trjfilenames, savepath, subtrjlength=None, bondvec=None, fitgroup=None):
    """
    Compute bond vector correlation functions from mdtraj trajectory.
    MPI parallelized batch version that stores the results on disk.

    Parameters
    ----------
    topfilename : string
        topology filename, e.g. PDB

    trjfilenames : sequence of strings
        trajectory filenames, e.g. XTC

    savepath : string
        path to store the results at

    subtrjlength: float, optional
        Subtrajetory length in ps.
        If none, whole trajectories are used.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.
    """

    # set up some timer and counter
    tc = {}
    tc['runtimer']  = time.time()
    tc['loadtimer'] = 0.0
    tc['corrtimer'] = 0.0
    tc['nsubtrjs' ] = 0
    tc['nframes'  ] = 0
 

    # set up MPI
    comm   = MPI.COMM_WORLD
    nranks = comm.Get_size()
    myrank = comm.Get_rank()
    root   = 0

    # to fit or not to fit?
    fit = "nofit"
    if type(fitgroup) != type(None):
        fit = "fit"

    # distribute tasks over ranks
    if myrank == root:
        ntasks        = len(trjfilenames)
        ntasksperrank = np.array(nranks * [ntasks / nranks], dtype=np.int)
        ntasksperrank[:ntasks%nranks] += 1

        for rank in range(0, nranks):
            taskstaken = sum(ntasksperrank[:rank])
            task = {}
            task['topfilename']  = topfilename
            task['trjindices']   = range(rank, len(trjfilenames), nranks)
            task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
            task['savepath']     = savepath

            comm.send(task, dest=rank, tag=rank)

    task = comm.recv(source=root, tag=myrank)
#    print ("rank {}: ".format(myrank), task['trjindices'], task['trjfilenames'])
#    sys.exit(0)

    # do the assigned piece of work
    for nf, trjfilename in enumerate(task['trjfilenames']):
        # determinde dt and chunksize
        trjs      = mdtraj.iterload(trjfilename, top=task['topfilename'], chunk=2)
        trj       = trjs.next()
        dt        = trj.timestep
        chunksize = int(subtrjlength / dt)
        trjindex  = task['trjindices'][nf]

        if myrank == root:
            print("\r", 100*" ", end="")
            print("\rProgress rank root: {:3.0f}% ".format(100.0*nf/len(task['trjfilenames']), nf), end="")
            sys.stdout.flush()
 
        loadstarttime = time.time()
        for ntrj, trj in enumerate(mdtraj.iterload(trjfilename, top=task['topfilename'], chunk=chunksize)):
            tc['loadtimer'] += time.time() - loadstarttime
            tc['nsubtrjs' ] += 1

            if trj.n_frames == chunksize:
                tc['nframes'] += trj.n_frames

                saveinfo = {}
                saveinfo['topfilename'] = task['topfilename']
                saveinfo['trjfilename'] = trjfilename
                saveinfo['frames'     ] = (ntrj * chunksize, (ntrj+1)*chunksize-1)
                saveinfo['zipfilename'] = savepath + '/' + ''.join(os.path.basename(trjfilename).split('.')[:-1]) + '_{}.zip'.format(fit)
                saveinfo['zipfilename'] = "{}/{}_no{}_frames{}to{}_{}.zip".format(savepath,
                                                                                  ''.join(os.path.basename(trjfilename).split('.')[:-1]),
                                                                                  trjindex,
                                                                                  saveinfo['frames'][0],
                                                                                  saveinfo['frames'][1],
                                                                                  fit)
            
                corrstarttime = time.time()
                bondvec_corr(trj, bondvec=bondvec, fitgroup=fitgroup, parallel=False, saveinfo=saveinfo, verbose=False)
                tc['corrtimer'] += time.time() - corrstarttime
            loadstarttime = time.time()

            if myrank == root:
                print(".", end="")
                sys.stdout.flush()

    print("")


    # report timers and counters to root
    if myrank == root:
        for rank in range(1, nranks):
            tci = comm.recv(source=rank, tag=rank)
            tc['loadtimer'] += tci['loadtimer'] 
            tc['corrtimer'] += tci['corrtimer'] 
            tc['nsubtrjs' ] += tci['nsubtrjs' ] 
            tc['nframes'  ] += tci['nframes'  ] 
    else:
        comm.send(tc, dest=root, tag=myrank)

    comm.barrier()

    tc['runtimer'] = time.time() - tc['runtimer']
    if myrank == root:
        print("Summary")
        print("-------\n")
        print("Loaded {} subtrajectories with a total of {} frames".format(tc['nsubtrjs'], tc['nframes']))
        print("Total runtime:                        {:8.0f} sec.".format(tc['runtimer' ]))
        print("Aggregated time for loading:          {:8.0f} sec.".format(tc['loadtimer']))
        print("Aggregated time for corr computation: {:8.0f} sec.".format(tc['corrtimer']))
        

# ============================================================================ #


def _check_corr_convergence(corr, diffThreshold=0.02, stdThreshold=0.02):
    """
    Check the convergence of the bond vector correlation functions

    Parameters
    ----------
    corr : (nbonds, nframes) array
        Correlation functions

    diffThreshold: float, optional
        Maximum mean difference for convergence check.

    stdThreshold: float, optional
        Maximum stdev difference for convergence check.

    Returns
    -------
    

        (convergence values, boolean array of converged correlation functions)
    """
    length            = corr.shape[1]
    quarter           = length/4
    thirdQuarter      = corr[:,2*quarter:3*quarter]
    fourthQuarter     = corr[:,3*quarter:4*quarter]
    fourthQuarterMean = fourthQuarter.mean(1)
    difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
    stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
    convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)
    return fourthQuarterMean, convergence
 


# ============================================================================ #


def order_parameter_mean(corr, converged=True, diffThreshold=0.02, stdThreshold=0.02):

    length            = corr.shape[1]
    quarter           = length/4
    thirdQuarter      = corr[:,2*quarter:3*quarter]
    fourthQuarter     = corr[:,3*quarter:4*quarter]
    fourthQuarterMean = fourthQuarter.mean(1)

    difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
    stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
    convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)

    return fourthQuarterMean, convergence 


# ============================================================================ #


def order_parameter(corrfilenames, method="mean", converged=True, verbose=True, **kwargs):
    """
    Compute bond vector order parameters from correlation functions.

    Parameters
    ----------
    corrfilenames : list
        List with *.zip filenames containing correlation functions.
        
    method: string, optional
        The method to use for order parameter computation.
        Options are:
            "mean"          Use the mean of the final quarter as order parameter
            "single exp"    Fit correlation functions to single exponential
            "double exp"    Fit correlation functions to double exponential
            "extLS"         Use the extended least squares method (method 3) from:
                            JPCB 2008, 112, 6155-6158, pubs.acs.org/doi/abs/10.1021/jp077018h
            "iRED"          Use the iRED method from:
                            Brüschweiler, JACS 2002, 124, 4522-4534, pubs.acs.org/doi/abs/10.1021/ja012750u

    converged : boolean, optional, default: True
        Use only converged correlation functions for averaging

    verbose : boolean, optional
        Report progress and other verbose output.

    **kwargs : optional keyword arguments
        All remaining keyword arguments are passed to the order parameter method.
    """

    starttime = time.time()

    # load correlation functions
    corrlist = []
    print("Loading {} set{} of correlation functions:".format(len(corrfilenames), *['s' if i>1 else '' for i in [len(corrlist)]]))
    for nf, filename in enumerate(corrfilenames):
        corr, corrstd, corrstdmean, info = load_corr(filename)
        corrlist.append(corrFunction(corr, corrstd, corrstdmean, info))
        if verbose:
            print("\rProgress: {:3.0f}%".format(100.0*nf/len(corrfilenames)), end="")
            sys.stdout.flush()

    # select order parameter estimation method
    if method == "mean":
        S2 = order_parameter_mean(corrlist, **kwargs)
    else:
        print("Order parameter estimation method unknown: {}".format(method))
        sys.exit(1)

    # report runtime
    print("\rRuntime: {:.2f} sec.".format(time.time() - starttime))

    return S2
 

    
