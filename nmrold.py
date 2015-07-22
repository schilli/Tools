# Compute NMR observables from MD trajectories

from __future__ import print_function

import sys, os
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpi4py import MPI


# ============================================================================ #


def _fit_trajectory(trajectory, fitframe=0, reference=None, selectiontext="name CA", parallel=True):
    """
    Fit MD trajectory to reference based on the specified set of atoms.
    Fitting is done in place.
    Input:
        trajectory:    MDTRAJ trajectory.
        fitframe:      Frame index of reference frame.
        reference:     Reference MDTRAJ trajectory. If none, the first trajectory is used as a reference.
        selectiontext: MDTRAJ selection string to select a set of atoms for RMSD fit.
    """
    if reference == None:
        reference = trajectory

    fit_atomndx = trajectory.topology.select(selectiontext)
    trajectory.superpose(trajectory, fitframe, atom_indices=fit_atomndx, parallel=parallel)


# ============================================================================ #


def _get_NH_bond_vectors(trajectory, Natomname='N', Hatomname='H', startresid=1):
    """
    Compute normalized NH bond vectors.
    Input:
        trajectory: MDTRAJ trajectory.
        Natomname:  Name of backbone N atoms.
        Hatomname:  Name of backbone H atoms.
        startresid: Index of first residue for bond vector computation (0 based).
    Return:
        NHvectors, NHbondlength, NHresids, NHresnames
    """

    # extract NH atom coordinates
    N_atomndx = trajectory.topology.select("name {} and not resname PRO and resid >= {}".format(Natomname, startresid))
    H_atomndx = trajectory.topology.select("name {} and not resname PRO and resid >= {}".format(Hatomname, startresid))
    N_trj     = trajectory.atom_slice(N_atomndx)
    H_trj     = trajectory.atom_slice(H_atomndx)

    # compute NH bond vectors
    NHvec = H_trj.xyz - N_trj.xyz

    # normalize NH bond vectors
    NHbondlength = (((NHvec**2).sum(2))**0.5)
    NHvec /= NHbondlength.reshape([NHvec.shape[0], NHvec.shape[1], 1])
    
    # get resids
    NHresids = []
    for res in list(N_trj.top.residues):
        NHresids.append([res.name, res.resSeq])
    NHresnames, NHresids = zip(*NHresids)

    return NHvec, NHbondlength, NHresids, NHresnames


# ============================================================================ #


def _compute_bond_vector_correlation_function(bondvec, verbose=True):
    """
    Compute 2nd order Legandre polynomial correlation functions of bond vectors.
    Input:
        nframes, nvectors, ndim  = bondvec.shape
    Return:
        correlation functions, standard deviations of correlation values, standard error of the mean
    """

    # compute NH bond vector correlation functions
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


def _check_corr_convergence(corr, diffThreshold=0.02, stdThreshold=0.02):
    """
    Check the convergence of the bond vector correlation functions
    Input:
        corr:          Correlation functions (shape = (nbonds, nframes)).
        diffThreshold: Maximum mean difference for convergence check.
        stdThreshold:  Maximum stdev difference for convergence check.
    Return:
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


def _n_exp(x, *args):
    """
    Evaluate the sum of n decaying exponentials depending on the length of *args.
    n coefficients and n exponential decay constants are expected.
    The coefficients are constraint to sum up to 1:
        y = sum_i(ai * exp(x/-ti)) + 1 - sum_i(ai)
    Example call:
        _n_exp(x, a1, a2, a3, t1, t2, t3)
    """
    coefficients = args[:len(args)/2]
    decayconst   = args[len(args)/2:]

    y = 0
    for i, c, t in zip(range(len(decayconst)), coefficients, decayconst):
        y += c * np.exp(x/(-1*t))
    y += 1 - sum(coefficients)

    return y


# ============================================================================ #


def _corr_fit(corr, dt=1, nexp=1, guess=None):
    """
    Fit a sum of exponentials to the correlation functions and report its parameters
    Input:
        corr:   Correlation functions of shape (nfunctions, timesteps)
        dt:     Delta t of subsequent correlation function values
        nexp:   Number of exponentials for fit
        guess:  Initial guess of coefficients and decay constants, if None, initialized to 1's.
    """
    xdata     = np.linspace(0, dt*corr.shape[1], corr.shape[1])
    coeffs    = np.zeros([corr.shape[0], nexp+1])
    coeffsstd = np.zeros([corr.shape[0], nexp+1])
    tconst    = np.zeros([corr.shape[0], nexp])
    tconststd = np.zeros([corr.shape[0], nexp])
    error     = np.zeros([corr.shape[0]])

    # Make np.exp() over and underflows throw an exception
    old_err_settings = np.seterr(over='ignore', under='ignore')

    for c in range(corr.shape[0]):
        print("\r", c, end="")
        sys.stdout.flush()

        # fit exponentials
        if guess == None:
            p0 = 2 * nexp * [1.0]
        else:
            p0 = guess
        successful     = False
        counter        = 1 
        nrandomguesses = 100
 
        while not successful:
            try: 
                popt, pcov = curve_fit(_n_exp, xdata, corr[c,:], p0=p0)
                if not float('NaN') in np.diag(pcov)**0.5:
                    successful = True
                else:
                    print("curve_fit failed with NaN error estimates")
            except:
                # change initial parameter guess
                if counter == 1 and c > 0:
                    # second try with optimized values from previous iteration
                    print("    trying previously optimized parameters as initial guess")
                    p0 = list(coeffs[c-1,:-1]) + list(tconst[c-1,:])
                elif counter == 2 and c > 0:
                    # third try with mean of optimized values from previous iterations
                    print("    trying mean of previously optimized parameters as initial guess")
                    p0 = list(coeffs[:c-1,:-1].mean(0)) + list(tconst[:c-1,:].mean(0)) 
                elif counter < 2 + nrandomguesses:
                    # continue trying random values
                    print("    trying random parameters as initial guess", counte)
                    p0 = np.random.rand(len(p0)) + 10**np.random.randint(0, 4, size=len(p0))
                else:
                    # failed
                    print("    failed to converge fitting")
                    popt = float('Nan') * np.ones([len(p0)         ])
                    pcov = float('Nan') * np.ones([len(p0), len(p0)])
                    successful = True

            counter += 1


        # extract fittet parameters and errors
        stdevs           = np.sqrt(np.diag(pcov))
        coeffs[c,:-1]    = popt[:len(popt)/2]
        coeffs[c, -1]    = 1 - sum(popt[:len(popt)/2])
        tconst[c,:  ]    = popt[len(popt)/2:]
        coeffsstd[c,:-1] = stdevs[:len(popt)/2]
        coeffsstd[c, -1] = (stdevs[:len(popt)/2]**2).sum()**0.5
        tconststd[c,:  ] = stdevs[len(popt)/2:]

        # sort by magnitude of decay constant
        srtIndx          = np.argsort(tconst[c,:])
        tconst   [c,:  ] = tconst   [c,srtIndx]
        coeffs   [c,:-1] = coeffs   [c,srtIndx]
        tconststd[c,:  ] = tconststd[c,srtIndx]
        coeffsstd[c,:-1] = coeffsstd[c,srtIndx]

        # compute error
        ncoeffs = coeffs.shape[1]-1
        prediction = _n_exp(xdata, *list(np.concatenate((coeffs[c,:-1].reshape(ncoeffs, 1), tconst[c,:].reshape(ncoeffs,1))).reshape(2*ncoeffs)))
        #error[c]   = abs(prediction - corr[c,:]).sum()
        error[c]   = ((prediction - corr[c,:])**2).mean()**0.5

    # restore error settings
    np.seterr(**old_err_settings)

    return coeffs, tconst, coeffsstd, tconststd, error


# ============================================================================ #


def running_std():
    pass


# ============================================================================ #


class SequencePlotter(object):

    def __init__(self, corr=None, coeffs=None, tconst=None, coeffsstd=None, tconststd=None, error=None, dt=1, indices=None, resnames=False, show=False, savedir=None):
        self.corr    = corr
        self.coeffs  = coeffs
        self.tconst  = tconst
        self.coeffsstd  = coeffsstd
        self.tconststd  = tconststd 
        self.error   = error
        self.dt      = dt
        self.resnames = resnames
        self.show    = show
        self.savedir = savedir
        self.c       = 0
        self.leftButtonPressed = False

        if type(indices) == type(None):
            self.indices = range(self.corr.shape[0])
        else:
            self.indices = indices

        self.lines   = []

        self.plot()


    def plot(self):
        self.c = 0

        # create figure and axis
        self.fig = plt.figure()
        self.axs = self.fig.add_subplot(111)  
        self.cid = self.fig.canvas.mpl_connect('button_press_event',   self.onclick)
        self.cid = self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.cid = self.fig.canvas.mpl_connect('motion_notify_event',  self.onmove)
        self.cid = self.fig.canvas.mpl_connect('scroll_event',         self.onscroll)

        self.redraw()
        plt.show()

    def redraw(self):
        index = self.indices[self.c]

        xdata  = np.linspace(0, self.dt*self.corr.shape[1], self.corr.shape[1])
        p      = list(self.coeffs[index,:-1]) + list(self.tconst[index,:])
        p_plus = list(self.coeffs[index,:-1] + self.coeffsstd[index,:-1]) + list(self.tconst[index,:] + self.tconststd[index,:])

        ydata = _n_exp(xdata, *p)
        ydata_plus = _n_exp(xdata, *p_plus)

        while len(self.lines) > 0:
            line = self.lines.pop()
            line.remove()
        self.lines += self.axs.plot(xdata, self.corr[index,:], 'b')
#        self.lines += self.axs.plot(xdata, self.corr[index,:], 'k')
        self.lines += self.axs.plot(xdata, ydata, 'r',
                label=r"$\tau$={:7.2f} ps, S2={:.2f}, error: {:.2e}".format(self.tconst[index,0], 1-self.coeffs[index,0], self.error[index]))
#        self.lines += self.axs.plot(xdata, ydata_plus, 'g')


        xmin, xmax = self.axs.get_xlim()
        self.lines += self.axs.plot([xmin, xmax], 2*[0.98], 'k', linewidth=2)
        self.lines += self.axs.plot(xmin + self.c*(xmax-xmin)/len(self.indices), 0.98, 'sk', markersize=15)

        resname = ""
        try:
            resname = self.resnames[index]
            self.axs.set_title(resname)
        except:
            self.axs.set_title("c={}, index={}".format(self.c, index))

        self.axs.set_xlabel("time [ps]")
        self.axs.set_ylabel("correlation")
        #self.axs.set_ylim(0, 1)
        self.axs.set_ylim(-0.5, 1)
        #self.axs.set_xscale('log')
        self.axs.legend(loc="lower left")
        self.fig.canvas.draw()

#        if savedir != None and os.path.isdir(savedir):
#            plt.savefig("{}/corr_{:03d}.png".format(savedir, c), dpi=300)
#
#        if show:
#            plt.show()
#        else:
#            plt.close()
 

    def onclick(self, event):
        if event.button == 1:
            self.leftButtonPressed = True
            try:
                xmin, xmax = self.axs.get_xlim()
                self.c = int(np.round((len(self.indices)-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.c < 0:
                    self.c = 0
                elif self.c >= len(self.indices):
                    self.c = len(self.indices) - 1
                self.redraw()
            except:
                pass

    def onrelease(self, event):
        if event.button == 1:
            self.leftButtonPressed = False

    def onmove(self, event):
        if self.leftButtonPressed:
            try:
                xmin, xmax = self.axs.get_xlim()
                self.c = int(np.round((len(self.indices)-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.c < 0:
                    self.c = 0
                elif self.c >= len(self.indices):
                    self.c = len(self.indices) - 1 
                self.redraw()
            except:
                pass
 


    def onscroll(self, event):
        self.c -= int(event.step)
        if self.c < 0:
            self.c = len(self.indices) - 1
        elif self.c >= len(self.indices):
            self.c = 0

        self.redraw()






# ============================================================================ #


def bondvec_corr(trj, bondvec=None, fitgroup=None, parallel=True, verbose=True):
    """
    Compute bond vector correlation functions from mdtraj trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to compute the order parameters of.

    bondvec : 2 element list/tuple, optional
        Bothe elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.

    parallel : boolean, optional
        Do the fitting in parallel on as many processors as available on the machine.

    verbose : boolean, optional
        Toggle verbosity level.


    Returns
    -------
    corr : (nvec, t) array
        Bond vector correlation functions.
    info :
        Dictionary with information on the bond vectors.
        Keys are: ['bondlength', 'resname', 'resids']
    """

# ============================================================================ #


def compute_S2(trajectory=None, filenames=None, verbose=True, parallel=True, fit=True):
    """
    Compute S2 order parameters from an MDTRAJ trajectory
    Input:
        Either an MDTRAJ trajectory or a dictionary: {"top": topfilename, "trj": trjfilename}
    Return:
        S2 values
        boolean array indicating convergence
        NHinfo dictionary with information on bondlength resids and resnames
        NHvec correlation functions
    """

    if trajectory == None:
        trajectory = mdtraj.load(filenames["trj"], top=filenames["top"])

    # superpose trajectory on first frame based on C alhpa atoms
    if fit:
        _fit_trajectory(trajectory, fitframe=0)

    # compute NH bond vectors
    NHvec, NHbondlength, NHresids, NHresnames = _get_NH_bond_vectors(trajectory)
    NHinfo = {"bondlength": NHbondlength, "resids": NHresids, "resnames": NHresnames}

    # compute correlation functions of NH bond vectors
    corr, corr_std, corr_stdmean = _compute_bond_vector_correlation_function(NHvec, verbose=verbose)

    # check convergence of correlation functions and compute S2 values
    S2, convergence = _check_corr_convergence(corr)

    return S2, convergence, NHinfo, corr


# ============================================================================ #


def compute_S2_batch(trajectorylist=None, filenamelist=None, only_converged=True, verbose=True):
    """
    Batch version of compute_S2().
    Processes a list of trajectories or a list of {top, trj} filename dictionaries
    and computes the average S2 and standard deviation.
    Input:
        Either a list of MDTRAJ trajectories or a list of {top, trj} filename dictionaries.
        only_converged: Use only converged bond vector correlation functions for averaging.
    Return:
        Mean S2 values
        standard deviations
        number of converged correlation funtions
        NHinfo dictionary with information on bondlength resids and resnames
        NHcorr correlation functions of shape (nsims, nvec, nframes/2)
    """
    
    S2list        = []
    convergedlist = []
    corrlist      = []

    # compute S2 for each individual MD Dataset
    if trajectorylist == None:
        for trj_ndx, filenames in enumerate(filenamelist):
            if verbose:
                print("Dataset {} of {} ({})".format(trj_ndx+1, len(filenamelist), filenames["trj"]))
            nextS2, nextConverged, NHinfo, corr = compute_S2(filenames=filenames, verbose=verbose)
            S2list.append(nextS2)
            convergedlist.append(nextConverged)
            corrlist.append(corr)
    else:
        for trj_ndx, trajectory in enumerate(trajectorylist):
            if verbose:
                print("Dataset {} of {}".format(trj_ndx+1, len(trajectorylist)))
            nextS2, nextConverged, NHinfo, corr = compute_S2(trajectory=trajectory, verbose=verbose)
            S2list.append(nextS2)
            convergedlist.append(nextConverged) 
            corrlist.append(corr)
        

    # Compute mean, std, etc.
    S2mean     = np.zeros_like(S2list[0])
    S2array    = np.zeros([S2mean.shape[0], len(S2list)], dtype=np.float)
    nconverged = np.zeros_like(convergedlist[0], dtype=np.int)
    NHcorr     = np.zeros([len(S2list), corrlist[0].shape[0], corrlist[0].shape[1]], dtype=np.float)

    for S2_idx, nextS2 in enumerate(S2list):
        nextConverged = convergedlist[S2_idx]
        if only_converged:
            S2mean[nextConverged] += nextS2[nextConverged]
        else:
            S2mean += nextS2
        nconverged += nextConverged
        S2array[:,S2_idx]  = nextS2
        NHcorr[S2_idx,:,:] = corrlist[S2_idx]

    largerZero = nconverged > 0
    S2mean[largerZero] /= nconverged[largerZero] 
    S2std = S2array.std(1)

    return S2mean, S2std, nconverged, NHinfo, NHcorr


# ============================================================================ #


def compute_S2_batch_mpi(trajectorylist=None, filenamelist=None, only_converged=True):
    """
    Batch version of compute_S2() with mpi4py parallelization.
    Processes a list of trajectories or a list of {top, trj} filename dictionaries
    and computes the average S2 and standard deviation.
    Input:
        Either a list of MDTRAJ trajectories or a list of {top, trj} filename dictionaries.
        only_converged: Use only converged bond vector correlation functions for averaging.
    Return:
        Mean S2 values
        standard deviations
        number of converged correlation funtions
        NHinfo dictionary with information on bondlength resids and resnames
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    S2list        = []
    convergedlist = []
    corrlist      = []

    # compute S2 for each individual MD Dataset
    if trajectorylist == None:
        for trj_ndx, filenames in enumerate(filenamelist):
            if trj_ndx % size == rank:
                print(rank, trj_ndx); sys.stdout.flush()
                nextS2, nextConverged, NHinfo, corr = compute_S2(filenames=filenames, verbose=False, parallel=False)
                S2list.append(nextS2)
                convergedlist.append(nextConverged)
                corrlist.append(corr)
            elif rank == 0:
                S2list.append(trj_ndx % size) # rank of process that rank 0 should receive this part of the data from
    else:
        for trj_ndx, trajectory in enumerate(trajectorylist):
            if trj_ndx % size == rank:
                nextS2, nextConverged, NHinfo, corr = compute_S2(trajectory=trajectory, verbose=False, parallel=False)
                S2list.append(nextS2)
                convergedlist.append(nextConverged) 
                corrlist.append(corr)
            elif rank == 0:
                S2list.append(trj_ndx % size) # rank of process that rank 0 should receive this part of the data from
        

    # Compute mean, std, etc.
    if rank == 0:
        S2mean     = np.zeros_like(S2list[0])
        S2array    = np.zeros([S2mean.shape[0], len(S2list)], dtype=np.float)
        nconverged = np.zeros_like(convergedlist[0], dtype=np.int)
        NHcorr     = np.zeros([len(S2list), corrlist[0].shape[0], corrlist[0].shape[1]], dtype=np.float)

        for S2_idx, nextS2 in enumerate(S2list):

            if type(nextS2) == type(int()):
                sendingRank   = nextS2
                nextS2        = np.zeros_like(S2list[0])
                nextConverged = np.zeros_like(convergedlist[0])
                nextCorr      = np.zeros_like(corrlist[0])
                comm.Recv(nextS2,        source=sendingRank, tag=0)
                comm.Recv(nextConverged, source=sendingRank, tag=1)
                comm.Recv(nextCorr,      source=sendingRank, tag=2)

            else:
                nextConverged = convergedlist[S2_idx]
                nextCorr      = corrlist[S2_idx]

            if only_converged:
                S2mean[nextConverged] += nextS2[nextConverged]
            else:
                S2mean += nextS2
            nconverged += nextConverged
            S2array[:,S2_idx] = nextS2
            NHcorr[S2_idx,:,:] = nextCorr

        largerZero = nconverged > 0
        S2mean[largerZero] /= nconverged[largerZero]
        S2std = S2array.std(1)

        return S2mean, S2std, nconverged, NHinfo, NHcorr

    else:
        for S2_idx, nextS2 in enumerate(S2list):
            comm.Send(nextS2,                dest=0, tag=0)
            comm.Send(convergedlist[S2_idx], dest=0, tag=1)
            comm.Send(corrlist[S2_idx],      dest=0, tag=2)
            

        return None, None, None, None, None
 



