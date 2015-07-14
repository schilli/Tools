# Compute NMR observables from MD trajectories

from __future__ import print_function

import sys
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


def _single_exp(x, a, tau):
    """
    Single exponential: y = (1-a) * exp(x/(-1*tau)) + a
    """
    return (1-a) * np.exp(x/(-1*tau)) + a


# ============================================================================ #


def _double_exp(x, a, b, tau1, tau2):
    """
    Double exponential: y = (1-a-b) * exp(x/(-1*tau1)) + a * exp(x/(-1*tau2)) + b
    """
    return (1-a-b) * np.exp(x/(-1*tau1)) + a * np.exp(x/(-1*tau2)) + b


# ============================================================================ #


def _double_exp2(x, a, tau1, tau2):
    return (0.75-a) * np.exp(x/(-1*tau1)) + a * np.exp(x/(-1*tau2)) + 0.25


# ============================================================================ #


def _n_exp(x, *args):
    """
    n exponentials depending on the length of *args.
    n coefficients and n exponential decay constants are expected.
    Example:
    _n_exp(x, a, b, c, t1, t2, t3)
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

    # Make np.exp() over and underflows throw an exception
    old_err_settings = np.seterr(over='ignore', under='ignore')

    for c in range(corr.shape[0]):
        print("\r", c, end="")

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

    # restore error settings
    np.seterr(**old_err_settings)

    return coeffs, tconst, coeffsstd, tconststd



#        if fit == 'single':
#            popt, pcov = curve_fit(_single_exp, xdata, corr[c,:], p0=[0.5, 100])
#            S2[c]      = popt[0]
#            tau1[c]    = popt[1]
#            S2std[c], tau1std[c] = np.sqrt(np.diag(pcov))
#        elif fit == 'double':
#            popt, pcov = curve_fit(_double_exp, xdata, corr[c,:], p0=[0.1, 0.1, 100, 5])
#            a[c]       = popt[0]
#            S2[c]      = popt[1]
#            tau1[c]    = popt[2]
#            tau2[c]    = popt[3]
#            astd[c], S2std[c], tau1std[c], tau2std[c] = np.sqrt(np.diag(pcov)) 
#        elif fit == 'double2':
#            popt, pcov = curve_fit(_double_exp2, xdata, corr[c,:], p0=[0.75/2, 100, 5])
#            a[c]       = popt[0]
#            tau1[c]    = popt[1]
#            tau2[c]    = popt[2]
#            astd[c], tau1std[c], tau2std[c] = np.sqrt(np.diag(pcov))  
#        else:
#            print("fit must be single or double")
#            sys.exit(1)
#
#    return a, S2, tau1, tau2, astd, S2std, tau1std, tau2std
 


# ============================================================================ #


def _plot_fit(corr, coeffs, tconst, dt=1, c=0):

    fig = plt.figure()
    axs = fig.add_subplot(111) 

    if c == None:
        xdata = np.linspace(0, dt*corr.shape[0], corr.shape[0])
        axs.plot(xdata, corr)
        p = list(coeffs[:-1]) + list(tconst)
        axs.plot(xdata, _n_exp(xdata, *p), label=r"$\tau$={:.2e} ps, S2={:.2f}".format(tconst[0], 1-coeffs[0]))
    else:
        xdata = np.linspace(0, dt*corr.shape[1], corr.shape[1])
        axs.plot(xdata, corr[c,:])
        p = list(coeffs[c,:-1]) + list(tconst[c,:])
        axs.plot(xdata, _n_exp(xdata, *p), label=r"c={}, $\tau$={:.2e} ps, S2={:.2f}".format(c, tconst[c,0], 1-coeffs[c,0]))

    axs.legend()
    #plt.show()
    plt.savefig("plots/corr_{:03d}.png".format(c), dpi=300)
    #plt.savefig("plots/corr.png", dpi=300)
    plt.close()


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
 



