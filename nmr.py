# Compute NMR observables from MD trajectories

from __future__ import print_function

import sys
import numpy as np

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

 
 


def S2(trajectory, verbose=True):
    """Compute S2 order parameters from an MDTRAJ trajectory"""

    # superpose trajectory on first frame based on C alhpa atoms
    _fit_trajectory(trajectory, fitframe=0)

    # compute NH bond vectors
    NHvec, NHbondlength, NHresids, NHresnames = _get_NH_bond_vectors(trajectory)

    # compute correlation functions of NH bond vectors
    corr, corr_std, corr_stdmean = _compute_bond_vector_correlation_function(NHvec, verbose=verbose)

    # check convergence of correlation functions and compute S2 values
    S2, convergence = _check_corr_convergence(corr)

    return S2, convergence





