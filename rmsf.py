# Compute rmsf of protein structures
from __future__ import print_function

import sys
from time import time as t

def _atom_indices_from_selection(topology, selection):
    """
    Compute atom Indices from selection text
    Input:
        topology:    A mdtraj.topology object.
        selection:   Either a selection string from mdtraj.topology.select() or atom indices.
                     If indices are given, these are simply returned unchanged.
    """
    if type(selection) in [type(str()), type(unicode())]:
        atomIndices = topology.select(selection)
    else:
        atomIndices = selection
    return atomIndices
 

def _fit_trajectory(trajectory, reference, selection, openMP=True):
    """
    Fit a trajecory to a reference frame
    Input:
        trajectory:   mdtraj trajectory instance
        reference:    mdtraj trajectory of which frame 0 is taken as reference
        selection:    Either a selection text for fitting or atom indices
        openMP:       Use openMP paralellization or not
    """
    atomIndices = _atom_indices_from_selection(trajectory.topology, selection)
    trajectory.superpose(reference, frame=0, atom_indices=atomIndices, parallel=openMP)



def _sumsquares(trajectory, reference, selection, refMean=False):
    """
    Compute sum of squares of atomic position difference around reference.
    Input:
        trajectory:   mdtraj trajectory instance
        reference:    mdtraj trajectory of which frame 0 is taken as reference
        selection:    Either a selection text for or atom indices 
        refMean:      If true, take the mean atom positions from trajectory as reference (reference can be set to None)
    """

    atomIndices = _atom_indices_from_selection(trajectory.topology, selection)

    if refMean:
        refpos = trajectory.xyz[:,atomIndices,:].mean(0)
    else:
        refpos = reference.xyz[0,atomIndices,:]

    diff = trajectory.xyz[:,atomIndices,:] - refpos
    sumsquares = ((diff**2).sum(2)).sum(0)

    return sumsquares
    



def rmsf(trajectories, reference, fitselection, rmsfselection, maxChunk=None, refMean=False, verbose=True):
    """
    Compute the RMSF of a number of trajectories to a reference frame
    Input:
        trajectories:    mdtraj.iterload() instance
        reference:       trajectory of which frame 0 is taken as a reference frame
        fitselection:    Selection text for rmsd fit or atom indices
        rmsdselection:   Selection text for rmsf computation or atom indices
        maxChunk:        Maximum number of trajectory chunks to process
        refMean:         If true, take the mean atom positions from trajectory as reference (reference can be set to None)
    Return:
        RMSF:            Root mean square fluctuations
        atomNdx:         Atom indices of selcted atoms
        resID:           Resiude IDs of selected atoms
    """

    # timers
    fitt    = 0.0
    sst     = 0.0
    loadt   = 0.0

    nframes = 0     # total number of frames

    loadst = t()

    for ntrj, trajectory in enumerate(trajectories):
        loadt += t() - loadst

        if verbose:
            print("\rProcessing trajectory chunk {} of maximally {}.".format(ntrj, maxChunk), end="")
            sys.stdout.flush()

        if type(maxChunk) != type(None) and ntrj >= maxChunk:
            break

        st = t()
        _fit_trajectory(trajectory, reference, fitselection, openMP=False)
        fitt += t() - st

        st = t()
        try:
            sumsquares += _sumsquares(trajectory, reference, rmsfselection, refMean=refMean)
        except:
            sumsquares  = _sumsquares(trajectory, reference, rmsfselection, refMean=refMean)
        sst += t() - st

        nframes += trajectory.n_frames
        loadst = t()

    if verbose:
        print("\rLoading trajectories with a total of {} frames took {:.2f} sec.".format(nframes, loadt))
        print("RMSD fit took {:.2f} sec.".format(fitt))
        print("RMSF computation took {:.2f} sec.".format(sst ))

    RMSF = (sumsquares / (3*nframes))**0.5

    atomIndices    = _atom_indices_from_selection(trajectory.topology, rmsfselection)
    atomIndicesTrj = trajectory.atom_slice(atomIndices)
    resIDs         = [atom.residue.resSeq for atom in atomIndicesTrj.topology.atoms]

    return RMSF, atomIndices, resIDs
 
