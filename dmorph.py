#!/usr/bin/env python

from __future__ import print_function

import sys, time
import numpy as np
import mdtraj as md
from scipy.constants import pi

# ============================================================================ #

def angular_diff(angle1, angle2):
    """
    Compute the shortest angle to rotate from angle1 to angle2 in radian.
    """
    diff = angle2 - angle1
    try:
        diff[diff > 180]  = diff[diff > 180] - 360
        diff[diff < -180] = diff[diff < -180] + 360 
        return diff
    except TypeError:
        diff = np.array([diff])
        diff[diff > 180]  = diff[diff > 180] - 360
        diff[diff < -180] = diff[diff < -180] + 360  
        return diff[0]
 
# ============================================================================ #

class Dihedral(object):
    """
    Contains all information on a dihedral necessary for a rotation.
    """
    def __init__(self, trj, indices, rotIndices, initAngle, finalAngle):
        """
        Contains all information on a dihedral necessary for a rotation.

        Parameter
        ---------
        trj: mdtraj trajectory

        indices: list
            List of atom indices in order of dihedral: 1-2-3-4
        initAngle: float
            initial rotation angle
        finalAngle: float
            final rotation angle 
        rotIndices: array (n,)
            Indices of all atoms that should be rotated by this dihedral
        """
        self.indices     = indices
        self.atomnames   = []
        self.rotIndices  = rotIndices
        self.initAngle   = initAngle
        self.finalAngle  = finalAngle
        self.diff        = angular_diff(initAngle, finalAngle)

        for idx in self.indices:
            self.atomnames.append(str(trj.top.atom(idx)))

    # ======================================================================== #

    def rot_angle(self, f, mode="linear"):
        """
        Compute rotation angle from fraction f of total rotation angle
        f = t/N

        Possible modes:
            {"linear", "cycle", "sin2"}
        """
        if mode == "linear":
            return f * self.diff
        elif mode == "cycle":
            return (1 - np.cos(f*pi)) / 2
        elif mode == "sin2":
            return np.sin(f*pi/2)**0.5
        else:
            print("\r", 100*" ", "\rinterpolation mode unknown:", mode)
            sys.exit(1)
 
    # ======================================================================== #

    def __str__(self):
        return str(self.atomnames)

# ============================================================================ #

def binSum(n):
    """
    Decompose integer into powers of 2.
    e.g.: 43 = 32 + 8 + 2 + 1

    Parameters
    ----------
    n: int
        number

    Returns
    -------
    a: array (m, )
        Array containing the components of n
    """

    a = []
    number   = n
    exponent = 0
    while sum(a) < n:
        if (number % 2**(exponent+1)) > 0:
            a.append(2**exponent)
            number -= 2**exponent
        exponent += 1

    return np.array(a)

# ============================================================================ #

def allocate_trj(first, nframes):
    """
    Allocate space for the morph trajectory and fill all frames with coordinates
    of first frame.
    Powers of 2 allocate particularly fast.
    """

    trj = None
    components = binSum(nframes)
    for c in components:
        nexttrj = md.load(first)
        while nexttrj.n_frames < c:
            nexttrj += nexttrj
        if trj is None:
            trj = nexttrj
        else:
            trj += nexttrj

    return trj

# ============================================================================ #

def ramachandran_angles(first, last):
    """
    Collect information on ramachandran dihedrals

    Parameters
    ----------
    first: string
        Filename of initial structure
    last: string
        Filename of final structure
    """
    dihedrals = []

    firsttrj = md.load(first)
    lasttrj  = md.load(last)

    phi_ndx, firstPhi = md.compute_phi(firsttrj)
    psi_ndx, firstPsi = md.compute_psi(firsttrj)
    phi_ndx,  lastPhi  = md.compute_phi(lasttrj)
    psi_ndx,  lastPsi  = md.compute_psi(lasttrj)

    for angle_idx in range(phi_ndx.shape[0]):
        resid       = firsttrj.top.atom(phi_ndx[angle_idx,1]).residue.index
        rotIndices  = list(firsttrj.top.select("resid {} and (sidechain or name C or name O) and not name H".format(resid)))
        rotIndices += list(firsttrj.top.select("resid > {}".format(resid)))
        dihedrals.append(Dihedral(firsttrj, phi_ndx[angle_idx,:], rotIndices, firstPhi[0,angle_idx], lastPhi[0,angle_idx]))

    for angle_idx in range(psi_ndx.shape[0]):
        resid       = firsttrj.top.atom(psi_ndx[angle_idx,1]).residue.index
        rotIndices  = list(firsttrj.top.select("resid {} and name O".format(resid)))
        rotIndices += list(firsttrj.top.select("resid > {}".format(resid)))
        dihedrals.append(Dihedral(firsttrj, psi_ndx[angle_idx,:], rotIndices, firstPsi[0,angle_idx], lastPsi[0,angle_idx]))
 
    return dihedrals

# ============================================================================ #

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
        rotation angle

    Returns
    -------
    M : array(4,4)
        4 dimensionl rotation and translation matrix
    """

    M = np.zeros([4,4], dtype=np.float64)
    a, b, c = r
    u, v, w = v
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
 
# ============================================================================ #

def rotate_dihedrals(xyz4, f, dihedrals, mode="linear"):
    """
    Rotate all dihedral angles of a frame.

    Parameter
    ---------
    xyz4: float array (natoms, 4)
        4 dimensional position vectors
    f: float
        Fraction between 0 (start) and 1 (end) of rotation
    dihedrals: list
        List of dihedrals angles to be rotated
    mode: string, optional
        Rotation mode
    """

    rottime = 0.0

    for dihedral in dihedrals:
        rotAngle   = dihedral.rot_angle(f, mode=mode)
        rotPoint   = xyz4[dihedral.indices[1],:3]
        rotVector  = xyz4[dihedral.indices[2],:3] - rotPoint
        rotVector /= np.linalg.norm(rotVector)

        M = rotation_matrix_arbitrary_axis(rotPoint, rotVector, rotAngle)

        start_t = time.time()
        for rot_idx in dihedral.rotIndices:
            xyz4[rot_idx,:] = np.dot(M, xyz4[rot_idx,:])
        rottime += time.time() - start_t


    return rottime

# ============================================================================ #

def dmorph(first, last, nframes, outfile, mode="linear"):
    """
    Linearly interpolate the dihedral angels from firststructure to laststructure.

    Parameters
    ----------
    first: string
        Starting structure for morph. PDB filename
    last: string
        Last structure for morph. PDB filename. 
    nframes: int
        Number of frames for morph.
    outfile: string
        Path and filename for the output trajectory
    mode: string
        Sets the interpolation mode between first and last.
        Mode is one of:
            {"linear", "cycle", "sin2"}
    """

    trj  = allocate_trj(first, nframes)
    xyz4 = np.ones([trj.n_frames, trj.n_atoms, 4], dtype=np.float64)
    xyz4[:,:,:3] = trj.xyz

    dihedrals = ramachandran_angles(first, last)

    phi_ndx, targetPhi = md.compute_phi(md.load(last))
    psi_ndx, targetPsi = md.compute_psi(md.load(last))

    rottime = 0.0
    start_t = time.time()

    error = np.zeros([nframes])
    for nf in range(1, nframes):
        print("\r", 100*" ", "\r", end="")
        print("frame {:5d} / {:5d}".format(nf+1, nframes), end="")
        rottime += rotate_dihedrals(xyz4[nf,:,:], 1.0*nf/nframes, dihedrals, mode=mode)
        sys.stdout.flush()

        trj.xyz = xyz4[:,:,:3]

        phi_ndx, phi = md.compute_phi(trj)
        psi_ndx, psi = md.compute_psi(trj)

        e = (((psi[nf,:] - targetPsi[0,:])**2).sum()/psi.shape[1])**0.5
        error[nf] = e
        print(" ", e)

    trj.superpose(trj)

    tottime = time.time() - start_t

    print()
    print("Runtime: {:6.2f} sec.".format(tottime))
#    print("rottime: {:6.2f}%".format(100*rottime/tottime))

#    lasttrj = md.load(last)
#    phi_ndx, targetPhi = md.compute_phi(lasttrj)
#    phi_ndx, targetPsi = md.compute_psi(lasttrj)
#    phi_ndx, phi = md.compute_phi(trj)
#    psi_ndx, psi = md.compute_phi(trj)
#
#    for nf in range(phi.shape[0]):
#        error[nf]  = 0.0
#        error[nf] += ((phi[nf,:] - targetPhi[0,:])**2).sum()
#        error[nf] += ((psi[nf,:] - targetPsi[0,:])**2).sum()
#        error[nf] /= 2*phi.shape[1]
#        error[nf]  = error[nf]**0.5
#    error = error

    trj.save(outfile)

    return error

# ============================================================================ #

if __name__ == "__main__":

    try:
        first   = sys.argv[1]
        last    = sys.argv[2]
        outfile = sys.argv[3]
    except IndexError:
        print("Should be called as: dmorph.py first.pdb last.pdb out.xtc [nframes] [mode]")
        sys.exit(1)
    
    nframes = 128
    mode    = "linear"

    if len(sys.argv) > 4:
        nframes = int(sys.argv[4])
    if len(sys.argv) > 5:
        mode = sys.argv[5]

    dmorph(first, last, nframes, outfile, mode=mode)







