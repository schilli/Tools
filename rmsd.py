import numpy as np
from scipy.linalg import block_diag

def make_Nx3(P):
    """reshape dataset of points to be of shape Nx3"""
    if len(P.shape) == 1 or P.shape[1] != 3:
        N = P.shape[0]/3
        P = P.reshape(N, 3) 
    return P


def center_points(P):
    """shift the geometric center of the set of points P to the origin
    returns (centeredPoints, center)"""
    P        = make_Nx3(P)
    center   = P.mean(0)
    centered = P - np.meshgrid(center, P.shape[0])[0]
    return centered, center


def rmsd_matrix(P, Q, centered=False):
    """Compute the rotation matrix that minimizes the
    rmsd between points P and Q with Kabsch algorithm.
    Before applying the rotation matrix, structures need to be
    centered at the origin.
    P and Q are either of shape Nx3 or 3Nx1.
    Refer to http://en.wikipedia.org/wiki/Kabsch_algorithm"""

    # ensure proper shape
    P = make_Nx3(P)
    Q = make_Nx3(Q)

    # center if necessary
    if not centered:
        P, Pcenter = center_points(P)
        Q, Qcenter = center_points(Q)
    else:
        Pcenter = np.zeros(3)
        Qcenter = np.zeros(3)

    # compute covariance matrix
    A = np.dot(Q.T, P)

    # singular value decomposition
    V,S,W = np.linalg.svd(A)

    # compute optimal rotation matrix
    d = np.linalg.det(np.dot(W.T, V.T))
    if d < 0:
        d = -1
    else:
        d = +1
    D = np.diag([1, 1, d])
    T = np.dot(D, V.T)
    U = np.dot(W.T, T)

    # done
    return U


def superimpose(P, Q):
    """superimpose the set of atoms Q on P"""
    # ensure proper shape
    P = make_Nx3(P)
    Q = make_Nx3(Q)
    N = P.shape[0]

    # center if necessary
    P, Pcenter = center_points(P)
    Q, Qcenter = center_points(Q) 

    U = rmsd_matrix(P, Q, centered=True)

    #R = block_diag(*(N*[U]))
    #Qsuper = np.dot(R, Q.reshape(3*N,1))

    Qsuper = Q.reshape(3*N, 1)
    for n in range(N):
        Qsuper[n*3:n*3+3] = np.dot(U, Qsuper[n*3:n*3+3])

    Qsuper = Qsuper.reshape(N,3)
    Qsuper += np.meshgrid(Pcenter, Q.shape[0])[0]
    Psuper  = P + np.meshgrid(Pcenter, Q.shape[0])[0]

    return Psuper, Qsuper



def superimpose_ref(P, Q, Pref, Qref):
    """Superimpose structure Q on P based on the coordinates in Pref and Qref"""

    # ensure proper shape
    P    = make_Nx3(P)
    Q    = make_Nx3(Q)
    Pref = make_Nx3(Pref)
    Qref = make_Nx3(Qref) 
    N    = P   .shape[0]
    Nref = Pref.shape[0]

    # center if necessary
    P,    Pcenter = center_points(P)
    Q,    Qcenter = center_points(Q) 
    Pref, Pcenter = center_points(Pref)
    Qref, Qcenter = center_points(Qref)  

    U = rmsd_matrix(Pref, Qref, centered=True)

    Qsuper = Q.reshape(3*N, 1)
    for n in range(N):
        Qsuper[n*3:n*3+3] = np.dot(U, Qsuper[n*3:n*3+3]) 

    Qsuper = Qsuper.reshape(N,3)
    Qsuper += np.meshgrid(Pcenter, Q.shape[0])[0]
    Psuper  = P + np.meshgrid(Pcenter, Q.shape[0])[0]

    return Psuper, Qsuper 



def rmsd(P, Q, lsfit=True):
    """Compute the RMSD of two sets of points
    lsfit determines if a least squares fitting is to be done"""

    N = P.size / 3 # number of particlesA

    if lsfit:
        Psuper, Qsuper = superimpose(P, Q)
    else:
        Psuper, Qsuper = (P, Q)

#    diff = P.reshape(3*N) - Qsuper.reshape(3*N)
#    RMSD = np.sqrt((diff**2).mean())
    diff = P.reshape(N,3) - Qsuper.reshape(N,3)
    RMSD = 0
    for i in range(N):
        RMSD += np.linalg.norm(diff[i,:])**2
    RMSD = np.sqrt(RMSD/N) 
    return RMSD
