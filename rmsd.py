import numpy as np
from scipy.linalg import block_diag

def make_Nx3(P):
    """reshape dataset of points to be of shape Nx3"""
    if len(P.shape) == 1 or P.shape[1] != 3:
        N = P.shape[0]/3
        P = P.reshape(N, 3) 
    return P


def center(P):
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
        P, Pcenter = center(P)
        Q, Qcenter = center(Q)
    else:
        Pcenter = np.zeros(3)
        Qcenter = np.zeros(3)

    # compute covariance matrix
    A = np.dot(P.T, Q)

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


def superimpose(P, Q, centered=False):
    """superimpose the set of atoms Q on P"""
    # ensure proper shape
    P = make_Nx3(P)
    Q = make_Nx3(Q)
    N = P.shape[0]

    # center if necessary
    if not centered: 
        P, Pcenter = center(P)
        Q, Qcenter = center(Q) 

    U = rmsd_matrix(P, Q, centered=True)

    #R = block_diag(*(N*[U]))
    #Qsuper = np.dot(R, Q.reshape(3*N,1))

    Qsuper = Q.reshape(3*N, 1)
    for n in range(N):
        Qsuper[n*3:n*3+3] = np.dot(U, Qsuper[n*3:n*3+3])

    Qsuper = Qsuper.reshape(N,3)
    Qsuper += np.meshgrid(Qcenter, Q.shape[0])[0]

    return Qsuper



def rmsd(P, Q, superimposed=False, centered=False):
    """Compute the RMSD of two sets of points"""

    N = P.size / 3 # number of particlesA

#    if not superimposed:
#        Qsuper = superimpose(P, Q, centered)
#    else:
#        Qsuper = Q
    Qsuper = Q

    Qsuper = Qsuper.reshape(3*N, 1)
#    print N
#    print P.shape
#    print Q.shape
#    print Qsuper.shape
    diff = P - Qsuper
#    diff = P - Qsuper.reshape(3*N, 1)
#    RMSD = np.linalg.norm(diff.reshape(3*N, 1))
#    return RMSD
