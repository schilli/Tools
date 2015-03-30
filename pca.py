import pca
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt


def linear_testdata(n=100):
    x = np.linspace(0,1,n)
    y = x + 0.01 * np.random.randn(n)
    data = np.zeros([2, len(x)])
    data[0,:] = x
    data[1,:] = y
    return data


def test(n=100):

    data = linear_testdata(n)

    w, v = pca(data)

    print v

    plt.plot(data[0,:], data[1,:], '.')
    plt.plot([0, v[0,0]], [0, v[1,0]], 'r')
    plt.plot([0, v[0,1]], [0, v[1,1]], 'g')
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 1.5)
    plt.show()




def pca(data):
    """Compute the principal components of the given data based on eigenvalue decomposition
    Input:
        data: First index (rows) are variables,
              second index (columns) are observations
    Returns:
        (eigvals, eigvecs)
        Eigenvalues are sorted in decreasing order,
        The corresponding eigenvectors are the columns of the second returned array.
    """

    # diagonalize covariance matrix
    covariance = np.cov(data)
    eigvals, eigvecs = np.linalg.eig(covariance)
    sortedeigvecs = np.copy(eigvecs)

    # sort eigenvectors by eigenvalue in decreasing order
    tmp = [[eigvals[i], eigvecs[:,i]] for i in range(len(eigvals))]
    tmp.sort()
    i = len(eigvals)
    for eigval, eigvec in tmp:
        i -= 1
        eigvals      [  i] = eigval
        sortedeigvecs[:,i] = eigvec

    return eigvals, sortedeigvecs



def dpca(dihedrals):
    """Perform a dihedral pca
    Input:
        Dihedral angle data: rows (first index) are the angles and
                             columns (second index) observations
    Returns:
        Array of coordinates in the space spanned by the dihedral principal components
    """

    # Create cartesian coordinate space of x = cos(phi), y = sin(phi)
    cartcoords = np.zeros([2 * dihedrals.shape[0], dihedrals.shape[1]], dtype=np.float)
    cosines = np.cos(const.pi / 180.0 * dihedrals)
    sines   = np.sin(const.pi / 180.0 * dihedrals)
    cos_idx = np.arange(0,2*dihedrals.shape[0],2)
    sin_idx = np.arange(1,2*dihedrals.shape[0],2)
    cartcoords[cos_idx,:] = cosines
    cartcoords[sin_idx,:] = sines

    # Compute pca
    eigvals, eigvecs = pca(cartcoords)

    # Project data on principal components
    projectedcoords = np.zeros_like(cartcoords)
    for point_idx in range(cartcoords.shape[1]):
        for eig_idx in range(eigvecs.shape[0]):
            projectedcoords[eig_idx, point_idx] = np.dot(eigvecs[:,eig_idx], cartcoords[:,point_idx])

    return eigvals, projectedcoords
    



def fes(projectedcoords, dim=2, T=300, nbins=5):
    """Estimate the free energy landscape for the given number of dimensions
    Ei = k*T * [ln(Nmax) - ln(Ni)]
    """
    if dim > 3 or dim > projectedcoords.shape[0]:
        print "dim must not be larger than 3 or the dimension of input coordinates"""
        sys.exit(1)

    coords = projectedcoords[:dim,:]
    maxima = coords.max(1)
    minima = coords.min(1)
    deltas = np.zeros_like(maxima)

    # construct coordinates
    coordRange = []
    for dim_idx in range(dim):
        coordRange.append(np.linspace(minima[dim_idx], maxima[dim_idx], nbins))
        deltas[dim_idx] = abs(coordRange[-1][0] - coordRange[-1][1])
    coordSpace = np.meshgrid(*coordRange)

    # count bins
    counts = np.zeros_like(coordSpace[0], dtype=np.int)
    binIndex = dim * [0]
    for point_idx in range(coords.shape[1]):
        for dim_idx in range(dim):
            binIndex[dim_idx] = int((coords[dim_idx, point_idx] - (minima[dim_idx] - deltas[dim_idx]/2.0)) / deltas[dim_idx])
            #print point_idx, dim_idx, binIndex[dim_idx]
        counts[tuple(binIndex)] += 1

    # estimate FES
    Nmax     = counts.max()
    lnCounts = np.log(counts+1)
    FES = const.k * T * (np.log(Nmax) - lnCounts)

    return coordSpace, FES


    



    


        
    
    











