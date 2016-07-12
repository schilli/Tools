from __future__ import print_function

import pca
import sys, time
import numpy as np
import sklearn.cluster as skc
import scipy.constants as const
import matplotlib.pyplot as plt
from collections import Counter


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

    print(v)

    plt.plot(data[0,:], data[1,:], '.')
    plt.plot([0, v[0,0]], [0, v[1,0]], 'r')
    plt.plot([0, v[0,1]], [0, v[1,1]], 'g')
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 1.5)
    plt.show()




def pca(data, project=False, verbose=False):
    """Compute the principal components of the given data based on eigenvalue decomposition
    Input:
        data: First index (rows) are variables,
              second index (columns) are observations
        project: Whether to project the input data on the principal components
    Returns:
        (eigvals, eigvecs)
        Eigenvalues are sorted in decreasing order,
        The corresponding eigenvectors are the columns of the second returned array.
        If project is True:
        (eigvals, eigvecs, projectedData)
    """

    # diagonalize covariance matrix
    if verbose:
        print("Constructing covariance matrix:", end="")
        starttime = time.time()
    covariance = np.cov(data)
    if verbose:
        print(" {:.2f} sec.".format(time.time() - starttime))

    if verbose:
        print("Diagonalising covariance matrix:", end="")
        starttime = time.time() 
    eigvals, eigvecs = np.linalg.eig(covariance)
    if verbose:
        print(" {:.2f} sec.".format(time.time() - starttime)) 

    # sort eigenvectors by eigenvalue in decreasing order
    if verbose:
        print("Sorting eigenvectors by eigenvalues:", end="")
        starttime = time.time()  
    sortedeigvecs = np.copy(eigvecs)
    tmp = [[eigvals[i], eigvecs[:,i]] for i in range(len(eigvals))]
    tmp.sort()
    i = len(eigvals)
    for eigval, eigvec in tmp:
        i -= 1
        eigvals      [  i] = eigval
        sortedeigvecs[:,i] = eigvec
    if verbose:
        print(" {:.2f} sec.".format(time.time() - starttime))  

    if project:
        if verbose:
            print("Projecting data in eigenspace:", end="")
            starttime = time.time()   
        # Project data on principal components
        projectedData = np.zeros_like(data)
        for point_idx in range(data.shape[1]):
            for eig_idx in range(eigvecs.shape[0]):
                projectedData[eig_idx, point_idx] = np.dot(eigvecs[:,eig_idx], data[:,point_idx]) 
        if verbose:
            print(" {:.2f} sec.".format(time.time() - starttime))   
        return eigvals, sortedeigvecs, projectedData

    else:
        return eigvals, sortedeigvecs



def dpca(dihedrals, unit='degree', verbose=False):
    """Perform a dihedral pca
    Input:
        Dihedral angle data: rows (first index) are the angles and
                             columns (second index) observations
        unit: degree or radian
    Returns:
        Array of coordinates in the space spanned by the dihedral principal components
    """

    # Create cartesian coordinate space of x = cos(phi), y = sin(phi)
    cartcoords = np.zeros([2 * dihedrals.shape[0], dihedrals.shape[1]], dtype=np.float)
    if unit == "degree":
        cosines = np.cos(const.pi / 180.0 * dihedrals)
        sines   = np.sin(const.pi / 180.0 * dihedrals)
    elif unit == "radian":
        cosines = np.cos(dihedrals)
        sines   = np.sin(dihedrals) 
    else:
        print("Angular unit must be degree or radian but not {}".format(unit))
        sys.exit(0)
    cos_idx = np.arange(0,2*dihedrals.shape[0],2)
    sin_idx = np.arange(1,2*dihedrals.shape[0],2)
    cartcoords[cos_idx,:] = cosines
    cartcoords[sin_idx,:] = sines

    # Compute pca
    eigvals, eigvecs = pca(cartcoords, verbose=verbose)

#    # Project data on principal components
#    if verbose:
#        print("Projecting data on principal components:", end="")
#        starttime = time.time()   
#    projectedcoords = np.zeros_like(cartcoords)
#    for point_idx in range(cartcoords.shape[1]):
#        for eig_idx in range(eigvecs.shape[0]):
#            projectedcoords[eig_idx, point_idx] = np.dot(eigvecs[:,eig_idx], cartcoords[:,point_idx])
#    if verbose:
#        print(" {:.2f} sec.".format(time.time() - starttime))   

    # Project data on principal components more efficiently
    if verbose:
        print("Projecting data on principal components:", end="")
        starttime = time.time()   
    projectedcoords = np.zeros_like(cartcoords)
    for eig_idx in range(eigvecs.shape[0]):
        product = cartcoords * eigvecs[:,eig_idx].reshape([eigvecs.shape[0],1])
        projectedcoords[eig_idx,:] = product.sum(0)
    if verbose:
        print(" {:.2f} sec.".format(time.time() - starttime))    

    return eigvals, projectedcoords
    



def fes(projectedcoords, dim=2, T=300, nbins=5):
    """Estimate the free energy landscape for the given number of dimensions
    Ei = k*T * [ln(Nmax) - ln(Ni)]
    """
    if dim > 3 or dim > projectedcoords.shape[0]:
        print("dim must not be larger than 3 or the dimension of input coordinates""")
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
            #print(point_idx, dim_idx, binIndex[dim_idx])
        counts[tuple(binIndex)] += 1

    # estimate FES
    Nmax     = counts.max()
    lnCounts = np.log(counts+1)
    FES = const.k * T * (np.log(Nmax) - lnCounts)

    return coordSpace, FES


    

def compute_distances(data, i):
    """Compute the distances of point i in data to all other points in data, including i"""
    diff   = data - data[:,i].reshape(data.shape[0], 1)
    diff **= 2
    dist   = diff.sum(0)
    dist **= 0.5
    return dist
    

def J_nearest_neighbours(data, i, J, sortmethod="heapsort"):
    """Return the indices of the J nearest neighbours of datapoint i"""
    dist = compute_distances(data, i)
    #dist = zip(dist, range(len(dist)))
    #dist.sort()
    #neighbours = zip(*dist)[1]
    neighbours = dist.argsort(kind=sortmethod)
    return set(neighbours[1:J+1])


def nearest_neighbours_cutoff(data, i, cutoff):
    """Return indices of all neighbours of datapoint i within cutoff"""
    dist = compute_distances(data, i)
    neighbours = set(np.where(dist < cutoff)[0])
    neighbours = neighbours.difference(set([i]))
    return neighbours


def construct_clusters_nonrecursive(neighbours, K, point):
    """Non recursively construct cluster from neighbourlists"""

    cluster   = set([point])
    me        = point                                 # current point index
    children  = neighbours[point].difference(cluster) # all yet unprocessed children
    parents   = []                                    # stack to hold all ancestors nodes, pop() returns the direct parent

    # process nodes until all children of starting node have been processed
    while (len(parents) > 0 or len(children) > 0):

        # if all children have been processed, go back to parent node
        if len(children) == 0:
            me       = parents.pop()
            children = neighbours[me].difference(cluster)

        # process next child
        elif len(children) > 0:
            child = children.pop()

            # if child is in cluster
            if me in neighbours[child] and len(neighbours[me].intersection(neighbours[child])) >= K and child not in cluster:   
           #if                             len(neighbours[me].intersection(neighbours[child])) >= K and child not in cluster:   
                # add it to cluster
                cluster.update([child]) 

                # make child current node
                parents.append(me)
                me       = child
                children = neighbours[me].difference(cluster) 

    return cluster


def construct_cluster(neighbours, K, point, cluster):
    """Recursively construct cluster from neighbourlist"""

    # add point to cluster
    cluster.update([point])

    for neighbour in neighbours[point]:

        # if two points are mutual neighbours and have K neighbours in common and the neighbour is not yet part of the cluster
        if point in neighbours[neighbour] and len(neighbours[point].intersection(neighbours[neighbour])) >= K and neighbour not in cluster: 

        # if have K neighbours in common and the neighbour is not yet part of the cluster
       #if                                    len(neighbours[point].intersection(neighbours[neighbour])) >= K and neighbour not in cluster: 

            construct_cluster(neighbours, K, neighbour, cluster)


def central_structures(clusters, coord):
    """return the index of the central structure of each cluster"""

    centers = np.zeros(len(clusters), dtype=np.int)

    for c, cluster in enumerate(clusters):

        coordSum = np.zeros([coord.shape[0]], dtype=np.float)
        for structure_idx in cluster:
            coordSum += coord[:,structure_idx]
        clusterMean = coordSum / len(cluster)

        minDist = float('inf')
        for structure_idx in cluster:
            diff = coord[:,structure_idx] - clusterMean
            dist = np.linalg.norm(diff)
            if dist < minDist:
                minDist = dist
                centerStructure = structure_idx

        centers[c] = centerStructure

    return centers


def cluster_jarvis_patrick(data, J=5, K=3, cutoff=None, dim="all", verbose=False):
    """Cluster data (rows = coordinates; columns = observations)
    with the Jarvis-Patrick algorithm.
    J is the number of nearest neighbours to consider,
    K is the number of common neighbours to be part of the same cluster
    dim is the number of dimensions to consider for computing the distance
    """

    # check dim argument validity
    if dim == "all":
        dim   = data.shape[1]
        coord = data
    elif dim <= data.shape[1]:
        coord = data[:dim,:]
    else:
        print("The number of dimensions for clustering must be <= the number of dimensions of the input data")
        sys.exit(1)

    # construct neighbour lists
    starttime = time.time()
    neighbours = []
    neighbourListSize = 0 # number of integers stored
    if cutoff == None:
        for i in range(data.shape[1]):
            neighbours.append(J_nearest_neighbours(data, i, J))
            neighbourListSize += len(neighbours[-1])
            if neighbourListSize * sys.getsizeof(int()) > 4*1024**3:
                print("Neighbour lists are getting too large! Too large cutoff?")
                sys.exit(1)
    else:
        for i in range(data.shape[1]):
            neighbours.append(nearest_neighbours_cutoff(data, i, cutoff))
            neighbourListSize += len(neighbours[-1])
            if neighbourListSize * sys.getsizeof(int()) > 4*1024**3:
                print("Neighbour lists are getting too large! Too large cutoff?")
                sys.exit(1) 
    endtime = time.time()
    if verbose:
        print("Neighbour list construction took {:.2e} sec.".format(endtime - starttime))

    # cluster
    starttime = time.time()
    clusters = []
    remainingPoints = set(range(data.shape[1]))
    while len(remainingPoints) > 0:
        cluster = set()
        point   = remainingPoints.pop()
        #construct_cluster(neighbours, K, point, cluster)
        cluster = construct_clusters_nonrecursive(neighbours, K, point)
        remainingPoints -= cluster
        clusters.append(cluster)
    endtime = time.time()
    if verbose:
        print("Clustering took {:.2e} sec.".format(endtime - starttime) )


    # sort clusters
    starttime = time.time()
    sizes = [len(c) for c in clusters]
    sizesClusters = zip(sizes, clusters)
    sizesClusters.sort()
    sizesClusters.reverse()
    sizes, clusters = zip(*sizesClusters)
    populations = np.array(sizes) / (1.0 * sum(sizes))

    trajectory = np.zeros(data.shape[1], dtype=np.int)
    for i in range(data.shape[1]):
        for c, cluster in enumerate(clusters):
            if i in cluster:
                trajectory[i] = c
                continue
    endtime = time.time()
    if verbose:
        print("Sorting clusters took {:.2e} sec.".format(endtime - starttime) )
 

    centers = central_structures(clusters, coord)

    return populations, clusters, centers, trajectory








def cluster_affinity_propagation(data, memlim=10):
    """Cluster data (rows = coordinates; columns = observations)
    with the affinity propagation algorithm from scikit learn.
    memlim:  memory limit in GB
    """ 
    ndim, nsamples = data.shape
    if nsamples**2*8.0/1024**3 > memlim:
        print("Distance matrix would be larger than {} GB. Aborting.".format(memlim))
        sys.exit(1)

    # compute distance matrix
    dist = np.zeros([nsamples, nsamples])
    for i in range(nsamples):
        dist[:,i] = compute_distances(data, i)
        dist[i,:] = dist[:,i]

    centers, featureTrj = skc.affinity_propagation(dist, copy=False, damping=0.5)
    nclusters = centers.shape[0]
    populations = Counter(featureTrj)
#    print(populations)
#    print(centers)
#    print(featureTrj)

    # sort clusters by size
    popsunsrt = np.array([int(i) for i in populations.values()]) / (1.0 * nsamples)
    indxunsrt = np.array(populations.keys())
    indxsrt   = np.argsort(popsunsrt)[::-1]
    popsrt    = popsunsrt[indxsrt]
    centersrt = centers[indxsrt]

    # renumber featureTrj
    featureTrjsrt = np.zeros_like(featureTrj, dtype=np.int)
    clusters = []
    for c in range(nclusters):
        where = featureTrj == c
        featureTrjsrt[where] = indxsrt[c]
        clusters.append(set(np.where(featureTrjsrt == c)[0]))

#    return popsrt, clusters, centersrt, featureTrjsrt
    return populations, clusters, centers, featureTrj


def plot_clusters(data, featureTrj, centers):
    nclusters = len(centers)
    clustercolors = plt.cm.rainbow(np.linspace(0,1,nclusters))
    fig = plt.figure()
    axs = fig.add_subplot(111)
    for p in range(data.shape[1]):
        axs.plot(data[0,p], data[1,p], 'o', color=clustercolors[featureTrj[p]])
    for i, c in enumerate(centers):
        axs.plot(data[0,c], data[1,c], '*',
                markerfacecolor=clustercolors[i], markersize=25, markeredgewidth=2, label="cluster {}".format(i))
    axs.legend(numpoints=1)
    #axs.set_xlim(0,1)
    #axs.set_ylim(0,1)
    plt.show()



def multiple_2d_gaussians(n=100, g=3):
    """Generate n points from g 2D gaussian distributions"""

    centers = np.random.rand(2,g)
    stdev   = 0.05 + 0.15 * np.random.rand(2,g)

    data   = np.zeros([2,n])
    labels = np.zeros(n, dtype=np.int)
    for i in range(n):
        d = centers[:,i%g].reshape([2,1]) + stdev[:,i%g].reshape([2,1]) * np.random.randn(2,1)
        data[:,i] = d.reshape(2)
        labels[i] = i%g

    np.random.shuffle(data.T)
    return data, labels


