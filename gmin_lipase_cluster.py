#! /usr/bin/env python
#
# GMIN lipase cluster analysis
# Author: Oliver Schillinger
# E-Mail: o.rusche@fz-juelich.de
#
# Usage:
#   gmin_lipase_cluster.py <configFile>
#
# The config file should contain folders of gmin runs

from gmin_lipase_cluster import *
import sys, os, time, math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cPickle as pickle

if __name__ == "__main__":
    main()  


## ============================================================================ #


def main():

    pickleFile = 'gmin_lipase_cluster.pickle'

    what = None
    if len(sys.argv) > 2:
        what = ' '.join(sys.argv[2:])

    directories = read_dirnames()

    # read structures
    print "Reading structures ...", 
    sys.stdout.flush()
    starttime = time.time()
    structures = read_all_structures(directories, pickleFile)
    runtime = time.time() - starttime
    print "{} structures read in {:.2f} seconds".format(len(structures), runtime)

    # compute collective coordinates
    print "Computing collective coordinates ...",
    sys.stdout.flush()
    starttime = time.time()
    for structure in structures:
        structure.comp_collective_coordinates()
    runtime = time.time() - starttime
    print "{:.2f} seconds".format(runtime)

    # cluster
    (structures, lowest) = mean_shift_clustering(structures)

#    structures = lowest

#    for structure in structures[:4]:
#        print structure

    ranking = []
    for structure in structures:
        ranking.append([structure.energy, structure.filename])

    ranking.sort()

    for r in ranking:
        print "{:.2f}".format(r[0]), r[1][80:]


    # print ranked structures filenames and paths to file
    filename = 'ranked_structures'
    outfile = open(filename, 'w')
    for r in ranking:
        outfile.write(r[1] + '.pdb\n')
    outfile.close()

    plot_CitA_collective_coordinates(structures, what)

       
# ============================================================================ #


def read_dirnames():
    """Read directory names of gmin runs
    Directories returned have a slash appended"""

    try:
        infilename = sys.argv[1]
    except IndexError:
        print "Usage: {} <configFile>".format(os.path.basename(__file__))
        sys.exit(1)

    infile = open(infilename, 'r')
    lines = infile.readlines()
    infile.close();

    directories = []

    for i, d in enumerate(lines):

        if d[0] != '#':

            if d[-1] == '\n':
                d = d[:-1]

            d = os.path.realpath(d)

            if d[-1] != '/':
                d = d + '/' 

            if not os.path.isdir(d):
                print "ERROR: Not a directory: {}".format(d)
                sys.exit(1)

            directories.append(d)

    return directories


# ============================================================================ #


def read_structures(directory):
    """Read all GMIN structures from the given directory"""

    structures = []

    # check if file 'lowest' exists
    if not os.path.isfile(directory + "lowest"):
        print "ERROR: {} does not contain file 'lowest'".format(directory)
        sys.exit(1)

    # read file 'lowest'
    lowestfile = open(directory + "lowest", 'r')
    lowest     = lowestfile.readlines()
    lowestfile.close()

    structureIndex = 0

    for line in lowest:

        if "Energy of minimum" in line:

            structureIndex += 1

            words  = line.split()
            energy = float(words[4])

            structureFilename = directory + "dbase." + str(structureIndex)

            structure = Structure(structureFilename)
            
            structure.set_energy(energy)

            structures.append(structure)

    return structures

# ============================================================================ #


def read_all_structures(directories, pickleFileName):
    """Given a list of directories, read all GMIN structures"""

    # if file with pickled data exists
    if os.path.isfile(pickleFileName):

        # read pickled data
        pickleFile = open(pickleFileName, 'r')
        (pickledDirectories, structures) = pickle.load(pickleFile)
        pickleFile.close()

        # check if list of directories is the same
        sameDirectories = True
        for directory in directories:
            if directory not in pickledDirectories:
                sameDirectories = False
        for directory in pickledDirectories:
            if directory not in directories:
                sameDirectories = False

        # if paths to structures are the same,
        # return them
        if sameDirectories:
            print 'Read structures from pickled File'
            return structures


    # if paths to structures are not the same,
    # read them in from respective files

    structures = []

    for directory in directories:

        structures += read_structures(directory)

    # pickle structures
    pickleFile = open(pickleFileName, 'w')
    pickle.dump((directories, structures), pickleFile)
    pickleFile.close() 

    return structures


# ============================================================================ #

def plot_CitA_collective_coordinates(structures, what=None):
    """The function name says it"""

    # Compile coordinates into array
    energy     = np.zeros(len(structures), dtype=np.float) 
    CitA_r     = np.zeros(len(structures), dtype=np.float)
    CitA_theta = np.zeros(len(structures), dtype=np.float)
    CitA_phi   = np.zeros(len(structures), dtype=np.float)
    Lip_r      = np.zeros(len(structures), dtype=np.float)
    Lip_theta  = np.zeros(len(structures), dtype=np.float)
    Lip_phi    = np.zeros(len(structures), dtype=np.float) 

    for i, structure in enumerate(structures):
        energy    [i] = structure.energy
        CitA_r    [i] = structure.CitA_r    
        CitA_theta[i] = structure.CitA_theta / math.pi * 180.0
        CitA_phi  [i] = structure.CitA_phi   / math.pi * 180.0
        Lip_r     [i] = structure.Lip_r     
        Lip_theta [i] = structure.Lip_theta  / math.pi * 180.0
        Lip_phi   [i] = structure.Lip_phi    / math.pi * 180.0

    # scale energies to 0-1 range for coloring
    # (1 = high energy / bad, 0 = low energy / good
    energy_colors = (energy - min(energy)) / max(energy - min(energy))
    #energy_colors = [ str(i) for i in energy_colors ]

    # (0 = high energy / bad, 1 = low energy / good
    energy_colors_invert = -1*(((energy - min(energy)) / max(energy - min(energy))) - 1) 

    # some plot parameters
    alpha       = 0.5
    energyScale = 'energy (kcal/mol)'
    lengthUnit  = 'Angstrom'
    angleUnit   = 'Degree'
    markerSize  = 200
    colorMap    = plt.cm.cool

    energySpan   = max(energy) - min(energy)
    energyMargin = 0.1
    energyRange  = [min(energy)-energyMargin*energySpan, max(energy)+energyMargin*energySpan]
    phiRange     = [-180.0, 180.0]
    #phiRange     = [100.0, 122.0]
    thetaRange   = [0.0, 180.0]
    #thetaRange   = [30.0, 75.0]
    rRange       = [24.0, 38.0]

    # Create plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    # Plot data

    if what is None:
        what = 'CitA phi'

    if what == 'CitA phi':
        sc = ax.scatter(CitA_phi, energy, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA phi (' + angleUnit + ')')
        ax.set_ylabel(energyScale)   
        plt.xlim(phiRange)
        plt.ylim(energyRange)
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale)  

    elif what == 'CitA theta':
        sc = ax.scatter(CitA_theta, energy, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA theta (' + angleUnit + ')')
        ax.set_ylabel(energyScale)   
        plt.xlim(thetaRange)
        plt.ylim(energyRange)
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale)   

    elif what == 'CitA r':
        sc = ax.scatter(CitA_r, energy, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA r (' + lengthUnit + ')')
        ax.set_ylabel(energyScale)   
        plt.xlim(rRange)
        plt.ylim(energyRange)
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale) 

    elif what == 'CitA phi theta':
        sc = ax.scatter(CitA_phi, CitA_theta, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA phi (' + angleUnit + ')')
        ax.set_ylabel('CitA theta (' + angleUnit + ')')
        plt.xlim(phiRange)
        plt.ylim(thetaRange) 
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale)

    elif what == 'CitA phi r':
#        ax = fig.add_subplot(111, polar=True)
#        s  = plt.scatter(CitA_phi, CitA_theta, c=energy_colors, s=10+100*energy_colors, cmap=plt.cm.cool)
#        s.set_alpha(alpha)
#        plt.colorbar()

        sc = ax.scatter(CitA_phi, CitA_r, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA phi (' + angleUnit + ')')
        ax.set_ylabel('CitA r (' + lengthUnit + ')')
        plt.xlim(phiRange)
        plt.ylim(rRange) 
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale) 

    elif what == 'CitA theta r':
        sc = ax.scatter(CitA_theta, CitA_r, c=energy, s=markerSize, cmap=colorMap)
        ax.set_xlabel('CitA theta (' + angleUnit + ')')
        ax.set_ylabel('CitA r (' + lengthUnit + ')')
        plt.xlim(thetaRange)
        plt.ylim(rRange) 
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale)  

    # ============ #

    elif what == 'Lip phi':
        ax.plot(Lip_phi, energy, 'o')
        ax.set_xlabel('Lip phi (rad)')
        ax.set_ylabel(energyScale) 

    elif what == 'Lip theta':
        ax.plot(Lip_theta, energy, 'o')
        ax.set_xlabel('Lip theta (rad)')
        ax.set_ylabel(energyScale)  

    elif what == 'Lip r':
        ax.plot(Lip_r, energy, 'o')
        ax.set_xlabel('Lip r (Angstrom)')
        ax.set_ylabel(energyScale)    

    # ============ #

    elif what == 'ultimate':
        edgecolors = []
        for angle in Lip_theta:
            angle = angle/180.0
            edgecolors.append(colorMap(angle))

        sc = ax.scatter(CitA_phi, CitA_theta, c=energy, s=200+Lip_phi, edgecolors=edgecolors, linewidths=2, cmap=colorMap)
        ax.set_xlabel('CitA phi (' + angleUnit + ')')
        ax.set_ylabel('CitA theta (' + angleUnit + ')')
        plt.xlim(phiRange)
        plt.ylim(thetaRange) 
        sc.set_alpha(alpha)
        cbar = plt.colorbar(sc)
        cbar.set_label(energyScale)   

    # ============ #

    else:
        print "Plot option not known: {}".format(what)
        sys.exit(1) 

    plt.show()

# ============================================================================ #

def mean_shift_clustering(structures):
    """Cluster the structures using mean shift clustering"""

#    structures = structures[:4]

    # parameters
    h = 0.2
    eps = 1e-5

    # Compile coordinates into array
    coordinates = np.zeros([6, len(structures)], dtype=np.float)
    for i, structure in enumerate(structures):
        coordinates[0,i] = structure.energy
        coordinates[1,i] = structure.CitA_r    
        coordinates[2,i] = structure.CitA_theta
        coordinates[3,i] = structure.CitA_phi  
        coordinates[4,i] = structure.Lip_theta 
        coordinates[5,i] = structure.Lip_phi   

#    np.savetxt("coordinates.dat", coordinates)


    (coordinates, stdev, means) = normalize_coordinates(coordinates)


#    fig = plt.figure()
#    ax  = fig.add_subplot(111) 
#    ax.scatter(coordinates[0,:], coordinates[1,:])

    endCoordinates = np.zeros_like(coordinates)

    # loop through points
    for i in range(coordinates.shape[1]):
        p = coordinates[:,i]
        M = m(p, coordinates, h);

        # find fixpoint
        oldNorm = 2*np.linalg.norm(M)
        while abs(np.linalg.norm(M) - oldNorm) > eps:
            p = p + M
            oldNorm = np.linalg.norm(M)
            M = m(p, coordinates, h)

        endCoordinates[:,i] = p

    
    (clusterCoords, representatives, lowest) = find_fix_points(endCoordinates, coordinates, 1e-1)

#    ax.scatter(clusterCoords[0,:], clusterCoords[1,:], c='r', marker='x')
#    ax.scatter(clusterCoords[0,2], clusterCoords[1,2], c='g', marker='^')
#    plt.show()

    coordinates   = denormalize_coordinates(coordinates, stdev, means)
    clusterCoords = denormalize_coordinates(clusterCoords, stdev, means)

    # compile list of representative sturctures (nearest to cluster center)
    repStructures = []
    for r in representatives:
        repStructures.append(structures[r])

    # compile list of representative sturctures (lowest energy in cluster)
    lowestStructures = []
    for l in lowest:
        lowestStructures.append(structures[l]) 

    
    return repStructures, lowestStructures

# ============================================================================ #

def m(x, points, h):
    """Compute negative density gradient"""
    A = np.zeros_like(x)
    B = 0
    for i in range(points.shape[1]):

        p = points[:,i]
        G = g(np.linalg.norm(x-p)**2 / h**2)

        A += p * G
        B +=     G

    import warnings
    warnings.filterwarnings('error')

    try:
        M = A/B - x
    except RuntimeWarning:
        M = -x

    return M

# ============================================================================ #

def g(x):
    """negative derivative of a gaussian"""
    G =  -2 * x * math.exp(-1 * x**2)
    return G

# ============================================================================ #

def normalize_coordinates(coordinates):
    """Normalize coordinates such that each coordinate has
    a stdev of 1 and mean 0.
    Return also a vector of stdevs and means for backtransformation"""
    
#    stdev = np.mat(np.std( coordinates, 1)).T
#    means = np.mat(np.mean(coordinates, 1)).T

    stdev = np.std( coordinates, 1)
    means = np.mean(coordinates, 1)

    normCoordinates = np.copy(coordinates)
    for col in range(normCoordinates.shape[1]):
        normCoordinates[:,col] = (normCoordinates[:,col] - means) / stdev

#    normCoordinates = coordinates - means
#    normCoordinates = normCoordinates / stdev

    return normCoordinates, stdev, means

# ============================================================================ #

def denormalize_coordinates(coordinates, stdev, means):
    """Reverse transformation of normalize_coordinates()"""

    coord = np.copy(coordinates)

    for col in range(coord.shape[1]):
        coord[:,col] = (coord[:,col] * stdev) + means
    
#    coord = coordinates / (1/stdev)
#    coord = coord + means

    return coord

# ============================================================================ #

def find_fix_points(points, originalP, eps):
    """Return a list of unique points in array points
    Two points are considered to be equal when their distance is
    smaller than eps"""

    cluster = -1 * np.ones(points.shape[1], dtype=np.int)
    currentCluster = -1

    # loop through all points
    for i in range(points.shape[1]):
        # if current point has not been assigned to a cluster yet
        if cluster[i] == -1:
            currentCluster += 1
            cluster[i] = currentCluster
            p = points[:,i]
            # find all points belonging to his cluster
            for j in range(i,points.shape[1]):
                if cluster[j] == -1 and np.linalg.norm(p-points[:,j]) <= eps:
                    cluster[j] = currentCluster

    # collect cluster coordinates and indices of representative points
    representatives = np.zeros(max(cluster)+1, dtype=np.int)
    lowest          = np.zeros(max(cluster)+1, dtype=np.int)
    distances       = (max(cluster)+1) * [float('Inf')]
    clusterCoords   = np.zeros([points.shape[0], max(cluster)+1], dtype=np.float)

    currentCluster = 0
    # loop through points
    for i in range(points.shape[1]):
        lowestEnergy = float('Inf');

        # only care about those belonging to the current cluster
        if cluster[i] == currentCluster:
            clusterCoords[:,currentCluster] = points[:,i]

            for j in range(originalP.shape[1]):

                if np.linalg.norm(clusterCoords[:,currentCluster] - originalP[:,j]) < distances[currentCluster]:
                    representatives[currentCluster] = j
                    distances[currentCluster] = np.linalg.norm(clusterCoords[:,currentCluster] - originalP[:,j])

                if cluster[j] == currentCluster and originalP[0,j] < lowestEnergy:
                    lowestEnergy = originalP[0,j]
                    lowest[currentCluster] = j

            currentCluster += 1


#    print representatives
#    print distances
#    
#    print 'Clusters found: ', max(cluster)+1

    return clusterCoords, representatives, lowest

# ============================================================================ #


class Structure:
    """Contains information about a protein structure"""

# ==================================== #

    def __init__(self, filename = None):
        """Initialize object
        Read data from GMIN dbase.x file, which has higher precision than PDB"""

        # general information
        self.filename  = filename
        self.energy    = 0.0
        self.length    = 0
        self.sequence  = []

        # per atom information
        self.resName    = []
        self.resNum     = []
        self.atomName   = []

        # 3*N array
        self.coordinates = [] 

        # Rotation matrices
        self.CitA_R = np.eye(3)
        self.Lip_R  = np.eye(3)

        # Collective coordinates
        self.CitA_r     = 0.0
        self.CitA_theta = 0.0
        self.CitA_phi   = 0.0
        self.Lip_r      = 0.0
        self.Lip_theta  = 0.0
        self.Lip_phi    = 0.0

        if filename != None:
            self.read_data(filename)

# ==================================== #

    def __str__(self):

        string = "Energy = {}\n".format(self.energy)
        string += "Length = {} Aminoacids".format(self.length)

        return string

# ==================================== #

    def show(self, n=None):
        """Print structure information up to atom n"""

        if n == None:
            n = self.length

        print "Energy = {}".format(self.energy)
        print "Length = {} Aminoacids".format(self.length)
        print "Sequence = ",
        for res in self.sequence:
            print res,
        print ""


# ==================================== #

    def read_data(self, filename):
        """Read data from GMIN dbase.x file, which has higher precision than PDB"""

        if not os.path.isfile(filename):
            print "ERROR: File does not exist: {}".format(filename)
            sys.exit(1)

        datafile = open(filename, 'r')
        data     = datafile.readlines()
        datafile.close()

        numberAtoms = int(data[0])

        self.coordinates = np.zeros(3 * numberAtoms, dtype=np.float)

        oldResNum = 0

        for i, line in enumerate(data[1:]):

            line = line.split()

            resNum   = int(line[1])
            resName  =     line[2]
            atomName =     line[3]

            self.resNum.append(resNum)
            self.resName.append(resName)
            self.atomName.append(atomName)

            if resNum > oldResNum:
                oldResNum = resNum
                self.sequence.append(resName)

            x = float(line[4])
            y = float(line[5])
            z = float(line[6])

            self.coordinates[i*3 + 0] = x
            self.coordinates[i*3 + 1] = y
            self.coordinates[i*3 + 2] = z

        self.length = len(self.sequence)

# ==================================== #

    def set_energy(self, energy):
        """set energy for structure"""

        self.energy = energy

# ==================================== #

    def comp_collective_coordinates(self):
        """Compute collective coordinates

        Definitions (atom counts start with 0)
        CitA:   Atom   78 - 2280 (ARG6   - LEU148)
        Lipase: Atom 2428 - 5009 (PRO159 - ASN328)
        
        Algorithm:
        Shift CitA center of mass (COM) to origin (each atom has mass 1).
        Align CitA COM - LYS34:CA (Atom 529) with positive x-axis.
        Rotate CitA PRO71:CA (Atom 1121) around x-axis into yx-plane (positive y-axis).
        Determine spherical coordinates of Lipase COM.
        Shift Lipase COM to origin.
        Determine spherical coordinates of THR199:CA (Atom 3050).
        """

        # First and last atoms of the two sub-molecules
        CitAstart = 78    # ARG6
        CitAend   = 2280  # LEU148
        LipStart  = 2428  # PRO159
        LipEnd    = 5009  # ASN328

        # Atom indices of refernece points for orientation determination
        CitArefi  = 529   # LYS34:CA
        CitAref2i = 1121  # PRO71:CA
        LipRefi   = 3050  # THR199:CA

        # basis vectors
        x = np.array([1,0,0], dtype=np.float)
        y = np.array([0,1,0], dtype=np.float)
        z = np.array([0,0,1], dtype=np.float)

        # Determine COM
        CitACOM = self.COM(CitAstart, CitAend)
        LipCOM  = self.COM(LipStart,  LipEnd )

        # Compute CitA alignment rotation matrix

        # CitA reference points coordinates in local CitA coordinate system
        CitAref  = self.coordinates[CitArefi*3:(CitArefi+1)*3]   - CitACOM 
        CitAref2 = self.coordinates[CitAref2i*3:(CitAref2i+1)*3] - CitACOM 

        # CitA reference point angle for rotation around y-axis into xy-plane
        tmpVec = np.copy(CitAref)
        tmpVec[1]   = 0.0  # project into xz-plane
        CitAyAngle  = math.acos(np.dot(tmpVec, x) / np.linalg.norm(tmpVec)) # angle with x-axis

        # Set up matrix for rotation around y-axis
        if CitAref[2] > 0:
            CitAyAngle = -CitAyAngle
        Ry = np.array([[math.cos(CitAyAngle), 0, -math.sin(CitAyAngle)],
                       [                   0, 1,                     0],
                       [math.sin(CitAyAngle), 0,  math.cos(CitAyAngle)]]) 

        # rotate around y-axis into xy-plane
        yRotated  = np.inner(Ry, CitAref )
        yRotated2 = np.inner(Ry, CitAref2)


        # CitA reference point angle for rotation around z-axis onto positive x-axis
        CitAzAngle  = math.acos(np.dot(yRotated, x) / np.linalg.norm(yRotated)) # angle with x-axis 

        # Set up matrix for rotation around z-axis
        if yRotated[1] < 0:
            CitAzAngle = -CitAzAngle
        Rz = np.array([[ math.cos(CitAzAngle), math.sin(CitAzAngle), 0],
                       [-math.sin(CitAzAngle), math.cos(CitAzAngle), 0], 
                       [                   0,                    0, 1]]) 

        # Rotate around z-axis onto positive x-axis
        zRotated  = np.inner(Rz, yRotated ) 
        zRotated2 = np.inner(Rz, yRotated2) 


        # CitA reference point 2 rotation around x-axis to point into positive y-direction
        tmpVec = np.copy(zRotated2)
        tmpVec[0]   = 0.0  # project into yz-plane
        CitAxAngle  = math.acos(np.dot(tmpVec, y) / np.linalg.norm(tmpVec)) # angle with y-axis 

        # Set up matrix for rotation around x-axis
        if zRotated2[2] < 0:
            CitAxAngle = -CitAxAngle
        Rx = np.array([[1,                     0,                    0],
                       [0,  math.cos(CitAxAngle), math.sin(CitAxAngle)],
                       [0, -math.sin(CitAxAngle), math.cos(CitAxAngle)]]) 

        # Rotate CitA reference point 2 around all 3 axis
        xRotated  = np.inner(Rx, zRotated )
        xRotated2 = np.inner(Rx, zRotated2)

        # Combine all 3 rotation matrices
        Rxz = np.dot(Rx , Rz)
        R   = np.dot(Rxz, Ry)
        self.CitA_R = np.copy(R)

        # Rotate both reference points
        CitArefRotated  = np.inner(R, CitAref)
        CitAref2Rotated = np.inner(R, CitAref2) 


        # Transform Lipase COM to CitA local coordinate system
        LipCOMlocal = LipCOM - CitACOM
        LipCOMlocal = np.inner(R, LipCOMlocal)
        

        # Compute sperical coordinates of Lipase COM in CitA local coordinate system
        self.CitA_r     = np.linalg.norm(LipCOMlocal)
        self.CitA_theta = math.acos(LipCOMlocal[2]/self.CitA_r)
        self.CitA_phi   = math.atan2(LipCOMlocal[1], LipCOMlocal[0])

#        print "{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(self.energy,
#                                                      self.CitA_r,
#                                                      self.CitA_theta/math.pi*180,
#                                                      self.CitA_phi/math.pi*180)


        # Transform Lipase refernece coordinates to CitA local coordinate system
        LipRef = self.coordinates[LipRefi*3:(LipRefi+1)*3] - CitACOM 
        LipRef = np.inner(R, LipRef)

        # Center coordinates on Lipase COM
        LipRefLocal = LipRef - LipCOMlocal

        # Compute sperical coordinates of Lipase reference in lipase local coordinate system
        self.Lip_r     = np.linalg.norm(LipRefLocal)
        self.Lip_theta = math.acos(LipRefLocal[2]/self.Lip_r)
        self.Lip_phi   = math.atan2(LipRefLocal[1], LipCOMlocal[0]) 

#        print "{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(self.energy,
#                                                      self.Lip_r,
#                                                      self.Lip_theta/math.pi*180,
#                                                      self.Lip_phi/math.pi*180) 
        

# ==================================== #

    def COM(self, first, last):
        """Determine center of mass
        All atoms are considered to have the same Mass.
        first and last are the indices of the first and last
        atoms to be considered"""

        COM = np.zeros(3, dtype=np.float)
        for i in range(first, last+1):
            COM += self.coordinates[i*3:(i+1)*3]

        return COM / (last-first+1)

        
# ============================================================================ #

