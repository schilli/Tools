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

if __name__ == "__main__":
    main()  


## ============================================================================ #


def main():

    directories = read_dirnames()

    # read structures
    print "Reading structures ...", 
    sys.stdout.flush()
    starttime = time.time()
    structures = read_all_structures(directories)
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

#    for structure in structures[:4]:
#        print structure

#    ranking = []
#    for structure in structures:
#        ranking.append([structure.energy, structure.filename])
#
#    ranking.sort()
#
#    for r in ranking:
#        print "{:.2f}".format(r[0]), r[1]

       
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


def read_all_structures(directories):
    """Given a list of directories, read all GMIN structures"""

    structures = []

    for directory in directories:

        structures += read_structures(directory)

    return structures


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
        LipRefi   = 3050  # THR199:CA

        # basis vectors
        x = np.array([1,0,0], dtype=np.float)
        y = np.array([0,1,0], dtype=np.float)
        z = np.array([0,0,1], dtype=np.float)

        # Determine COM
        CitACOM = self.COM(CitAstart, CitAend)
        LipCOM  = self.COM(LipStart,  LipEnd )

        # Compute CitA alignment rotation matrix

        # CitA reference point coordinates in local CitA coordinate system
        CitAref = self.coordinates[CitArefi*3:(CitArefi+1)*3] - CitACOM 

        # CitA reference point angle for rotation around x-axis
        tmpVec      = np.copy(CitAref)
        tmpVec[0]   = 0.0  # project into yz-plane
        CitAxAngle  = math.acos(np.dot(tmpVec, y) / np.linalg.norm(tmpVec)) # angle with y axis

        # Set up matrix for rotation around x-axis
        if CitAref[2] < 0:
            CitAxAngle = -CitAxAngle
        Rx = np.array([[1,                     0,                    0],
                       [0,  math.cos(CitAxAngle), math.sin(CitAxAngle)],
                       [0, -math.sin(CitAxAngle), math.cos(CitAxAngle)]])

        # Rotate around x-axis
        xRotated = np.inner(Rx, CitAref)

        # Determine angle with positive x-axis
        rotatedxAngle = math.acos(np.dot(xRotated, x) / np.linalg.norm(xRotated))

        # Set up rotation matrix around z-axis
        if xRotated[1] < 0:
            rotatedxAngle = -rotatedxAngle
        Rz = np.array([[ math.cos(rotatedxAngle), math.sin(rotatedxAngle), 0],
                       [-math.sin(rotatedxAngle), math.cos(rotatedxAngle), 0], 
                       [                       0,                       0, 1]])

        # Rotate around z-axis
        xzRotated = np.inner(Rz, xRotated)


        # Transform Lipase COM to CitA local coordinate system
        LipCOMlocal = LipCOM - CitACOM
        LipCOMlocal = np.inner(Rx, LipCOMlocal)
        LipCOMlocal = np.inner(Rz, LipCOMlocal)
        

        # Compute sperical coordinates of Lipase COM in CitA local coordinate system
        r     = np.linalg.norm(LipCOMlocal)
        theta = math.acos(LipCOMlocal[2]/r)
        phi   = math.atan2(LipCOMlocal[1], LipCOMlocal[0])

        print r, theta, phi

#        tmpVec      = CitAref / np.linalg.norm(CitAref)
#        print CitAxAngle, tmpVec
        

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

