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
import sys, os, time
import numpy as np

if __name__ == "__main__":
    main()  


# ============================================================================ #


def main():

    directories = read_dirnames()

    # read structures
    starttime = time.time()
    structures = read_all_structures(directories)
    runtime = time.time() - starttime
    print "Reading {} structures took {} seconds".format(len(structures), runtime)

#    print "Read {} structures".format(len(structures))

#    for structure in structures[:4]:
#        print structure

       
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
    directories = infile.readlines()
    infile.close();

    for i, d in enumerate(directories):

        if d[0] == '#':
            directories.remove(d)

        else:

            if d[-1] == '\n':
                d = d[:-1]

            d = os.path.realpath(d)

            if d[-1] != '/':
                d = d + '/' 

            if not os.path.isdir(d):
                print "ERROR: Not a directory: {}".format(d)
                sys.exit(1)

            directories[i] = d

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


        
# ============================================================================ #

