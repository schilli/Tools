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

    plot_CitA_collective_coordinates(structures)

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

def plot_CitA_collective_coordinates(structures):
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
        CitA_theta[i] = structure.CitA_theta
        CitA_phi  [i] = structure.CitA_phi  
        Lip_r     [i] = structure.Lip_r     
        Lip_theta [i] = structure.Lip_theta 
        Lip_phi   [i] = structure.Lip_phi   

    # scale energies to 0-1 range for coloring
    # (0 = high energy / bad, 1 = low energy / good
    energy_colors = -1*(((energy - min(energy)) / max(energy - min(energy))) - 1)

    # Create plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    # Plot data

    what = 'CitA phi r'

    if what == 'CitA phi':
        ax.plot(CitA_phi, energy, 'o')
        ax.set_xlabel('CitA phi (rad)')
        ax.set_ylabel('energy (a.u.)') 

    elif what == 'CitA theta':
        ax.plot(CitA_theta, energy, 'o')
        ax.set_xlabel('CitA theta (rad)')
        ax.set_ylabel('energy (a.u.)')  

    elif what == 'CitA r':
        ax.plot(CitA_r, energy, 'o')
        ax.set_xlabel('CitA r (a.u.)')
        ax.set_ylabel('energy (a.u.)')   

    elif what == 'CitA phi theta':
        ax.plot(CitA_phi, CitA_theta, 'o')
        ax.set_xlabel('CitA phi (a.u.)')
        ax.set_ylabel('CitA theta (a.u.)')

    elif what == 'CitA phi r':
        ax = fig.add_subplot(111, polar=True)
        s  = plt.scatter(CitA_phi, CitA_theta, c=energy_colors, s=10+100*energy_colors, cmap=plt.cm.cool)
        s.set_alpha(0.75)
        plt.colorbar()

    else:
        print "Plot option not known: {}".format(what)
        sys.exit(1)

    plt.show()

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

