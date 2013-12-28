#! /usr/bin/env python
#
# Author: Oliver Schillinger
# E-Mail: o.rusche@fz-juelich.de
#
#
# The config file should contain folders of gmin runs

from complex_traj_collective_coords import *
import sys, os, time, math
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

if __name__ == "__main__":
    main()  


# ============================================================================ #


def main():

#    filenames = ["/home/oliver/SiSc/Courses/Thesis/iGRASP/liPASe/jubiofs/iGRASP/Complex/folded/citrate_bound/structure1/prod01.pdb",
#                 "/home/oliver/SiSc/Courses/Thesis/iGRASP/liPASe/jubiofs/iGRASP/Complex/folded/citrate_bound/structure2/prod01.pdb",
#                 "/home/oliver/SiSc/Courses/Thesis/iGRASP/liPASe/jubiofs/iGRASP/Complex/folded/citrate_bound/structure4/prod01.pdb",
#                 "/home/oliver/SiSc/Courses/Thesis/iGRASP/liPASe/jubiofs/iGRASP/Complex/folded/citrate_bound/structure6/prod01.pdb"]

    filenames = ["/home/oliver/BigData/citrate_bound/structure1/prod01.pdb",
                 "/home/oliver/BigData/citrate_bound/structure2/prod01.pdb",
                 "/home/oliver/BigData/citrate_bound/structure4/prod01.pdb",
                 "/home/oliver/BigData/citrate_bound/structure6/prod01.pdb",
                 "/home/oliver/BigData/citrate_free/structure1/prod01.pdb",
                 "/home/oliver/BigData/citrate_free/structure2/prod01.pdb",
                 "/home/oliver/BigData/citrate_free/structure4/prod01.pdb",
                 "/home/oliver/BigData/citrate_free/structure6/prod01.pdb"] 

#    filenames = ["/home/oliver/BigData/prod01.pdb"]


    pickleFileName = "complex_traj_collective_coords.pickle"

    if os.path.isfile(pickleFileName):
        start = time.time()
        pickleFile = open(pickleFileName, 'r')
        MDs = pickle.load(pickleFile)
        pickleFile.close() 
        print "Read MDs from pickled file in {:.2f} sec.".format(time.time()-start)

    else:

        MDs = []
        for filename in filenames:

            start = time.time()
            MDs.append(Trajectory(filename))
            MDs[-1].get_collective_coordinates()
            MDs[-1].frames = []
            print "Done with file {}\n in {:.2f} sec.".format(filename, time.time()-start)

        pickleFile = open(pickleFileName, 'w')
        pickle.dump(MDs, pickleFile)
        pickleFile.close()
    
    
    a = 10
    for i in range(len(MDs)):
        print "MDs[{}]".format(i)
        print "\tCitA_theta = [",
        for n in MDs[i].CitA_theta[:a]:
            print "{:.2f}, ".format(n),
        print "{:.2f}]".format(MDs[i].CitA_theta[a])

    i = 0
    plt.plot(MDs[i  ].CitA_theta, MDs[i  ].CitA_phi, 'r')
    plt.plot(MDs[i+4].CitA_theta, MDs[i+4].CitA_phi, 'b')
    plt.show()





# ============================================================================ #

class Trajectory:
    """Holds structures of a whole trajectory"""

# ==================================== #

    def __init__(self, filename = None):
        """Initialize object"""

        # general information
        self.filename  = filename
        self.nframes   = 0
        self.frames    = []

        # collective coordinates
        self.CitA_r     = []
        self.CitA_theta = []
        self.CitA_phi   = []
        self.Lip_r      = []
        self.Lip_theta  = []
        self.Lip_phi    = []

        if (filename != None):
            self.read_trajectory()
 
# ==================================== #

    def plot_collective_coordinates(self):

        plt.plot(self.CitA_theta, self.CitA_phi)

# ==================================== #

    def get_collective_coordinates(self):

        for frame in self.frames:
            self.CitA_r    .append(frame.CitA_r    )
            self.CitA_theta.append(frame.CitA_theta)
            self.CitA_phi  .append(frame.CitA_phi  )
            self.Lip_r     .append(frame.Lip_r     )
            self.Lip_theta .append(frame.Lip_theta )
            self.Lip_phi   .append(frame.Lip_phi   )


# ==================================== #

    def read_trajectory(self):
        "read whole trajectory from pdb file"

        self.count_frames()

        self.frames = []

        for m in range(self.nframes):
            self.frames.append(Structure(self.filename, m))
            self.frames[-1].comp_collective_coordinates()


# ==================================== #

    def count_frames(self):
        """Count frames in pdb file"""

        if not os.path.isfile(self.filename):
            print "ERROR: File does not exist: {}".format(self.filename)
            sys.exit(1)

        datafile = open(self.filename, 'r')
        data     = datafile.readlines()
        datafile.close()

        self.nframes = 0

        for line in data:
            if line[:5] == "MODEL":
                self.nframes += 1


# ============================================================================ #


class Structure:
    """Contains information about a protein structure"""

# ==================================== #

    def __init__(self, filename = None, n = None):
        """Initialize object
        Read data from GMIN dbase.x file, which has higher precision than PDB"""

        # general information
        self.filename  = filename
        self.energy    = 0.0
        self.length    = 0
        self.frame     = 0
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

        if filename != None and n == None:
            self.read_data(filename)
        elif filename != None and n != None:
            self.read_pdb_traj(filename, n)

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


    def read_pdb_traj(self, filename, n):
        """Read data from pdb file that contains multiple models
        Read model n only"""

        if not os.path.isfile(filename):
            print "ERROR: File does not exist: {}".format(filename)
            sys.exit(1)

        datafile = open(filename, 'r')
        data     = datafile.readlines()
        datafile.close()

        # seek model n
        modelCount = -1
        lineNumber = 0
        while (modelCount < n and lineNumber < len(data)):
            if data[lineNumber][:5] == "MODEL":
                modelCount += 1
            lineNumber += 1

        if (lineNumber >= len(data)):
            print "ERROR: Cannot read model {} from File {}, contains only {} models".format(n, filename, modelCount+1)
            sys.exit(1)

        self.frame = modelCount

        modelStart = lineNumber

        # find end of model n
        while (data[lineNumber][:4] == "ATOM"):
            lineNumber += 1
        modelEnd = lineNumber - 1

        numberAtoms = modelEnd - modelStart + 1

        self.coordinates = np.zeros(3 * numberAtoms, dtype=np.float)

        oldResNum = 0


        for i, line in enumerate(data[modelStart:modelEnd+1]):

            line = line.split()

            shift = 0

            try:
                resNum   = int(line[4+shift])
            except ValueError:
                shift = 1
                resNum   = int(line[4+shift])

            resName  =     line[3]
            atomName =     line[1]

            self.resNum.append(resNum)
            self.resName.append(resName)
            self.atomName.append(atomName)

            if resNum > oldResNum:
                oldResNum = resNum
                self.sequence.append(resName)

            x = float(line[5+shift])
            y = float(line[6+shift])
            z = float(line[7+shift])

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
        CitAstart = 26   #78    # ARG6
        CitAend   = 1118 #2280  # LEU148
        LipStart  = 1198 #2428  # PRO159
        LipEnd    = 2485 #5009  # ASN328

        # Atom indices of refernece points for orientation determination
        CitArefi  = 248  # 529   # LYS34:CA
        CitAref2i = 544  # 1121  # PRO71:CA
        LipRefi   = 1514 # 3050  # THR199:CA

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


