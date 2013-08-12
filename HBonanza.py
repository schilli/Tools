# HBonanza is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me,
# Jacob Durrant, at jdurrant [at] ucsd [dot] edu.
#
# If you use HBonanza in your work, please cite [REFERENCE HERE]

import math
import sys
import operator
import os
from multiprocessing import Process, Queue, Value, Lock, cpu_count
import platform

def program_info(): # still need to update this
    global program_output
    
    program_output = program_output + "# HBonanza Version 1.0\n"
    program_output = program_output + "# =================\n\n"
    program_output = program_output + "# If you use HBonanza in your research, please cite the following reference:\n"
    program_output = program_output + "# Durrant, J. D., ....\n\n"


# a geometric point
class point:
    x=99999.0
    y=99999.0
    z=99999.0
    
    def __init__ (self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z

    def print_coors(self):
        print str(self.x)+"\t"+str(self.y)+"\t"+str(self.z)
        
    def dist_to(self,apoint):
        return math.sqrt(math.pow(self.x - apoint.x,2) + math.pow(self.y - apoint.y,2) + math.pow(self.z - apoint.z,2))

    def CopyOf(self):
        return point(self.x, self.y, self.z)

    def average_with(self, other_point):
	return point((self.x + other_point.x) / 2.0, (self.y + other_point.y) / 2.0, (self.z + other_point.z) / 2.0)
	
    def dot_product_with(self, other_point):
	return self.x * other_point.x + self.y * other_point.y + self.z * other_point.z
    
    def cross_product_with(self, other_point):
	return point(self.y*other_point.z - self.z*other_point.y, self.z*other_point.x - self.x*other_point.z, self.x*other_point.y - self.y*other_point.x)
    
    def length(self):
	return self.dist_to(point(0.0,0.0,0.0))
	
    def minus(self, other_point):
	return point(self.x - other_point.x, self.y - other_point.y, self.z - other_point.z)
	
    def plus(self, other_point):
	return point(self.x + other_point.x, self.y + other_point.y, self.z + other_point.z)
    
    def scale(self, scalar):
	return point(self.x*scale,self.y*scale,self.z*scale)

# an atom object
class atom:
        
    def __init__ (self):
        self.atomname = ""
        self.resid = 0
        self.chain = ""
        self.resname = ""
        self.coordinates = point(99999, 99999, 99999)
        self.undo_coordinates = point(99999, 99999, 99999)
        self.element = ""
        self.PDBIndex = ""
        self.line=""
        self.IndicesOfAtomsConnecting=[]
	self.beta = 0.00
	self.occupancy = 0.00

    # function to determine if the atom belongs to a protein
    # returns true or false
    def belongs_to_protein(self):
        protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        if self.resname.strip() in protein_residues: return True
	else: return False

    def CopyOf(self):
        newatom = atom()
	newatom.chain = self.chain
	newatom.resid = self.resid
        newatom.atomname = self.atomname
        newatom.resname = self.resname
        newatom.coordinates = self.coordinates.CopyOf()
        newatom.undo_coordinates = self.undo_coordinates.CopyOf()
        newatom.element = self.element
        newatom.PDBIndex = self.PDBIndex
        newatom.line = self.line
        for index in self.IndicesOfAtomsConnecting:
            newatom.IndicesOfAtomsConnecting.append(index)
        return newatom

    # Reads text (PDB format) into an atom object
    # Requires: A string containing the PDB line
    def ReadPDBLine(self, Line):
        self.line = Line
        self.atomname = Line[11:16].strip()
        self.chain = Line[21:22]
        self.resid = int(Line[22:26])
        
        if len(self.atomname)==1: # redo using rjust
            self.atomname = self.atomname + "  "
        elif len(self.atomname)==2:
            self.atomname = self.atomname + " "
        elif len(self.atomname)==3:
            self.atomname = self.atomname + " " # This line is necessary for babel to work, though many PDBs in the PDB would have this line commented out
        
        self.coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))
        
        if len(Line) >= 79: self.element = Line[76:79].strip().upper() # element specified explicitly at end of life
	if len(Line) < 79 or self.element == "":
	    #elif self.element == "": # try to guess at element from name
	    
	    two_letters = self.atomname.strip().upper()[0:2]
            
	    # Any number needs to be removed from the element name
	    two_letters = two_letters.replace('0','')
	    two_letters = two_letters.replace('1','')
	    two_letters = two_letters.replace('2','')
	    two_letters = two_letters.replace('3','')
	    two_letters = two_letters.replace('4','')
	    two_letters = two_letters.replace('5','')
	    two_letters = two_letters.replace('6','')
	    two_letters = two_letters.replace('7','')
	    two_letters = two_letters.replace('8','')
	    two_letters = two_letters.replace('9','')

            if two_letters=='BR':
                self.element='BR'
            elif two_letters=='CL':
                self.element='CL'
            elif two_letters=='BI':
                self.element='BI'
            elif two_letters=='AS':
                self.element='AS'
            elif two_letters=='AG':
                self.element='AG'
            elif two_letters=='LI':
                self.element='LI'
            #elif two_letters=='HG':
            #    self.element='HG'
            elif two_letters=='MG':
                self.element='MG'
            elif two_letters=='RH':
                self.element='RH'
            elif two_letters=='ZN':
                self.element='ZN'
            else: #So, just assume it's the first letter.
		self.element = self.atomname.strip().upper()
		
		# Any number needs to be removed from the element name
		self.element = self.element.replace('0','')
		self.element = self.element.replace('1','')
		self.element = self.element.replace('2','')
		self.element = self.element.replace('3','')
		self.element = self.element.replace('4','')
		self.element = self.element.replace('5','')
		self.element = self.element.replace('6','')
		self.element = self.element.replace('7','')
		self.element = self.element.replace('8','')
		self.element = self.element.replace('9','')
		
                self.element = self.element[0:1]

        self.PDBIndex = Line[6:12].strip()
        self.resname = Line[16:20]
        if self.resname.strip() == "": self.resname = " MOL"

    # Creates a PDB line from the atom object
    # Returns: PDB String
    def CreatePDBLine(self):

        #if len(self.atomname) > 1: self.atomname = self.atomname[:1].upper() + self.atomname[1:].lower()

        output = "ATOM "
        #output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
        output = output + self.PDBIndex.rjust(6) + self.atomname.rjust(5) + self.resname.rjust(4)
	output = output + self.chain.rjust(2)
	output = output + str(self.resid).rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(12)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + ("%.2f" % self.occupancy).rjust(6) # assumes occupancy <=9.99, positive
	output = output + ("%.2f" % self.beta).rjust(6) # assumes beta <=9.99, positive
	
        output = output + self.element.rjust(12) # + "   " + str(uniqueID) #This last part must be removed
        return output
    
    def AddNeighborAtomIndex(self, index):
        if not (index in self.IndicesOfAtomsConnecting):
            self.IndicesOfAtomsConnecting.append(index)
	    
    def NumberOfNeighbors(self):
        return len(self.IndicesOfAtomsConnecting)

    # Sets the undo point for later undoing
    def SetUndoPoint(self):
        self.undo_coordinates = self.coordinates.CopyOf()
        
    # Resets coordinate values after translations or rotations ("Undo")
    def Undo(self):
        self.coordinates = self.undo_coordinates.CopyOf()
    
    # Should be sufficient to say that equals if coordinates, number of neighbors, and element are the same (for lig)
    def __eq__(self,other_atom):
	return self.coordinates.x == other_atom.coordinates.x & self.coordinates.y == other_atom.coordinates.y & self.coordinates.z == other_atom.coordinats.z & self.element == other_atom.element & self.NumberOfNeighbors() == other_atom.NumberOfNeighbors()
    
    def __ne__(self,other_atom):
	return not self.__eq__(other_atom)
   
# PDB class
class PDB:

    def __init__ (self):
        self.AllAtoms={}
	self.N_O_F_S = None
	self.hydrogen_indicies = None
	self.resids = []

    # Loads a PDB from a file
    # Requires: trajectory_filename, a string containing the trajectory_filename
    def LoadPDB_FromFile(self, trajectory_filename):
        # Now load the file into a list
        file = open(trajectory_filename,"r")
        lines = file.readlines()
        file.close()
	self.LoadPDB_FromArray(lines)
    
    def LoadPDB_FromArray(self, lines):

        autoindex = 1

        self.__init__()
        
        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    TempAtom = atom()
                    TempAtom.ReadPDBLine(line)
		    
		    if not TempAtom.resid in self.resids: self.resids.append(TempAtom.resid)
                    
                    self.AllAtoms[autoindex] = TempAtom # because points files have no indices
                    autoindex = autoindex + 1

    # Saves the PDB object to a PDB file
    # Requires: trajectory_filename to be saved
    def SavePDB(self, trajectory_filename):

        if len(self.AllAtoms) > 0: # so the pdb is not empty (if it is empty, don't save)

            file = open(trajectory_filename,"w")

            # write coordinates
            for atomindex in self.AllAtoms:
                file.write(self.AllAtoms[atomindex].CreatePDBLine() + "\n")

            file.close()

    # identifies the greatest distance between any two atoms of the PDB
    # Used for quickly comparing two PDBs so RMSD alignment not necessary if they're different
    # Returns float = distance
    def farthest_away_atoms_dist(self):
	dist_max = 0.0
	keys = self.AllAtoms.keys()
	for key_index1 in range(0,len(keys)-1):
	    key1 = keys[key_index1]
	    atom1 = self.AllAtoms[key1]
	    for key_index2 in range(key_index1+1,len(keys)):
		key2 = keys[key_index2]
		atom2 = self.AllAtoms[key2]
		dist = atom1.coordinates.dist_to(atom2.coordinates)
		if dist > dist_max: dist_max = dist
	self.max_inter_atom_distance = dist_max
	return dist_max

    # Print out info about the PDB
    def print_out_info(self):
        for index in self.AllAtoms:
            print self.AllAtoms[index].CreatePDBLine()

    def NumberOfNeighorsOfElement(self, index, the_element):
        num = 0
        for index in self.AllAtoms[index].IndicesOfAtomsConnecting:
            if self.AllAtoms[index].element == the_element: num = num + 1
        return num
    
    def IndexOfNeighborOfElement(self,index,the_element): # returns the index of the first neighbor of specified element
        for index in self.AllAtoms[index].IndicesOfAtomsConnecting:
            if self.AllAtoms[index].element == the_element: return index
        return -1 # returns -1 if no match
    
    def weight(self):
        total = 0.0
        for atomindex in self.AllAtoms:
            atom = self.AllAtoms[atomindex].element
            mass = 0.0
            if atom == 'H': mass = 1.00794
            elif atom == 'C': mass = 12.0107
            elif atom == 'CL': mass = 35.453
            elif atom == 'N': mass = 14.0067
            elif atom == 'O': mass = 15.9994
            elif atom == 'P': mass = 30.973762
            elif atom == 'S': mass = 32.065
            elif atom == 'BR': mass = 79.904
            elif atom == 'I': mass = 126.90447
            elif atom == 'F': mass = 18.9984032
            elif atom == 'B': mass = 24.3051
            elif atom == 'HG': mass = 200.59
            elif atom == 'BI': mass = 208.98040
            elif atom == 'AS': mass = 74.92160
            elif atom == 'AG': mass = 107.8682
            elif atom == 'K': mass = 39.0983
            elif atom == 'LI': mass = 6.941
            elif atom == 'MG': mass = 24.3050
            elif atom == 'RH': mass = 102.90550
            elif atom == 'ZN': mass = 65.38
            total += mass
        return total
    
    def id_hydrogen_bonds(self, arguments):
	
	# get possible donors/acceptors.
	N_O_F_S = []
	if self.N_O_F_S is None: # so not predefined
	    for index in self.AllAtoms:
		atom = self.AllAtoms[index]
		if atom.element == "O" or atom.element == "N" or atom.element == "F" or atom.element == "S":
		    N_O_F_S.append(index)
	    self.N_O_F_S = N_O_F_S
	else: # so it is predefined
	    N_O_F_S = self.N_O_F_S
		
	# get hydrogen atoms  ### This will be the same for all frames. Why not calculate once?
	hydrogen_indicies = []
	if self.hydrogen_indicies is None: # so not predefined
	    for index in self.AllAtoms:
		atom = self.AllAtoms[index]
		if atom.element == "H":
		    hydrogen_indicies.append(index)
	    self.hydrogen_indicies = hydrogen_indicies
	else: # so is predefined
	    hydrogen_indicies = self.hydrogen_indicies
		
	# now do pairwise comparison to see if any are close to each other
	possible_bonds = []
	for index1 in range(len(N_O_F_S)-1):
	    atom_index1 = N_O_F_S[index1]
	    #if atom_index1 not in self.AllAtoms:
		#print index1, atom_index1
		#print self.AllAtoms.keys()
		#self.SavePDB('test.pdb')
	    #print len(self.AllAtoms.keys())
	    atom1 = self.AllAtoms[atom_index1]
	    for index2 in range(index1 + 1, len(N_O_F_S)):
		atom_index2 = N_O_F_S[index2]
		atom2 = self.AllAtoms[atom_index2]
		dist = atom1.coordinates.dist_to(atom2.coordinates)
		if dist < arguments['HYDROGEN_BOND_DISTANCE_CUTOFF']:
		    possible_bonds.append((atom_index1, atom_index2))
	
	# now check eack of the possible hydrogen bonds, and see if there's a hydrogen in between them.
	possible_bonds_refined = []
	for possible_bond in possible_bonds:
	    heavy1_index = possible_bond[0]
	    heavy2_index = possible_bond[1]
	    heavy1 = self.AllAtoms[heavy1_index]
	    heavy2 = self.AllAtoms[heavy2_index]
	    for hydrogen_index in hydrogen_indicies:
		hydrogen = self.AllAtoms[hydrogen_index]
		dist1 = heavy1.coordinates.dist_to(hydrogen.coordinates)
		if dist1 < arguments['HYDROGEN_BOND_DISTANCE_CUTOFF']:
		    dist2 = heavy2.coordinates.dist_to(hydrogen.coordinates)
		    if dist2 < arguments['HYDROGEN_BOND_DISTANCE_CUTOFF']:
			if dist1 < 1.3 or dist2 < 1.3: # so it's bonded to at least one
			    if dist1 < 1.3: # so the order is H-HD-HA
				possible_bonds_refined.append((hydrogen_index, heavy1_index, heavy2_index))
			    else: 
				possible_bonds_refined.append((hydrogen_index, heavy2_index, heavy1_index))
	# now check the angles
	hydrogen_bonds = []
	for possible_bond_refined in possible_bonds_refined:
	    hydrogen_index = possible_bond_refined[0]
	    bond_donor_index = possible_bond_refined[1]
	    bond_acceptor_index = possible_bond_refined[2]
	    
	    bond_donor = self.AllAtoms[bond_donor_index]
	    hydrogen = self.AllAtoms[hydrogen_index]
	    bond_acceptor = self.AllAtoms[bond_acceptor_index]
	    
	    vector1 = point(hydrogen.coordinates.x - bond_donor.coordinates.x, hydrogen.coordinates.y - bond_donor.coordinates.y, hydrogen.coordinates.z - bond_donor.coordinates.z)
	    vector2 = point(bond_acceptor.coordinates.x - bond_donor.coordinates.x, bond_acceptor.coordinates.y - bond_donor.coordinates.y, bond_acceptor.coordinates.z - bond_donor.coordinates.z)
	    
	    angle = math.acos(vector1.dot_product_with(vector2)/(vector1.length() * vector2.length())) * 180 / math.pi # at most 180
	    
	    if angle <= arguments['HYDROGEN_BOND_ANGLE_CUTOFF']: hydrogen_bonds.append(possible_bond_refined)
	
	return hydrogen_bonds

# a funtion that reads the command-line parameters and stores them in a dictionary called "arguments"
def get_commandline_parameters():
    
    global program_output, extra_file_string
    
    arguments = {}
    
    # set defaults
    arguments['TRAJECTORY_FILENAME'] = ''
    arguments['OUTPUT_BASENAME'] = ''
    arguments['SINGLE_FRAME_FILENAME'] = ''
    arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] = 0.75
    arguments['SEED_RESIDUES'] = []
    arguments['WRITE_COLUMN'] = 'O'
    arguments['HYDROGEN_BOND_DISTANCE_CUTOFF'] = 3.5
    arguments['HYDROGEN_BOND_ANGLE_CUTOFF'] = 30
    
    arguments['DASHED_BONDS'] = "TRUE"
    
    arguments['LOW_FREQUENCY_COLOR_R'] = 255
    arguments['LOW_FREQUENCY_COLOR_G'] = 255
    arguments['LOW_FREQUENCY_COLOR_B'] = 255
    
    arguments['HIGH_FREQUENCY_COLOR_R'] = 0
    arguments['HIGH_FREQUENCY_COLOR_G'] = 255
    arguments['HIGH_FREQUENCY_COLOR_B'] = 0
    
    arguments['NUMBER_PROCESSORS'] = cpu_count()

    arguments['JUST_IMMEDIATE_CONNECTIONS'] = "FALSE" # an option to not identify hydrogen bonds recursively, but to only look at immediate connections to the seed residues.
    
    # get commandline parameters
    for t in range(1,len(sys.argv)):
	sys.argv[t] = sys.argv[t].replace('--','-')
	if sys.argv[t].upper() == '-HELP':
	    program_output = ''
	    program_info()
	    print program_output.rstrip() + "\n"
	    print "INTRODUCTION"
	    print "============\n"
	    print "XXXX is a computer program that can be used to perform hydrogen-bond analyses of molecular dynamics trajectories. It accepts a PDB file of the trajectory as input, and outputs a TCL visualization script that can be loaded into VMD using either the Tk Console or the -e option from the command line. Trajectory files in other formats can be converted to the PDB format easily using VMD's \"Save Coordinates...\" option. The generated TCL file also contains comments to aid the user if he or she wishes to use other molecular visualization software.\n"
	    
	    print "COMAND-LINE PARAMETERS"
	    print "======================\n"
	    print "A number of command-line options are available:\n"
	    
	    print "\t-help\t\t\t\t\tDisplays this help documentation."
	    print "\t-trajectory_filename\t\t\tThe filename containing the trajectory, in PDB format."
	    print "\t-output_basename\t\t\tThe initial string used to begin the names of all output\n\t\t\t\t\t\tfiles."
	    print "\t-single_frame_filename\t\t\tThe filename of a PDB containing a single frame, onto\n\t\t\t\t\t\twhich hydrogen bond visualization will be mapped."
	    print "\t-number_processors\t\t\tThe number of processors to use when performing the\n\t\t\t\t\t\thydrogen-bond analysis of the trajectory. By default, \n\t\t\t\t\t\tall processors on the system are used."
	    print "\t-hydrogen_bond_frequency_cutoff\t\tIn creating visualization files, discard all hydrogen\n\t\t\t\t\t\tbonds that are formed less frequently than this value."
	    print "\t-hydrogen_bond_distance_cutoff\t\tConsider only hydrogen bonds whose donor-acceptor\n\t\t\t\t\t\tdistances are less than or equal to this value."
	    print "\t-hydrogen_bond_angle_cutoff\t\tConsider only hydrogen bonds whose hydrogen-donor-\n\t\t\t\t\t\tacceptor angles are less than or equal to this value."
	    print "\t-dashed_bonds\t\t\t\tIf true (default), hydrogen bonds are represented by\n\t\t\t\t\t\tdashed lines. If false, solid lines are used."
	    print "\t-just_immediate_connections\t\tRather than branching out recursively to identify the\n\t\t\t\t\t\tentire hydrogen-bond network that supports the seed \n\t\t\t\t\t\tresidues, consider only those residues that are\n\t\t\t\t\t\timmediately connected to the seeds via hydrogen bonds."
	    print "\t-seed_residues or -seed_residue\t\tConsider only hydrogen bonds connected to residues with\n\t\t\t\t\t\tthese IDs, or connecting residues that connect to these\n\t\t\t\t\t\tresidues through hydrogen bonds. Multiple seed residues\n\t\t\t\t\t\tcan be specified by evoking this tag multiple times (see\n\t\t\t\t\t\texample below). If no seed residues are specified, the\n\t\t\t\t\t\tprogram will visualize all hydrogen bonds, regardless of\n\t\t\t\t\t\tconnectivity."
	    print "\t-write_column\t\t\t\tThe program writes a PDB file containing the hydrogen-\n\t\t\t\t\t\tbond frequencies in either the occupancy or beta columns.\n\t\t\t\t\t\tIf this tag is set to \"O\", the frequencies are written to\n\t\t\t\t\t\tthe occupancy column. If set to \"B\", they are written to\n\t\t\t\t\t\tthe beta column."
	    print "\t-low_frequency_color_r\t\t\tHydrogen bonds are colored according to frequency. Colors\n\t\t\t\t\t\tare specified using RGB (red, green, blue) values. This\n\t\t\t\t\t\ttag specifies the R value of the color associated with\n\t\t\t\t\t\tthe lowest frequency displayed (given by\n\t\t\t\t\t\t-hydrogen_bond_frequency_cutoff). Valid values range from\n\t\t\t\t\t\t0 to 255."
	    print "\t-low_frequency_color_g\t\t\tThe G value of the color associated with the lowest\n\t\t\t\t\t\tfrequency displayed. Valid values range from 0 to 255."
	    print "\t-low_frequency_color_b\t\t\tThe B value of the color associated with the lowest\n\t\t\t\t\t\tfrequency displayed. Valid values range from 0 to 255."
	    print "\t-high_frequency_color_r\t\t\tThe R value of the color associated with 100% persistent\n\t\t\t\t\t\thydrogen bonds. Valid values range from 0 to 255."
	    print "\t-high_frequency_color_g\t\t\tThe G value of the color associated with 100% persistent\n\t\t\t\t\t\thydrogen bonds. Valid values range from 0 to 255."
	    print "\t-high_frequency_color_b\t\t\tThe B value of the color associated with 100% persistent\n\t\t\t\t\t\thydrogen bonds. Valid values range from 0 to 255.\n"
	    print "DESCRIPTION OF OUTPUT FILES"
	    print "===========================\n"
	    print "Assuming the command-line parameter -output_basename is set to \"output.\", the following files will be created:\n"
	    print "\toutput.average_hbonds\t\t\t\tA two-column text file describing all hydrogen bonds in the trajectory,\n\t\t\t\t\t\t\tregardless of the -hydrogen_bond_frequency_cutoff tag. The first column\n\t\t\t\t\t\t\tis a triplet indicating the atom indices of the atoms that form each\n\t\t\t\t\t\t\thydrogen bond, and the second column is the frequency with which that\n\t\t\t\t\t\t\tparticular hydrogen bonds occurs across all the frames of the trajectory."
	    print "\toutput.frame_by_frame_hbonds.csv\t\tA comma-separated-values file (CSV file) indicating which hydrogen bonds\n\t\t\t\t\t\t\tare present in each of the frames of the trajectory. CSV files can be\n\t\t\t\t\t\t\tloaded by a number of popular spreadsheet programs, including Microsoft\n\t\t\t\t\t\t\tExcel and Open Office Calc."
	    print "\toutput.hbond_averages_in_occupancy_column.pdb\tA pdb file identical to the file specified by the -single_frame_filename"
	    print "\t\t\tOR\t\t\t\ttag, except that the average frequencies with which relevant hydrogen-"
	    print "\toutput.hbond_averages_in_beta_column.pdb\tbond atoms appear in hydrogen bonds across the trajectory are written in\n\t\t\t\t\t\t\tthe occupancy or beta column, depending on the value given by the\n\t\t\t\t\t\t\t-write_column tag. This file may facilitate visualization in programs\n\t\t\t\t\t\t\tother than VMD, if desired.\n"
	    print "EXAMPLES OF USAGE"
	    print "=================\n"
	    print "Example 1: Load in a trajectory named \"traj.pdb\". Define hydrogen bonds to be those that 1) have donor-acceptor distances less than or equal to 3.0 angstroms, 2) have hydrogen-donor-acceptor angles less than or equal to 30 degrees, and 3) are formed in at least 50% of the trajectory frames. Show the visualization based on the coordinates in the single-frame PDB file \"first_frame.pdb\". Write output files using filenames that start with \"output.\" One of these output filenames will be a PDB file with the hydrogen-bond frequencies listed in the occupancy column. Because no seed residues are specified, all hydrogen bonds in the trajectory will be visualized. Hydrogen bonds that appear in 50% of the trajectory frames will appear white, and those that appear in 100% of the trajectories frames will be green. The program output will be directed to the file \"visualize.tcl\" (in UNIX-based systems), which can then be visualized in VMD.\n"
	    print "python HBonanza.py -trajectory_filename traj.pdb -hydrogen_bond_distance_cutoff 3.0 -hydrogen_bond_angle_cutoff 30 -hydrogen_bond_frequency_cutoff 0.5 -single_frame_filename first_frame.pdb -output_basename output. -write_column O -low_frequency_color_r 255 -low_frequency_color_g 255 -low_frequency_color_b 255 -high_frequency_color_r 0 -high_frequency_color_g 255 -high_frequency_color_b 0 > visualize.tcl\n"
	    print "/path/to/VMD/executable/VMD -e visualize.tcl\n"
	    print "Example 2: Same as above, except only hydrogen bonds connected to residues with IDs 156 and 234, or connecting residues that connect to these residues through hydrogen bonds, will be visualized.\n"
	    print "python HBonanza.py -trajectory_filename traj.pdb -hydrogen_bond_distance_cutoff 3.0 -hydrogen_bond_angle_cutoff 30 -hydrogen_bond_frequency_cutoff 0.5 -single_frame_filename first_frame.pdb -output_basename output. -write_column O -low_frequency_color_r 255 -low_frequency_color_g 255 -low_frequency_color_b 255 -high_frequency_color_r 0 -high_frequency_color_g 255 -high_frequency_color_b 0 -seed_residue 156 -seed_residues 234 > visualize.tcl\n"
	    print "Example 3: If the -single_frame tag is not specified, the program will automatically use the first frame from the trajectory specified by the -trajectory_filename tag.\n"
	    print "python HBonanza.py -trajectory_filename traj.pdb -hydrogen_bond_distance_cutoff 3.0 -hydrogen_bond_angle_cutoff 30 -hydrogen_bond_frequency_cutoff 0.5 -output_basename output. -write_column O -low_frequency_color_r 255 -low_frequency_color_g 255 -low_frequency_color_b 255 -high_frequency_color_r 0 -high_frequency_color_g 255 -high_frequency_color_b 0 -seed_residue 156 -seed_residues 234 > visualize.tcl\n"
	    print "Example 4: If an analysis has been performed previously, it's not necessary to specify the -single_frame and -trajectory_filename tags, as long as the same -output_basename and -write_column tags are used. The trajectory hydrogen-bond information will be read from the previously created file \"{output_basename}.average_hbonds\", and the -single_frame tag will be automatically set to the previously created \"{output_basename}.hbond_averages_in_{write_column}*_column.pdb\" file.\n"
	    print "python HBonanza.py -hydrogen_bond_distance_cutoff 3.0 -hydrogen_bond_angle_cutoff 30 -hydrogen_bond_frequency_cutoff 0.5 -output_basename output. -write_column O -low_frequency_color_r 255 -low_frequency_color_g 255 -low_frequency_color_b 255 -high_frequency_color_r 0 -high_frequency_color_g 255 -high_frequency_color_b 0 -seed_residue 156 -seed_residues 234 > visualize.tcl\n"
	    
	    sys.exit()
	if sys.argv[t].upper() == '-TRAJECTORY_FILENAME':
	    arguments[sys.argv[t].upper().replace('-','')] = sys.argv[t+1]
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-JUST_IMMEDIATE_CONNECTIONS':
	    val = sys.argv[t+1].upper()
	    if val != "TRUE" and val !='FALSE':
		program_output = program_output + '# ERROR: The command-line parameter "-just_immediate_connections" must be either "true" or "false," not "' + sys.argv[t+1] + '." Setting to "false" by default.' + "\n\n"
		val = "FALSE"
	    arguments[sys.argv[t].upper().replace('-','')] = val
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-DASHED_BONDS':
	    val = sys.argv[t+1].upper()
	    if val != "TRUE" and val !='FALSE':
		program_output = program_output + '# ERROR: The command-line parameter "-dashed_bonds" must be either "true" or "false," not "' + sys.argv[t+1] + '." Setting to "true" by default.' + "\n\n"
		val = "TRUE"
	    arguments[sys.argv[t].upper().replace('-','')] = val
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-OUTPUT_BASENAME':
	    arguments[sys.argv[t].upper().replace('-','')] = sys.argv[t+1]
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-SINGLE_FRAME_FILENAME':
	    arguments[sys.argv[t].upper().replace('-','')] = sys.argv[t+1]
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HYDROGEN_BOND_FREQUENCY_CUTOFF':
	    arguments[sys.argv[t].upper().replace('-','')] = float(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HYDROGEN_BOND_DISTANCE_CUTOFF':
	    arguments[sys.argv[t].upper().replace('-','')] = float(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HYDROGEN_BOND_ANGLE_CUTOFF':
	    arguments[sys.argv[t].upper().replace('-','')] = float(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-SEED_RESIDUES' or sys.argv[t].upper() == '-SEED_RESIDUE' or sys.argv[t].upper() == '-SEED' or sys.argv[t].upper() == '-SEEDS':
	    arguments['SEED_RESIDUES'].append(int(sys.argv[t+1]))
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-WRITE_COLUMN':
	    if sys.argv[t+1] == "O" or sys.argv[t+1] == "B":
		arguments[sys.argv[t].upper().replace('-','')] = sys.argv[t+1]
		sys.argv[t] = ''
		sys.argv[t+1] = ''
	    else:
		program_output = program_output + '# ERROR: The command-line parameter "-write_column" must be either "O" or "B".' + "\n"
	if sys.argv[t].upper() == '-LOW_FREQUENCY_COLOR_R':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-LOW_FREQUENCY_COLOR_G':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-LOW_FREQUENCY_COLOR_B':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HIGH_FREQUENCY_COLOR_R':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HIGH_FREQUENCY_COLOR_G':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-HIGH_FREQUENCY_COLOR_B':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
	if sys.argv[t].upper() == '-NUMBER_PROCESSORS':
	    arguments[sys.argv[t].upper().replace('-','')] = int(sys.argv[t+1])
	    sys.argv[t] = ''
	    sys.argv[t+1] = ''
    
    # OUTPUT_BASENAME needs to end in a period
    if arguments['OUTPUT_BASENAME'] == '':
	arguments['OUTPUT_BASENAME'] = arguments['TRAJECTORY_FILENAME'] + "."
    if arguments['OUTPUT_BASENAME'][-1] != '.': arguments['OUTPUT_BASENAME'] = arguments['OUTPUT_BASENAME'] + '.'
    
    # check to see if anything is missing
    #print arguments['TRAJECTORY_FILENAME']
    if arguments['TRAJECTORY_FILENAME'] == '' and not os.path.exists(arguments['OUTPUT_BASENAME'] + 'average_hbonds'):
	program_output = program_output + '# ERROR: You have not specified a trajectory filename using the -trajectory_filename tag, and the file ' + arguments['OUTPUT_BASENAME'] + 'average_hbonds' + ',\n# which might contain a hydrogen-bond analysis of the trajectory, does not exist. Consequently, there are no hydrogen-bond networks\n# to analyze! Aborting!' + "\n"
	print program_output
	sys.exit(0)
    
    # SINGLE_FRAME_FILENAME, if not specified, will be set to the PDB file where frequencies are written
    extra_file_string = "occupancy"
    if arguments['SINGLE_FRAME_FILENAME'] == '':
	if arguments['WRITE_COLUMN'] == "O":
	    extra_file_string = "occupancy"
	elif arguments['WRITE_COLUMN'] == "B":
	    extra_file_string = "beta"
	arguments['SINGLE_FRAME_FILENAME'] = arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb"

    # report all command-line parameters not used
    notused = []
    for t in range(1,len(sys.argv)):
	if sys.argv[t] != "": notused.append(sys.argv[t])
    if len(notused) > 0:
	program_output = program_output + "# WARNING: The following command-line parameters were not used: "
	for item in notused: program_output = program_output + item + " "
	program_output = program_output + "\n\n"
	
    # now report values used
    program_output = program_output + "# The following parameters were used:\n"
    program_output = program_output + "#\ttrajectory_filename:\t\t" + arguments['TRAJECTORY_FILENAME'] + "\n"
    program_output = program_output + "#\toutput_basename:\t\t" + arguments['OUTPUT_BASENAME'] + "\n"
    program_output = program_output + "#\tsingle_frame_filename:\t\t" + arguments['SINGLE_FRAME_FILENAME'] + "\n"
    program_output = program_output + "#\tnumber_processors:\t\t" + str(arguments['NUMBER_PROCESSORS']) + "\n"
    program_output = program_output + "#\thydrogen_bond_frequency_cutoff:\t" + str(arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) + "\n"
    
    if len(arguments['SEED_RESIDUES']) != 0:
	seeds = ""
	for aseed in arguments['SEED_RESIDUES']:
	    seeds = seeds + str(aseed) + " "
	program_output = program_output + "#\tseed_residues:\t\t\t" + seeds + "\n"
    else:
	program_output = program_output + "#\tseed_residues:\t\t\tall (i.e., because no seed residues were specified, all hydrogen bonds will be visualized.)\n"
	
    program_output = program_output + "#\twrite_column:\t\t\t" + arguments['WRITE_COLUMN'] + "\n"
    program_output = program_output + "#\thydrogen_bond_distance_cutoff:\t" + str(arguments['HYDROGEN_BOND_DISTANCE_CUTOFF']) + "\n"
    program_output = program_output + "#\thydrogen_bond_angle_cutoff:\t" + str(arguments['HYDROGEN_BOND_ANGLE_CUTOFF']) + "\n"
    program_output = program_output + "#\tdashed_bonds:\t\t\t" + str(arguments['DASHED_BONDS']) + "\n"
    program_output = program_output + "#\tjust_immediate_connections:\t" + str(arguments['JUST_IMMEDIATE_CONNECTIONS']) + "\n"
    program_output = program_output + "#\tlow_frequency_color_r:\t\t" + str(arguments['LOW_FREQUENCY_COLOR_R']) + "\n"
    program_output = program_output + "#\tlow_frequency_color_g:\t\t" + str(arguments['LOW_FREQUENCY_COLOR_G']) + "\n"
    program_output = program_output + "#\tlow_frequency_color_b:\t\t" + str(arguments['LOW_FREQUENCY_COLOR_B']) + "\n"    
    program_output = program_output + "#\thigh_frequency_color_r:\t\t" + str(arguments['HIGH_FREQUENCY_COLOR_R']) + "\n"
    program_output = program_output + "#\thigh_frequency_color_g:\t\t" + str(arguments['HIGH_FREQUENCY_COLOR_G']) + "\n"
    program_output = program_output + "#\thigh_frequency_color_b:\t\t" + str(arguments['HIGH_FREQUENCY_COLOR_B']) + "\n"

    program_output = program_output + "\n"

    if arguments['SINGLE_FRAME_FILENAME'] == '' or arguments['SINGLE_FRAME_FILENAME'] == arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb": # no single_frame_filename specified
	if arguments['TRAJECTORY_FILENAME'] != "": # so a trajectory file does exist, we'll be getting the first frame from there.
	    program_output = program_output + "# The -single_frame_filename tag was not specified, so I'll use the first frame from the trajectory file, " + arguments['TRAJECTORY_FILENAME'] + ".\n\n"
	elif os.path.exists(arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb"): # no trajectory, but a frame had been saved previously using the same OUTPUT_BASENAME tag. Use that.
	    program_output = program_output + "# The -single_frame_filename tag was not specified, but the file " + arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb exists.\n# I'll use that filename as the -single_frame_filename tag."
	    arguments['SINGLE_FRAME_FILENAME'] == arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb"
    
    
    return arguments

class athread:
    #def __init__(self):
	#self.hbond_analyses = {}

    #def value_func(self, running, mutex, pdbs_lines, total_pdbs_in_traj, queue, counting_queue):

    # Takes a list of lists, where each list has the lines of a pdb file. It also accepts the total number of frames in the trajectory for updating purposes.
    # Returns a dictionary that describes the hydrogen bonds
    def value_func(self, running, mutex, pdb_filename, queue, counting_queue): # calculate the hydrogen bonds of the pdb file
	
	hbonds_per_frame = {}
	
	N_O_F_S = None
	hydrogen_indicies = None
	
	frame_index = -1
	
	current_pdb_lines = []
	file = open(pdb_filename,'r')
	while True: 
	    line = file.readline()
	    if len(line) == 0: break # EOF
	    if line[0:4]=="ATOM" or line[0:6]=="HETATM":
		current_pdb_lines.append(line)
	    if line[:13] == "REMARK FRAME ":
		frame_index = int(line[13:])
		#print frame_index
		
	    if (line[:3] == "END" or len(line) == 0) and len(current_pdb_lines) != 0: # "END" or end of file. Also, not an empty frame.
		val = counting_queue.get()  # an array containing all the frame_indices
		if not frame_index in val: val.append(frame_index)
		counting_queue.put(val)
		
		print "#\t\tProcessing frame " + str(len(val)) + "..." #+ " of " + str(count_used_frames) + "..."
	
		# create a PDB object, loading in from the list of lines
		if len(current_pdb_lines) > 0:
		
		    pdb = PDB() # make pdb object
		    pdb.LoadPDB_FromArray(current_pdb_lines)
		    
		    # if this is not the first frame, pass the lists that will be the same for each frame to the PDB object, so it doesn't need to be recalculated
		    if N_O_F_S is not None: 
			pdb.N_O_F_S = N_O_F_S
			pdb.hydrogen_indicies = hydrogen_indicies
	    
		    # now do the hydrogen bond analysis
		    hbonds_for_this_pdb = pdb.id_hydrogen_bonds(arguments) # id hydrogen bonds
		    
		    # if this is the first frame, define the atom lists that will be the same for all frames
		    if N_O_F_S is None:
			N_O_F_S = pdb.N_O_F_S
			hydrogen_indicies = pdb.hydrogen_indicies
	    
		    # now add that analysis into the larger, cumulative data structure (averages) and  the hbonds_per_frame
		    for hbond in hbonds_for_this_pdb: # add to the dictionary containing the list of hydrogen bonds
			key = str(hbond[0]) + "_" + str(hbond[1]) + "_" + str(hbond[2])
			
			# add into hbonds_per_frame
			if not key in hbonds_per_frame.keys(): # so the key doesn't exists
			    hbonds_per_frame[key] = [] # {}
			hbonds_per_frame[key].append(frame_index) # the first one should be the frame *******
    
		current_pdb_lines = []

	
	file.close()
    
	queue.put(hbonds_per_frame)
    
	#return hbonds_per_frame
	#queue.put(calculate_hydrogen_bonds_from_pdb_filename(pdb_filename, counting_queue))
	
        mutex.acquire()
        running.value -= 1
        mutex.release()

def control_multiple_processors(pdbs_divided_files):
    
    # count the total number of executions, across all chunks
    #total_pdbs_in_traj = 0
    #for pdb_lines in pdbs_lines:
	#total_pdbs_in_traj = total_pdbs_in_traj + len(pdb_lines)

    TASKS = len(pdbs_divided_files)
    running = Value('i', TASKS)
    mutex = Lock()

    #arrays = []
    threads = []
    queues = []
    for i in range(TASKS):
	threads.append(athread())
	queues.append(Queue())
	#arrays.append(multiprocessing.Array('i',[0, 1]))
	#arrays.append(multiprocessing.Array('_objects',[{}]))
	
    counting_queue = Queue()
    counting_queue.put([])

    processes = []
    for i in range(TASKS):
	#p = Process(target=threads[i].value_func, args=(running, mutex, pdbs_divided_files[i], total_pdbs_in_traj, queues[i], counting_queue))
	p = Process(target=threads[i].value_func, args=(running, mutex, pdbs_divided_files[i], queues[i], counting_queue))
	p.start()
	processes.append(p)

    while running.value > 0:
	mutex.acquire()
	sys.stdout.flush()
	mutex.release()
    
    hbond_analyses = []
    for i in range(TASKS):
	hbond_analyses.append(queues[i].get())
    
    return hbond_analyses, counting_queue
    
# Split a list into sublists of more or less same size
# Accepts as input a list and a number
# As output, returns a lists of lists.
def split_list(alist, num):
    lists = []
    for idx in range(num):
        lists.append(alist[idx::num])
    return lists

# a function that identifies the hydrogen bonds from a PDB trajectory and saves the frequency to a file
# returns a PDB object of the first frame, in case it's needed later
def identify_hydrogen_bonds_from_trajectory(pdb_filenames):
    
    global arguments, program_output
    
    # load pdb data, separating it into chunks (separated by END)
    program_output = program_output + "\n# Identify hydrogen bonds from the trajectory " + arguments['TRAJECTORY_FILENAME'] + "\n"
    
    #*****
    
    # now calculate the hydrogen bonds for each frame
    program_output = program_output + "#\tIdentifying hydrogen bonds...\n"
    print program_output.rstrip()
    program_output = ""
    
    if "WINDOWS" in platform.system().upper(): # so it's running on windows
        single_thread = athread()
        counting_queue = Queue()
        counting_queue.put([])
        aqueue = Queue()
        single_thread.value_func(Value('i', 1), Lock(), pdb_filenames[0], aqueue, counting_queue)
        hbond_analyses = [aqueue.get()]
    else:
        hbond_analyses,counting_queue = control_multiple_processors(pdb_filenames)
    
    used_frame_indices = counting_queue.get()
    
    #hbond_analyses = []
    #for chunk in all_pdbs_chunks:
	#hbond_analyses.append(calculate_hydrogen_bonds_from_list_of_pdb_lines(chunk, len(used_frame_indices)))

    # now merge all this hydrogen bond information into a single super dictionary
    hbonds_per_frame = {}
    for hbond_analysis in hbond_analyses:
	for hbond_key in hbond_analysis.keys():
	    
	    # if the hbond_key is not in the hbonds_per_frame key set, add it
	    if not hbond_key in hbonds_per_frame.keys(): hbonds_per_frame[hbond_key] = []
	    
	    # now go through the information associated with that key and add it into the larger hbonds_per_frame
	    for frame_index in hbond_analysis[hbond_key]: hbonds_per_frame[hbond_key].append(frame_index)

    # now write out the frame-by-frame hydrogen bonds to a csv file
    averages = {}
    lines = {}
    used_frame_indices.sort()
    
    #******

    for key in hbonds_per_frame: # for each hydrogen bond description
	
	theline = ""
	
	theline = theline + key + ","
	
	# write the hydrogen bonds
	hbonds_in_this_frame = hbonds_per_frame[key]
	for frame_index in used_frame_indices:
	    if frame_index in hbonds_in_this_frame:
		theline = theline + "1,"
		#val = hbonds_in_this_frame[frame_index]
		#theline = theline + str(val) + ","
	    else:
		theline = theline + "0,"
	
	anaverage = float(len(hbonds_in_this_frame)) / float(len(used_frame_indices)) # average(nums)
	theline = theline + "," + str(anaverage)
	
	lines[theline] = anaverage
	averages[key] = anaverage # keep track of averages for later output
    
    # now sort by the averages
    lines = sorted(lines.iteritems(), key=operator.itemgetter(1), reverse=True) # sort by percentage to make it easier to read

    # write to the file
    # first, write out the header
    program_output = program_output + "#\tWriting frame-by-frame hydrogen bond analysis to the file " + arguments['OUTPUT_BASENAME'] + "frame_by_frame_hbonds.csv\n"
    file = open(arguments['OUTPUT_BASENAME'] + "frame_by_frame_hbonds.csv",'w')
    file.write('H_HD_HA,')
    for frame_index in used_frame_indices: file.write("frame_" + str(frame_index) + ",")
    file.write(",Average")
    file.write("\n")
    # now write the meat
    for key in lines: file.write(key[0] + "\n")
    file.close()
    
    # now sort the average results
    averages = sorted(averages.iteritems(), key=operator.itemgetter(1), reverse=True) # sort by percentage to make it easier to read
    
    # now write the average results to a file
    program_output = program_output + "#\tWriting the average hydrogen bonds to the file " + arguments['OUTPUT_BASENAME'] + "average_hbonds\n"
    file = open(arguments['OUTPUT_BASENAME'] + "average_hbonds",'w')
    for hbond in averages:
	file.write(hbond[0] + "\t" + str(hbond[1]) + "\n")
    file.close()
    
    program_output = program_output + "\n"
    

# A nice function to identify the VMD color id based on a numerical value (hbond frequency).
# Returns the VMD color index
def get_color_id(val):
    global arguments
    if val == arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']: return 22 # blue
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 23
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 2 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 24
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 3 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 25
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 4 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 26
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 5 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 27
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 6 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 28
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 7 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 29
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 8 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 30
    elif val <= arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] + 9 * (1.0 - arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']) / 10: return 31
    elif val <= 1.0: return 32 # red

# Accepts as input a list of numbers
# Returns the average of the values
def average(list):
    sum = 0.0
    count = 1
    for item in list:
	sum = sum + item
	count = count + 1
    return sum / (float(count)-1)

# A function to set up the tcl visualization. Sets up things like colors, atom radii, etc.
def set_up_tcl_visualization():
    global arguments, program_output
    program_output = program_output + "# load in the pdb\n"
    program_output = program_output + 'mol new ' + arguments['SINGLE_FRAME_FILENAME'] + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all' + "\n"
    program_output = program_output + "\n"
    
    program_output = program_output + "# set the color scale\n"
    program_output = program_output + 'color scale method BGR' + "\n"
    program_output = program_output + "\n"
    program_output = program_output + '# redifine colors 22 - 32' + "\n"
    
    start_r_frac = arguments['LOW_FREQUENCY_COLOR_R'] / 255.0
    start_g_frac = arguments['LOW_FREQUENCY_COLOR_G'] / 255.0
    start_b_frac = arguments['LOW_FREQUENCY_COLOR_B'] / 255.0
    
    end_r_frac = arguments['HIGH_FREQUENCY_COLOR_R'] / 255.0
    end_g_frac = arguments['HIGH_FREQUENCY_COLOR_G'] / 255.0
    end_b_frac = arguments['HIGH_FREQUENCY_COLOR_B'] / 255.0
    
    delta_r_frac = (end_r_frac - start_r_frac) / 10.0
    delta_g_frac = (end_g_frac - start_g_frac) / 10.0
    delta_b_frac = (end_b_frac - start_b_frac) / 10.0
    
    index = 22
    for t in range(11):
	program_output = program_output + 'color change rgb ' + str(int(index)) + ' ' + str(start_r_frac + delta_r_frac * t) + ' ' + str(start_g_frac + delta_g_frac * t) + ' ' + str(start_b_frac + delta_b_frac * t) + "\n" # blue 
	index = index + 1
	
    program_output = program_output + "\n"
    program_output = program_output + "# variables to control appearance\n"
    program_output = program_output + 'set hbond_radius 0.1' + "\n"
    program_output = program_output + 'set atom_radius_scale 1.0' + "\n"
    
    program_output = program_output + "\n"
    program_output = program_output + "# calculate the radii of the atoms\n"
    program_output = program_output + '# first, the default radii' + "\n"
    #program_output = program_output + 'set H_radius 0.25' + "\n"
    program_output = program_output + 'set H_radius 0.5' + "\n"
    program_output = program_output + 'set N_radius 0.65' + "\n"
    program_output = program_output + 'set O_radius 0.60' + "\n"
    program_output = program_output + 'set F_radius 0.50' + "\n"
    program_output = program_output + "\n"
    program_output = program_output + '# now, calculate the drawing radii' + "\n"
    program_output = program_output + 'set actual_H_radius [expr $H_radius*$atom_radius_scale]' + "\n"
    program_output = program_output + 'set actual_N_radius [expr $N_radius*$atom_radius_scale]' + "\n"
    program_output = program_output + 'set actual_O_radius [expr $O_radius*$atom_radius_scale]' + "\n"
    program_output = program_output + 'set actual_F_radius [expr $F_radius*$atom_radius_scale]' + "\n"

# Accepts as input a list of the lines from the hbond file created by identify_hydrogen_bonds_from_trajectory
# Returns a dictionary containing the residues that are connected to each other through hydrogen bonds
def determine_residue_connections(hbonds_lines):
    global arguments
    
    # determine which residues are connected to each other through hydrogen bonds
    residue_connections = {}
    for line in hbonds_lines:
	line = line.strip()
	line = line.split("\t")
	percent = float(line[1])
	
	if arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] <= percent: # so ignore all bonds with frequency less than cutoff
	    key = line[0]
	    indicies = key.split("_")
	    #index1 = int(indicies[0]) # this is the hydrogen atom
	    index2 = int(indicies[1])
	    index3 = int(indicies[2])
	    #atom1 = single_pdb.AllAtoms[index1] # this is the hydrogen atom
	    atom2 = single_pdb.AllAtoms[index2]
	    atom3 = single_pdb.AllAtoms[index3]
	    
	    # update the residue_connections dictionary
	    if atom2.resid in residue_connections.keys():
		if not atom3.resid in residue_connections[atom2.resid]: residue_connections[atom2.resid].append(atom3.resid)
	    else:
		residue_connections[atom2.resid] = [atom3.resid]
	    
	    if atom3.resid in residue_connections.keys():
		if not atom2.resid in residue_connections[atom3.resid]: residue_connections[atom3.resid].append(atom2.resid)
	    else:
		residue_connections[atom3.resid] = [atom2.resid]
    return residue_connections

# Accepts as input a dictionary specifying which residues are connected to each other through hydrogen bonds
# Returns a list containing all the residues that are connected to the user-specified seed residues
def residues_connected_to_the_seeds(residue_connections):
    global arguments, program_output
    
    # first we need to see if there are any residues
    connections_to_seed = False
    if len(arguments['SEED_RESIDUES']) == 0: # so none specified, doing all.
	connections_to_seed = True
    else:
	for resid in arguments['SEED_RESIDUES']:
	    if resid in residue_connections.keys() > 0:
		connections_to_seed = True
		break
    if connections_to_seed == False:
	program_output = program_output + "\n# ERROR: Given your hydrogen-bond criteria (hydrogen_bond_frequency_cutoff, hydrogen_bond_distance_cutoff, and hydrogen_bond_angle_cutoff),\n# there are no residues connected through hydrogen bonds to the seed residues. Aborting!"
	print program_output
	sys.exit()
    
    # now we need to determine which residues are connected by hydrogen bonds to the seed residues
    connected_residues_growing_list = []
    if len(arguments['SEED_RESIDUES']) > 0: # so a seed as been specified
	connected_residues_growing_list = arguments['SEED_RESIDUES'][:]
	if arguments['JUST_IMMEDIATE_CONNECTIONS'] == "FALSE": # so you're going to "recursively" walk through all neighbors
	    index = 0
	    while index < len(connected_residues_growing_list):
		theresid = connected_residues_growing_list[index]
		for connecteds in residue_connections[theresid]:
		    if not connecteds in connected_residues_growing_list: connected_residues_growing_list.append(connecteds)
		index = index + 1
	else: # so you only care about the immediate neighbors to the seed residues
	    tmp = []
	    for theresid in connected_residues_growing_list:
		for connecteds in residue_connections[theresid]:
		    if not connecteds in tmp: tmp.append(connecteds)
	    for t in tmp:
		if t not in connected_residues_growing_list:
		    connected_residues_growing_list.append(t)
    else: # so no seed specified
	for index in single_pdb.AllAtoms:
	    atom = single_pdb.AllAtoms[index]
	    if atom.resid not in connected_residues_growing_list: connected_residues_growing_list.append(atom.resid)
	    
    return connected_residues_growing_list

# A function to draw cylinders. x1, y1, z1 is the coordinate of the start of the cylinder, and x2, y2, z2 is the cylinder's end. dashed is a boolean variable to draw a dashed cylinder (vs solid).
def draw_cylinder(x1, y1, z1, x2, y2, z2, dashed="TRUE"):
	if dashed=="FALSE":
		return 'graphics top cylinder {' + str(x1) + ' ' + str(y1) + ' ' + str(z1) + '} {' + str(x2) + ' ' + str(y2) + ' ' + str(z2) + '} radius $hbond_radius resolution 10 filled yes' + "\n"
	else:
		length_dash = 0.5

		# figure out which x1 is lower
		if x1<x2:
			lower_x = x1
			lower_y = y1
			lower_z = z1
			higher_x = x2
			higher_y = y2
			higher_z = z2
		else:
			lower_x = x2
			lower_y = y2
			lower_z = z2
			higher_x = x1
			higher_y = y1
			higher_z = z1

		deltax = higher_x-lower_x
		deltay = higher_y-lower_y
		deltaz = higher_z-lower_z
		dist = math.sqrt(math.pow(deltax,2) + math.pow(deltay,2) + math.pow(deltaz,2))
		deltax = (deltax / dist) * length_dash
		deltay = (deltay / dist) * length_dash
		deltaz = (deltaz / dist) * length_dash
		current_x = lower_x
		current_y = lower_y
		current_z = lower_z
		
		toprintout = ""
		while current_x < higher_x:
			toprintout=toprintout+'graphics top cylinder {' + str(current_x) + ' ' + str(current_y) + ' ' + str(current_z) + '} {' + str(current_x + deltax/2) + ' ' + str(current_y + deltay/2) + ' ' + str(current_z + deltaz/2) + '} radius $hbond_radius resolution 10 filled yes' + "\n"
			current_x = current_x + deltax
			current_y = current_y + deltay
			current_z = current_z + deltaz
		return toprintout
		
# Accepts as input a pdb object containing one frame and a list containing all the residues that are connected through hydrogen bonds to the user-specified seed residues
# Aside from outputing the tcl commands to draw the hydrogen bonds, returns a dictionary mapping the hydrogen-bond atoms to a list of all the frequencies for their associated hydrogen bonds.
# This dictionary will ultimately be used to calculate an average frequency for the relevant atoms, so they can be colored appropriately even if they participate in multiple hydrogen bonds.
def tcl_draw_hydrogen_bonds(single_pdb, connected_residues_growing_list):
    global arguments, program_output
    
    program_output = program_output + "\n"
    program_output = program_output + '# draw hydrogen bonds' + "\n"
    
    atom_vals = {}
    
    for line in hbonds_lines:
	line = line.strip()
	line = line.split("\t")
	percent = float(line[1])
	
	if arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF'] <= percent: # so ignore all bonds with frequency less than cutoff
	    
	    key = line[0]
	    indicies = key.split("_")
	    index1 = int(indicies[0])
	    index2 = int(indicies[1])
	    index3 = int(indicies[2])
	    
	    #print line, index1, index2, index3, percent
	    #atom1 = single_pdb.AllAtoms[index1] # this is the hydrogen atom
	    atom2 = single_pdb.AllAtoms[index2]
	    atom3 = single_pdb.AllAtoms[index3]
		    
	    if atom2.resid in connected_residues_growing_list and atom3.resid in connected_residues_growing_list: # so only if somehow connected to the original seed residues
    
		# tally up values for each atom (to be averaged later)
		if index1 in atom_vals.keys():
		    atom_vals[index1].append(percent)
		else:
		    atom_vals[index1] = [percent]
		    
		if index2 in atom_vals.keys():
		    atom_vals[index2].append(percent)
		else:
		    atom_vals[index2] = [percent]
		    
		if index3 in atom_vals.keys():
		    atom_vals[index3].append(percent)
		else:
		    atom_vals[index3] = [percent]
	    
		program_output = program_output + 'graphics top color ' + str(get_color_id(percent)) + "\n"
		program_output = program_output + draw_cylinder(atom2.coordinates.x, atom2.coordinates.y, atom2.coordinates.z, atom3.coordinates.x, atom3.coordinates.y, atom3.coordinates.z, arguments['DASHED_BONDS'])
		#program_output = program_output + 'graphics top cylinder {' + str(atom2.coordinates.x) + ' ' + str(atom2.coordinates.y) + ' ' + str(atom2.coordinates.z) + '} {' + str(atom3.coordinates.x) + ' ' + str(atom3.coordinates.y) + ' ' + str(atom3.coordinates.z) + '} radius $hbond_radius resolution 10 filled no' + "\n"
    
    return atom_vals

# Accepts as input a pdb object of a single frame, and a dictionary specifying the frequencies of all the hydrogen bonds each relevant atom participates in.
# Aside form updating the occupancy or beta column of the pdb object with the average frequency of the relevant atoms, returns a dictionary containing the average values for the relevant atoms.
def calculate_average_frequency_relevant_atoms_participate_in_hydrogen_bonds(single_pdb, atom_vals):
    
    global arguments
    
    atom_average_vals = {}
    for key in atom_vals.keys():
	atom_average_vals[key] = average(atom_vals[key])
	if arguments['WRITE_COLUMN'] == "O":
	    single_pdb.AllAtoms[key].occupancy = atom_average_vals[key]
	    extra_file_string = "occupancy"
	elif arguments['WRITE_COLUMN'] == "B":
	    single_pdb.AllAtoms[key].beta = atom_average_vals[key]
	    extra_file_string = "beta"
    
    if extra_file_string != "": single_pdb.SavePDB(arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb")
    
    return atom_average_vals
    
# Accepts as input a pdb object containing a single frame, the average frequency values for the atoms, and a list of the residues that are connected via hydrogen bonds to the seed residues.
# Aside from outputing the tcl commands to draw the atoms, also outputs a list of all the residues connected to the seed residue via hydrogen bonds
def draw_atoms(single_pdb, atom_average_vals, connected_residues_growing_list):

    global program_output

    program_output = program_output + "\n"
    program_output = program_output + '# draw atoms that participate in hydrogen bonds' + "\n"
    
    resids = []
    
    for index in atom_average_vals.keys():
	percent = atom_average_vals[index]
	atom = single_pdb.AllAtoms[index]
	if atom.resid in connected_residues_growing_list: # so only if somehow connected to the original seed residues
	    if atom.resid not in resids: resids.append(atom.resid)
	    program_output = program_output + 'graphics top color ' + str(get_color_id(percent)) + "\n"
	    program_output = program_output + 'graphics top sphere {' + str(atom.coordinates.x) + ' ' + str(atom.coordinates.y) + ' ' + str(atom.coordinates.z) + '} radius $actual_' + atom.element + '_radius' + "\n"
    
    return resids

# Accepts as input a list of all the residues connected to the seed residue via hydrogen bonds.
# Writes tcl commands to improve the visualization in VMD
def final_tcl_representations(resids):
    
    global program_output, resid_text, arguments
    
    # now add in some representations
    resid_text = ""
    for resid in resids: resid_text = resid_text + str(resid) + ' '
    program_output = program_output.replace('{CONNECTED_RESID}',resid_text.strip())
    
    seeds_text = ""
    for resid in arguments['SEED_RESIDUES']: seeds_text = seeds_text + str(resid) + ' '
    
    program_output = program_output + "\n"
    program_output = program_output + "# add some representations\n"
    program_output = program_output + 'mol delrep 0 top' + "\n"
    
    program_output = program_output + "\n"
    program_output = program_output + "# make all residues that have hydrogen bonds appear in licorice" + "\n"
    
    program_output = program_output + 'mol representation Licorice 0.100000 10.000000 10.000000' + "\n"
    program_output = program_output + 'mol color Name' + "\n"
    program_output = program_output + 'mol selection {noh and resid ' + resid_text.strip() + '}' + "\n"
    program_output = program_output + 'mol material Opaque' + "\n"
    program_output = program_output + 'mol addrep top' + "\n"
    program_output = program_output + 'mol selupdate 0 top 0' + "\n"
    program_output = program_output + 'mol colupdate 0 top 0' + "\n"
    
    program_output = program_output + "\n"
    program_output = program_output + "# draw the protein in NewCartoon" + "\n"
    
    program_output = program_output + 'mol representation NewCartoon 0.300000 6.000000 4.100000 0' + "\n"
    program_output = program_output + 'mol color Structure' + "\n"
    program_output = program_output + 'mol selection {all}' + "\n"
    program_output = program_output + 'mol material Opaque' + "\n"
    program_output = program_output + 'mol addrep top' + "\n"
    program_output = program_output + 'mol selupdate 0 top 0' + "\n"
    program_output = program_output + 'mol colupdate 0 top 0' + "\n"
    
    
    program_output = program_output + "\n"
    program_output = program_output + "# emphasize the seed residues" + "\n"
    
    program_output = program_output + 'mol representation Licorice 0.500000 10.000000 10.000000' + "\n"
    program_output = program_output + 'mol color Name' + "\n"
    #program_output = program_output + 'mol selection {not protein and not nucleic and not water and not mass 23 39}' + "\n" # note that this doesn't get rid of Cl. Why? Because it could be in the
    program_output = program_output + 'mol selection {resid ' + seeds_text.strip() + '}' + "\n" 
    
    program_output = program_output + 'mol material Opaque' + "\n"
    program_output = program_output + 'mol addrep top' + "\n"
    program_output = program_output + 'mol selupdate 0 top 0' + "\n"
    program_output = program_output + 'mol colupdate 0 top 0' + "\n"
    
    program_output = program_output + "\n"
    program_output = program_output + '# additional modifications to improve visualization' + "\n"
    program_output = program_output + 'display projection Orthographic' + "\n"
    program_output = program_output + 'display depthcue on' + "\n"


def divide_trajectory_file():
    
    global arguments, program_output
    
    # first, get the total number of frames, look for any errors
    print "\n# Counting the number of frames in the trajectory..."
    file = open(arguments['TRAJECTORY_FILENAME'],'r')
    count_frames = 0
    num_atoms = -1
    count_atoms = 0
    lines_of_first_pdb = []
    printnumreduce = 0

    while True:
	line = file.readline()
	if (line[:5] == "ATOM " or line[:7] == "HETATM "):
	    count_atoms = count_atoms + 1
	    if num_atoms == -1: lines_of_first_pdb.append(line.strip())
	if line[:3] == "END" or len(line) == 0: # so END tag or EOF
	    count_frames = count_frames + 1
	    if num_atoms == -1: # get information about the number of atoms per frame
		num_atoms = count_atoms
		print "#\tNumber of atoms per frame: " + str(num_atoms)
	    elif count_atoms != 0: # no different atoms per frame errors if empty, in which case it would be skipped (see below)
		if count_atoms != num_atoms: # so not the same number of atoms as the first frame
		    print "#\tERROR: Frame " + str(count_frames) + " contains " + str(count_atoms) + " atoms. This is a different number of atoms than the first frame, which had " + str(num_atoms) + " atoms. Aborting program."
		    sys.exit()
	    else: # so an empty frame
		if len(line) <> 0: 
			print "#\tWARNING: Frame " + str(count_frames) + " appears to contain no atoms. Skipping..."
		else:
			printnumreduce = -1 # the last frame

	    count_atoms = 0
	if len(line) == 0: break # EOF
	
    file.close()
    count_frames = count_frames + 1 # in case the last frame didn't end in END, just add one

    # if the number of frames is less than the nnumber of processors, reset.
    if count_frames < arguments['NUMBER_PROCESSORS']:
	print "#\tWARNING: Too few frames to use " + str(arguments['NUMBER_PROCESSORS']) + " procesors. Using " + str(count_frames-1) + " instead."
	arguments['NUMBER_PROCESSORS'] = count_frames - 1
    
    # on windows, multiprocessing must be turned off
    if "WINDOWS" in platform.system().upper():
	print "#\tWARNING: You seem to be using Microsoft Windows. Multiple processors not supported."
	arguments['NUMBER_PROCESSORS'] = 1

    print "#\tThere are " + str(count_frames-1+printnumreduce) + " frames." # I think this gives a correct count in all cases, but not sure. It's definitely close, though.
    #if count_frames-1+printnumreduce == 1: # so just one frame # this is not necessary
        #print "#\tWARNING: Because there is only one frame in the specified trajectory, the -hydrogen_bond_frequency_cutoff tag has been set to 1.0, regardless of the user-specified value."
        #arguments['HYDROGEN_BOND_FREQUENCY_CUTOFF']  = 1.0

    # now save the first frame in case it's needed.
    single_frame = PDB()
    single_frame.LoadPDB_FromArray(lines_of_first_pdb)
        
    # now divide the trajectory into smaller files
    frame_per_file = int(count_frames/arguments['NUMBER_PROCESSORS'])
    count_frames = 0
    total_frame_count = 1
    file_index = 1
    
    divided_files = []
    file_input = open(arguments['TRAJECTORY_FILENAME'],'r')
    file_output = open(arguments['OUTPUT_BASENAME'] + str(file_index) + ".tmp",'w')
    divided_files.append(arguments['OUTPUT_BASENAME'] + str(file_index) + ".tmp")
    next_frame_remark = "REMARK FRAME " + str(total_frame_count) + "\n"
    while True:
	line = file_input.readline()
	if len(line) == 0: break # EOF
	if next_frame_remark != "":
	    file_output.write(next_frame_remark)
	    next_frame_remark = ""
	if line[:5] == "ATOM " or line[:7] == "HETATM " or line[:3] == "END": file_output.write(line)
	if line[:3] == "END":
	    count_frames = count_frames + 1
	    total_frame_count = total_frame_count + 1
	    next_frame_remark = "REMARK FRAME " + str(total_frame_count) + "\n"
	if count_frames == frame_per_file and file_index < arguments['NUMBER_PROCESSORS']: # so time to switch files
	    count_frames = 0
	    file_index = file_index + 1
	    file_output.close()
	    file_output = open(arguments['OUTPUT_BASENAME'] + str(file_index) + ".tmp",'w')
	    divided_files.append(arguments['OUTPUT_BASENAME'] + str(file_index) + ".tmp")
    
    if not "END" in line: file_output.write("END\n") # in case the PDB doesn't end in "END"
    
    file_output.close()
    file_input.close()

    return (single_frame,divided_files)

# get a variable to store protein text output in
program_output = ""

# initial program info and some additional useful info
program_info()
program_output = program_output + "# This TCL output file was designed to be loaded into the Visual Molecular Dynamics (VMD) computer program using the -e tag from the\n# command line, or through the in-program Tk Console. See below for instructions describing how to use this program's output in other\n# visualization programs.\n\n"
program_output = program_output + "# Run this python script using the \"-help\" tag to get a detailed help file, including examples demonstrating usage.\n\n"

# get the commandline parameters
extra_file_string = ""
arguments = get_commandline_parameters()

# print out and reset program output
print program_output.rstrip()
program_output = ""

# identify the hydrogen bonds from the trajectory if necessary
single_frame_from_trajectory = None
if os.path.exists(arguments['OUTPUT_BASENAME'] + "average_hbonds"):
    program_output = program_output + "\n# The file " + arguments['OUTPUT_BASENAME'] + "average_hbonds, which may contain information about the hydrogen bonds of the selected trajectory, already exists.\n# I'm not going to recalculate the hydrogen bonds from the trajectory file, so the -hydrogen_bond_distance_cutoff and\n# -hydrogen_bond_angle_cutoff tags you specified will be ignored. If you'd like to recalculate, perhaps with new\n# -hydrogen_bond_distance_cutoff and -hydrogen_bond_angle_cutoff values, move or delete this file, or select a different value for\n# the -output_basename parameter.\n\n"
else: # so the file doesn't exist, you need to create it.
    single_frame_from_trajectory,divided_files = divide_trajectory_file() # divide the trajectory PDB file into smaller files of equal size.
    identify_hydrogen_bonds_from_trajectory(divided_files)
    for file in divided_files: os.remove(file) # delete temporary files

# description of how to use output with another visualization program
program_output = program_output + "# Visualizing the program's output in VMD is easy using the -e option from the command line or the in-program Tk Console. If you wish\n# to visualize the program's output in a different molecular visualization program, the following information may be helpful:\n"
program_output = program_output + "#\t1) The average frequency with which relevant atoms participate in hydrogen bonds has been stored in the {EXTRA_FILE_STRING}\n#\t   column of the PDB file " + arguments['OUTPUT_BASENAME'] + "hbond_averages_in_{EXTRA_FILE_STRING}_column.pdb\n"
program_output = program_output + "#\t2) Residues with the following IDs are connected to the user-specified seed residues via hydrogen bonds:\n#\t   {CONNECTED_RESID}\n\n"

# write tcl commands to set up the visualization
program_output = program_output + "# BEGIN VMD TCL FILE COMMANDS\n"
program_output = program_output + "# ===========================\n\n"
set_up_tcl_visualization()

# load single-frame pdb file
single_pdb = PDB()

if arguments['SINGLE_FRAME_FILENAME'] != "" and arguments['SINGLE_FRAME_FILENAME'] != arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb": # so a single-frame filename was specified, use that
    single_pdb.LoadPDB_FromFile(arguments['SINGLE_FRAME_FILENAME']) # must have same atom order and indicies.
elif single_frame_from_trajectory is not None: # so a trajectory has already been loaded, and the first frame was recorded
    single_pdb = single_frame_from_trajectory
elif arguments['TRAJECTORY_FILENAME'] != '': # so a trajectory has been specified, but not loaded, probably because the hydrogen bonds had already been calcuated from a previous run.
    # we need to load in just the first frame.
    file = open(arguments['TRAJECTORY_FILENAME'],'r')
    pdb_lines = []
    line = ""
    while line[:3] != "END":
	line = file.readline()
	if len(line) == 0: break # EOF
	line = line.strip()
	if (line[:5] == "ATOM " or line[:7] == "HETATM "):
	    pdb_lines.append(line)
	if line[:3] == "END": break # so it's at the end of the pdb frame
    file.close()
    single_pdb.LoadPDB_FromArray(pdb_lines)
elif os.path.exists(arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb"): # so use this file if it's been previously saved
    single_pdb.LoadPDB_FromFile(arguments['OUTPUT_BASENAME'] + "hbond_averages_in_" + extra_file_string + "_column.pdb") # must have same atom order and indicies.

# now let's make sure the seed residues actually exist in the PDB file
for seed_index in range(len(arguments['SEED_RESIDUES'])-1,-1,-1):
    seed = arguments['SEED_RESIDUES'][seed_index]
    if not seed in single_pdb.resids:
	print "# WARNING: One of the seed residues, residue " + str(seed) + ", does not exist in the file specified by the \"-single_frame_filename\" tag. I'll ignore that seed."
	arguments['SEED_RESIDUES'].pop(seed_index)

# load the hbonds file
file = open(arguments['OUTPUT_BASENAME'] + "average_hbonds", 'r')
hbonds_lines = file.readlines()
file.close()

# determine which residues are connected to each other through hydrogen bonds
residue_connections = determine_residue_connections(hbonds_lines)

# now we need to determine which residues are connected by hydrogen bonds to the seed residues
connected_residues_growing_list = residues_connected_to_the_seeds(residue_connections)

# write the tcl code to draw hydrogen bonds
atom_vals = tcl_draw_hydrogen_bonds(single_pdb, connected_residues_growing_list)

# calculate the average frequency relevant atoms participate in hydrogen bonds
atom_average_vals = calculate_average_frequency_relevant_atoms_participate_in_hydrogen_bonds(single_pdb, atom_vals)
program_output = program_output.replace('{EXTRA_FILE_STRING}',extra_file_string)

# draw the atoms that participate in the hydrogen bonds
resids = draw_atoms(single_pdb, atom_average_vals, connected_residues_growing_list)

# draw final tcl representations
final_tcl_representations(resids)

print program_output.rstrip()
