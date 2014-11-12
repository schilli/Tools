import sys, os
import numpy as np


class Molecule(object):

    def __init__(self):
        self.chainid     = []
        self.atomid      = []
        self.atomname    = []
        self.resname     = []
        self.resid       = []
        self.sequence    = []
        self.coordinates = []
        self.occupancy   = []
        self.Bfactor     = []
        self.element     = []
        self.charge      = []

    def postprocess(self):
        # create sequence
        self.sequence = []
        for i, ID in enumerate(self.resid):
            if i == 0 or ID != self.resid[i-1]:
                self.sequence.append(self.resname[i])

        # make numpy arrays
        self.coordinates = np.array(self.coordinates)
        self.atomid      = np.array(self.atomid, dtype=np.int)
        self.resid       = np.array(self.resid,  dtype=np.int)
        self.occupancy   = np.array(self.occupancy)
        self.Bfactor     = np.array(self.Bfactor)
        self.charge      = np.array(self.charge)


    def add_atom(self, PDBline):
        if PDBline[0:6].find("ATOM") or PDBline[0:6].find("HETATM"):
            try:
                self.atomid     .append( int( PDBline[ 6:11])       )
                self.atomname   .append(      PDBline[12:16].strip())
                self.resname    .append(      PDBline[17:20].strip())
                self.chainid    .append(      PDBline[21:22].strip())
                self.resid      .append( int( PDBline[22:26])       )
                self.coordinates.append(float(PDBline[30:38])       )
                self.coordinates.append(float(PDBline[38:46])       )
                self.coordinates.append(float(PDBline[46:54])       )

                try:
                    self.occupancy  .append(float(PDBline[54:60])   )
                except ValueError:
                    self.occupancy  .append(1.0)

                try:
                    self.Bfactor    .append(float(PDBline[60:66])   )
                except ValueError:
                    self.Bfactor    .append(1.0)

                self.element    .append(      PDBline[76:78].strip())
                if self.element[-1] == "":
                    # find first non digit character of atomname
                    self.element[-1] = self.atomname[-1].strip(string.digits)[0]

                try:
                    self.charge     .append(      PDBline[79:80].strip())
                except ValueError:
                    self.charge     .append('')

            except Exception as e:
                print(e.message)
                raise PDBreadError("Something wrong with atom on line: \"{}\"".format(PDBline.rstrip()))
        else:
            raise PDBreadError("Expected ATOM line, found: \"{}\"".format(PDBline.rstrip()))


    def select_index(self, indices):
        """Return molecule of only atoms with given indices"""
        molecule = Molecule()
        for i in range(len(self.atomid)):
            if i in indices:
                molecule.chainid  .append(self.chainid  [i])
                molecule.atomid   .append(self.atomid   [i])
                molecule.atomname .append(self.atomname [i])
                molecule.resname  .append(self.resname  [i])
                molecule.resid    .append(self.resid    [i])
                molecule.occupancy.append(self.occupancy[i])
                molecule.Bfactor  .append(self.Bfactor  [i])
                molecule.element  .append(self.element  [i])
                molecule.charge   .append(self.charge   [i])

                molecule.coordinates.append(self.coordinates[i*3+0])
                molecule.coordinates.append(self.coordinates[i*3+1])
                molecule.coordinates.append(self.coordinates[i*3+2]) 

        molecule.postprocess()
        return molecule



    def select_CA(self):
        """Return a molecule of only CA atoms"""
        indices = []
        for i in range(len(self.atomid)):
            if self.atomname[i].strip() == "CA":
                indices.append(i)
        molecule = self.select_index(indices)
        return molecule


    def select_noh(self):
        """Return a molecule of only non hydrogen atoms"""
        indices = []
        for i in range(len(self.atomid)):
            if self.element[i] != 'H':
                indices.append(i)
        molecule = self.select_index(indices)
        return molecule 

    def select_protein(self):
        """Return a molecule of only protein atoms."""
        indices = []
        aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                      "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                      "THR", "TRP", "TYR", "VAL", "ACE", "NME"]
        for i in range(len(self.atomid)):
            if self.resname[i] in aminoacids:
                indices.append(i)
        molecule = self.select_index(indices)
        return molecule  


    def select_elements(self, *elements):
        """Return a molecule of only atoms of the given elements
        Pass elements as list, or separate arguments"""
        selectedElements = []
        for item in elements:
            if type(item) == list or type(item) == tuple:
                for subitem in item:
                    if type(subitem) == str:
                        selectedElements.append(subitem)
                    else:
                        raise TypeError("Selection elements should be strings, lists or touples of strings")
            elif type(item) == str:
                selectedElements.append(item)
            else:
                raise TypeError("Selection elements should be strings, lists or touples of strings")

        indices = []
        for i in range(len(self.atomid)):
            if self.element[i] in selectedElements:
                indices.append(i)
        molecule = self.select_index(indices)
        return molecule


    def select_chains(self, *chains):
        """Return a molecule of only atoms of the given chains
        Pass chains as list, or separate arguments"""
        selectedChains = []
        for item in chains:
            if type(item) == list or type(item) == tuple:
                for subitem in item:
                    if type(subitem) == str:
                        selectedChains.append(subitem)
                    else:
                        raise TypeError("Selection chains should be strings, lists or touples of strings")
            elif type(item) == str:
                selectedChains.append(item)
            else:
                raise TypeError("Selection chains should be strings, lists or touples of strings")

        indices = []
        for i in range(len(self.atomid)):
            if self.chainid[i] in selectedChains:
                indices.append(i)
        molecule = self.select_index(indices)
        return molecule 


def read_pdb(filename):

    structures = []
    fh = open(filename, 'r')
    model    = 0
    readingStructure = False
    line     = fh.readline()
    lineNumber = 1

    while len(line) > 0:
        # if at the start of a new model
        if not readingStructure and (line.find("ATOM") == 0 or line.find("HETATM") == 0):
            molecule = Molecule()
            readingStructure = True

        # read atom
        if readingStructure:
            if line.find("ATOM") == 0 or line.find("HETATM") == 0:
                molecule.add_atom(line)
            elif line.find("TER") == 0:
                pass
            # finish up model
            else:
                molecule.postprocess()
                structures.append(molecule)
                readingStructure = False

        line = fh.readline()
        lineNumber += 1

        # finish up model
        if len(line) == 0 and readingStructure:
            molecule.postprocess()
            structures.append(molecule)
            readingStructure = False
 
    fh.close()

    return structures
    


def write_pdb(structures, filename):

    with open(filename, 'w') as fh:

        template = "{:<6s}{:5d} {:<4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}{:2s}\n"

        try:
            l = len(structures)
        except TypeError:
            structures = [structures]

        for s, structure in enumerate(structures):
            if len(structures) > 1:
                fh.write("MODEL        {}\n".format(s+1))
            for a in range(len(structure.atomid)):
                atomidOffset = 0
                if a > 0 and structure.chainid[a] != structure.chainid[a-1]:
                    fh.write("TER   {:5d} {:<4s} {:3s} {:1s}{:4d}\n".format(structure.atomid[a-1]+1,
                                                                            "",
                                                                            structure.resname[a-1],
                                                                            structure.chainid[a-1],
                                                                            structure.resid[a-1]))
                    if len(structure.atomid) > a-1 and structure.atomid[a] == structure.atomid[a-1]+1:
                        atomidOffset += 1
                line = template.format("ATOM", structure.atomid[a] + atomidOffset,
                                               structure.atomname[a],
                                               structure.resname[a],
                                               structure.chainid[a],
                                               structure.resid[a],
                                               structure.coordinates[3*a+0],
                                               structure.coordinates[3*a+1],
                                               structure.coordinates[3*a+2],
                                               structure.occupancy[a],
                                               structure.Bfactor[a],
                                               structure.element[a],
                                               structure.charge[a])
                fh.write(line)
            if len(structures) > 1:
                fh.write("ENDMDL\n")



class PDBreadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
