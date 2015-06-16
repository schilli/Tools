#!/usr/bin/env python
# Resort and renumber PDB files based on a template
#
# usage: pdbrenum.py template.pdb structure.pdb

from __future__ import print_function
import sys, time, copy
import numpy as np

# ============================================================================ #


if __name__ == "__main__":
    from pdbrenum import *
    pdb1 = readPDB(sys.argv[1])
    pdb2 = readPDB(sys.argv[2])

    pdb2sorted = []
    for model in pdb2:
        pdb2sorted.append(resortAtoms_name(pdb1[0], model))

    
# ============================================================================ #


class AtomEntry(object):
    """Class to hold PDB Atom information"""

# ==================================== #

    def __init__(self, line=None):
        self.line    = line        # original PDB entry text
        self.atom    = "ATOM"      # or "HETATM"
        self.num     = 0           # atom serial number
        self.name    = ""          # atom name
        self.alt     = ""          # alternate location indicator
        self.resname = ""          # residue name
        self.chain   = ""          # chain ID
        self.resid   = 0           # residue ID
        self.ins     = ""          # code for insertion of residues
        self.coords  = np.zeros(3) # coordinates
        self.occ     = 1.0         # occupancy
        self.Bfac    = 0.0         # B-Factor
        self.element = ""          # element symbol
        self.charge  = 0.0         # charge on element

        if line != None:
            self.read(line)

# ==================================== #

    def read(self, line):
        """Read in atom data from PDB line.
        Format taken from: http://cupnet.net/pdb-format/"""

        self.line    =       line
        self.atom    =       line[0:6].strip()
        self.num     =   int(line[6:11])
        self.name    =       line[12:16].strip()
        self.alt     =       line[16:17].strip()
        self.resname =       line[17:20].strip()
        self.chain   =       line[21:22].strip()
        try:
            self.resid =   int(line[22:26])
        except:
            self.resid = 1
        self.ins     =       line[26:27].strip()
        self.coords  = np.array([float(i) for i in [line[30:38], line[38:46], line[46:54]]])
        self.occ     = float(line[54:60])
        self.Bfac    = float(line[60:66])
        self.element =       line[76:78].strip()
        self.charge  =       line[78:80].strip()

# ==================================== #

    def write(self):
        """Write atom data in PDB format to string"""

        PDBformat  = "{:<6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "
        PDBformat += "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}{:2s}"

        atomLine = PDBformat.format(self.atom   , \
                                    self.num    , \
                                    self.name   , \
                                    self.alt    , \
                                    self.resname, \
                                    self.chain  , \
                                    self.resid  , \
                                    self.ins    , \
                                    self.coords[0], self.coords[1], self.coords[2], \
                                    self.occ    , \
                                    self.Bfac   , \
                                    self.element, \
                                    self.charge)
        return atomLine


# ============================================================================ #


def readAtoms(rawLines):
    """Read atom entries from a list of raw PDB atom line entries"""
    atomList = []
    for line in rawLines:
        if max(line[0:6].find("ATOM"), line[0:6].find("HETATM")) >= 0:
            atomList.append(AtomEntry(line)) 

    return atomList


# ============================================================================ #


def readPDB(filename):
    """Read PDB file and return a list of lists with all atoms per model"""

    with open(filename, 'r') as f:
        data = f.readlines()

    models = []

    # check if there is a model entry
    multipleModels = False
    for line in data:
        if line[:5].find("MODEL") >= 0:
            multipleModels = True

    if not multipleModels:
        # read all in one shot
        print("Hi there!")
        models.append(readAtoms(data))

    else:
        # read separate models
        rawLines = []
        for line in data:

            rawLines.append(line)

            if line[:5].find("MODEL") >= 0:
                rawLines = []
            elif line[:6].find("ENDMDL") >= 0:
                models.append(readAtoms(rawLines))


    return models


# ============================================================================ #


def writePDB(models, filename=None):
    """Write all models to PDB file"""
    if filename == None:
        outfile = sys.stdout
    else:
        outfile = open(filename, 'w')
    with open(fil
    model_ndx = 1
    for model in models:



# ============================================================================ #


def resortAtoms_name(template, target):
    """Sort atoms of a target structure based on template structure.
    Resorting is entirely based on atom names which are assumed to be unique!"""
    sortedAtoms = []
    atom_ndx = 1
    for atom in template:
        # find corresponding atom in target if present
        for targetAtom in target:
            if atom.name == targetAtom.name:
                sortedAtoms.append(copy.deepcopy(targetAtom))
                sortedAtoms[-1].num = atom_ndx
                sortedAtoms[-1].resid = atom.resid
                atom_ndx += 1

    return sortedAtoms


# ============================================================================ #
