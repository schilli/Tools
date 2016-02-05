#!/usr/bin/env python

from __future__ import print_function

import sys
import numpy  as np
import mdtraj as md

infile = sys.argv[1]

trj    = md.load(infile)
masses = np.array([atom.element.mass for atom in trj.top.atoms])

rg = md.compute_rg(trj, masses=masses)

print("{:.2f} nm".format(rg[0]))
