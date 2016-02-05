#!/usr/bin/env python

from __future__ import print_function

import sys
import numpy  as np
import mdtraj as md

infile = sys.argv[1]


trj    = md.load(infile) 
dmax   = 0.0
for atom_ndx in range(trj.n_atoms):
    diff = trj.xyz[0,:,:] - trj.xyz[0,atom_ndx,:]
    dist = (diff**2).sum(1)**0.5
    distmax = dist.max()
    if distmax > dmax:
        dmax = distmax

print("{:.2f} nm".format(dmax))
