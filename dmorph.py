from __future__ import print_function

import sys
import numpy as np


def dmorph(firststructure, laststructure, nframes):
    """
    Linearly interpolate the dihedral angels from firststructure to laststructure.
    Input:
        firststructure: Starting structure for morph. Single frame mdtraj trajectory.
                        Will contain the final trajectory upon exit
        laststructure:  Last structure for morph. Single frame mdtraj trajectory. 
        nframes:        Number of frames for morph.
    """

    firststructure.n_frames         = nframes
    firststructure.time             = np.arange(nframe, dtype=np.float32)
    firststructure.timestep         = 1.0
    firststructure.unitcell_angles  = np.tile(firststructure.unitcell_angels[0,:],    (nframes, 1)   )
    firststructure.unitcell_lengths = np.tile(firststructure.unitcell_lengths[0,:],   (nframes, 1)   )
    firststructure.unitcell_vectors = np.tile(firststructure.unitcell_vectors[0,:,:], (nframes, 1, 1))
    firststructure.unitcell_volumes = np.tile(firststructure.unitcell_volumes[0],     (nframes)      )
    firststructure.xyz              = np.tile(firststructure.xyz[0,:,:],              (nframes, 1, 1))

    return firststructure



