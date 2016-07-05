#!/usr/bin/env python

from __future__ import print_function, division

import sys, os
import numpy as np
#import matplotlib.pyplot as plt
import scipy.constants as cn


def data2fes(data, T=298, bins=10, range=None):
    """
    Compute free energy from arbitrary data in kcal/mol.
    data.shape = [npoints, nvariables]
    """

    H, edges = np.histogramdd(data, bins=bins, range=range, normed=False, weights=None)

    E = cn.k*T * np.log(H.max() / H)
    E[E == float('inf')] = float('nan')
    E = E * cn.N_A / (cn.kilo * cn.calorie)

    return E, edges
 
