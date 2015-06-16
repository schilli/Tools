#!/usr/bin/env python
# Usage:
# resplica_hist.py replica_temp.xvg

import sys
import numpy as np
import matplotlib.pyplot as plt

temptrj    = np.loadtxt(sys.argv[1])
timestamps = temptrj[:,0]
temptrj    = temptrj[:,1:]
nreps      = temptrj.shape[1]
repIDs     = range(nreps)
repIDarray = np.ones_like(temptrj) * (np.arange(nreps)+1)

# sort temptrj by mean temperature
meantemp   = temptrj.mean(0)
meantempsrt_ndx = np.argsort(meantemp)
temptrj = temptrj[:,meantempsrt_ndx]

# plot temperature distribution
plt.hist2d(repIDarray.T.flatten(), temptrj.T.flatten(), bins=(nreps, nreps))
plt.xlabel("Replica ID")
plt.ylabel("Temperature ID")
plt.title("HREX - Temperature distribution")
plt.colorbar()
#plt.savefig("plots/HREX_temp_hist2d.png")
plt.show()
