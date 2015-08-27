# general Lipari Szabo model fits

from __future__ import print_function

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ============================================================================ #

class LS(object):
    """
    Fit data to the general Lipari Szabo model.
    """

    # ================================ #

    def __init__(self, t, C, sigma=None):
        """
        Parameters
        ----------
        t : (N, ) array
            timeseries
        C : (M, N) array
            single correlation function or set of correlation functions
        """

        self.t = t 
        self.C = C
        if sigma is not None:
            self.sigma = sigma
        else:
            self.sigma = np.ones_like(C)
        self.internal = False

    # ================================ #

    def fit(self, n, fast=True, internal=False):
        """
        Fit the data and return optimal parameters
        """

        nS   = n
        ntau = n
        self.internal = internal
        if self.internal:
            ntau -= 1

        # set initial parameters
        p0  = list(np.linspace(0.1, 0.9, nS))
        p0 += ntau*[1.0]

        # set bounds on parameters
        bounds  = nS   * [(0,1)]    # S
        bounds += ntau * [(0,None)] # tau

        # set inequality constraints on S parameters
        constraints = []
        for idx in range(nS-1):
            constraints += [{"type": "ineq", "fun":  lambda x,i=idx: x[i+1] - x[i]}]  # Si+1 > Si

        # set inequality constraints on tau parameters
        for idx in range(nS, nS+ntau-1):
            constraints += [{"type": "ineq", "fun":  lambda x,i=idx: x[i+1] - x[i]}]  # tau_{i+1} > tau_{i}

        # if not fast, constrain Sf to 1.0
        if not fast:
            constraints += [{"type": "ineq", "fun":  lambda x,i=nS: x[i-1] - 1.0}] # constrain Sf to 1.0

        # fit
        res = minimize(self.generalLS_obj, x0=p0, bounds=bounds, constraints=constraints)
        
        # extract optimal parameters
        result = {}
        result["p"]    = res.x # all parameters for function input
        result["S"]    = res.x[:nS]
        result["tau"]  = res.x[nS:]
        if self.internal:
            result["tau"] = np.array(list(result["tau"]) + [float('inf')])

        return result

    # ================================ #

    def fit_global(self, n, fast=True):

        self.nC = self.C.shape[0]

        # set initial parameters (n S parameters, n-1 tau parameters and 1 global rotation time)
        p0 = np.ones((self.nC * (2*n-1)) + 1)

        # set bounds
        bounds = []
        for ic in range(self.nC):
            bounds +=  n   * [(0,1)]
            bounds += (n-1)* [(0,None)]
        bounds += [(0,None)]

        # set constraints
        nconst = 2*n-3 # constraints per correlation function
        constraints = []
        for ic in range(self.nC):
            for idx in range(n-1):
                constraints += [{"type": "ineq", "fun":  lambda x,i=ic*nconst+idx: x[i+1] - x[i]}]  # Si+1 > Si
            for idx in range(n, nconst-1):
                constraints += [{"type": "ineq", "fun":  lambda x,i=ic*nconst+idx: x[i+1] - x[i]}]  # tau_i+1 > tau_i

        if not fast:
            for ic in range(self.nC):
                constraints += [{"type": "ineq", "fun":  lambda x,i=ic*nconst+n-1: x[i] - 1.0}]

        # fit
        res = minimize(self.generalLS_global_obj, x0=p0, bounds=bounds, constraints=constraints)

        # extract optimal parameters
        result = {}

        return res.x

    # ================================ #

    def generalLS_internal(self, p):
        """
        Evaluate the general internal Lipari Szabo model.
        """

        n   = len(p)+1
        S   = p[:n/2]
        tau = p[n/2:] 

        C_I = S[0]
        for i, taui in enumerate(tau):
            C_I += (S[i+1] - S[i]) * np.exp(-self.t/taui)
        return C_I

    # ================================ #

    def generalLS(self, p):
        """
        General LS model
        
        Parameters
        ----------
        t : (N, ) array
            timeseries
        p : (2*n, ) array
        """

        if self.internal:
            C_I = self.generalLS_internal(p)
            return C_I
        else:
            C_I = self.generalLS_internal(p[:-1])
            C   = C_I * np.exp(-self.t/p[-1])
            return C
        
    # ================================ #

    def generalLS_obj(self, p):

        C          = self.generalLS(p)
        diff       = C - self.C
        squared    = diff**2
        weighted   = squared / self.sigma**2
        sumsquares = weighted.sum()
        return sumsquares

    # ================================ #

    def generalLS_global_obj(self, p):

        sumsquares  = 0.0
        pperC       = (len(p)-1) / self.C.shape[0] # parameters per correlation function

        for ic in range(self.nC):
            para        = list(p[ic*pperC:(ic+1)*pperC]) + [p[-1]]
            C           = self.generalLS(para)
            diff        = C - self.C[ic,:]
            squared     = diff**2
            weighted    = squared / self.sigma[ic,:]**2
            sumsquares += weighted.sum()

        print(sumsquares, p[-1])
        return sumsquares

 
    # ================================ #

