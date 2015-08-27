# -*- coding: UTF-8 -*-
# Compute NMR observables from MD trajectories

from __future__ import print_function

import sys, os, time, string, copy
import mdtraj
import numpy as np
import cPickle as pickle
import zipfile, bz2
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from mpi4py import MPI

import warnings


# ============================================================================ #


class corrFunction(object):
    """
    correlation function data, additional information and manipulation methods

    Parameters
    ----------
    corr : (nvec, nframes) array
        Correlation functions

    std : (nvec, nframes) array
        Correlation function standard deviations

    error : (nvec, nframes) array
        Correlation function standard error of the mean

    info : dict
        Dictionary with information on correlation functions
    """

    def __init__(self, corr=None, std=None, error=None, info=None):
        self.corr        = corr
        self.std         = std
        self.error       = error

        if info is not None:
            self.resid       = info['bondvecinfo']['resid'     ]
            self.resindex    = info['bondvecinfo']['resindex'  ]
            self.resname     = info['bondvecinfo']['resnames'  ]
            self.atomindex   = info['bondvecinfo']['atomindex' ]
            self.atomname    = info['bondvecinfo']['atomnames' ]
            self.element     = info['bondvecinfo']['element'   ]
            self.chain       = info['bondvecinfo']['chain'     ]
            self.bondlength  = info['bondvecinfo']['bondlength']
            self.bondvec     = info['bondvecinfo']['bondvec'   ]
            self.fitgroup    = info['bondvecinfo']['fitgroup'  ]
            try:
                self.fit     = info['bondvecinfo']['fit'       ]
            except KeyError:
                self.fit     = False
            try:
                self.S2direct = np.array(info['bondvecinfo']['S2'])
            except KeyError:
                self.S2direct = None
            self.dt          = info['bondvecinfo']['dt'        ]
            self.topfilename = info['topfilename']
            self.npzfilename = info['npzfilename']
            self.trjfilename = info['trjfilename'] 
            self.frames      = info['frames'     ]
        else:
            self.resid       = None
            self.resindex    = None
            self.resname     = None
            self.atomindex   = None
            self.atomname    = None
            self.element     = None
            self.chain       = None
            self.bondlength  = None
            self.bondvec     = None
            self.fitgroup    = None
            self.fit         = None
            self.dt          = None
            self.topfilename = None
            self.npzfilename = None
            self.trjfilename = None
            self.frames      = None


# ============================================================================ #


class OrderParameter(object):
    """
    Bond vector order parameters.

    Parameters
    ----------
    corrfilenames : list
        List with *.zip filenames containing correlation functions.

    converged : boolean, optional, default: True
        Use only converged correlation functions for averaging

    verbose : boolean, optional
        Report progress and other verbose output.

    **kwargs : optional keyword arguments
        All remaining keyword arguments are passed to the order parameter method.
    """

# ==================================== #

    def __init__(self, corrfilenames=None, converged=True, verbose=True, **kwargs):
        self.corrfilenames = corrfilenames
        self.converged     = converged
        self.verbose       = verbose

        self.corrlist      = []

        self.method_kwargs = kwargs

        # mean values
        self.avgcorr = None

        # for plotting
        self.figure  = None
        self.axs     = None
        self.cids    = {}   # connection ids for plot event handlers
        self.lines   = []   # line object handels
        self.corrset = 0    # which set of correlation functions from the corrlist should be plottet
        self.corridx = 0    # which correlatin function to plot next
        self.leftButtonPressed = False
        self.plotfit = None

        # load data
        if self.corrfilenames is not None:
            self.load()
 
# ==================================== #

    def load(self, corrfilenames=None):
        """
        Load sets of correlation functions.

        Parameters
        ----------
        corrfilenames : list
            List with *.zip filenames containing correlation functions.
 
        """

        if corrfilenames is not None:
            self.corrfilenames = corrfilenames

        starttime = time.time()

        # load correlation functions
        self.corrlist = []
        print("Loading {} set{} of correlation functions:".format(len(self.corrfilenames), *['s' if i>1 else '' for i in [len(self.corrfilenames)]]))
        for nf, filename in enumerate(self.corrfilenames):
            corr, corrstd, corrstdmean, info = load_corr(filename)
            self.corrlist.append(corrFunction(corr, corrstd, corrstdmean, info))
            if self.verbose:
                print("\rProgress: {:3.0f}%".format(100.0*nf/len(self.corrfilenames)), end="")
                sys.stdout.flush()

        # compute mean
        self.average_corr()

        # report runtime
        if self.verbose:
            print("\rLoading took: {:.2f} sec.".format(time.time() - starttime))

# ==================================== #

    def average_corr(self):
        """
        Compute the average correlation function and its standard deviation
        """
        nvec, nframes = self.corrlist[0].corr.shape
        ncorr         = len(self.corrlist)
        allcorr = np.zeros([ncorr, nvec, nframes], dtype=self.corrlist[0].corr.dtype)

        for c, corr in enumerate(self.corrlist):
            allcorr[c,:,:] = corr.corr

        self.avgcorr = corrFunction(corr=allcorr.mean(0), std=allcorr.std(0), error=allcorr.std(0)/allcorr.shape[0]**0.5)
        self.avgcorr.resid       = copy.copy(self.corrlist[0].resid      )
        self.avgcorr.resindex    = copy.copy(self.corrlist[0].resindex   )
        self.avgcorr.resname     = copy.copy(self.corrlist[0].resname    )
        self.avgcorr.atomindex   = copy.copy(self.corrlist[0].atomindex  )
        self.avgcorr.atomname    = copy.copy(self.corrlist[0].atomname   )
        self.avgcorr.element     = copy.copy(self.corrlist[0].element    )
        self.avgcorr.chain       = copy.copy(self.corrlist[0].chain      )
        self.avgcorr.bondlength  = copy.copy(self.corrlist[0].bondlength )
        self.avgcorr.bondvec     = copy.copy(self.corrlist[0].bondvec    )
        self.avgcorr.fitgroup    = copy.copy(self.corrlist[0].fitgroup   )
        self.avgcorr.fit         = copy.copy(self.corrlist[0].fit        )
        self.avgcorr.dt          = copy.copy(self.corrlist[0].dt         )
        self.avgcorr.topfilename = copy.copy(self.corrlist[0].topfilename)
        self.avgcorr.npzfilename = None
        self.avgcorr.trjfilename = None
        self.avgcorr.frames      = None

# ==================================== #

    def estimate(self, method="mean", converged=True, **kwargs):
        """
        Estimate bond vector order parameters from correlation functions.

        Parameters
        ----------
        method: string, optional
            The method to use for order parameter computation.
            Options are:
                "direct"        If provided, use S2 values directly approximated from the bond vectors as described in: 
                                Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
                "mean"          Use the mean of the final quarter as order parameter
                "single exp"    Fit correlation functions to single exponential
                "double exp"    Fit correlation functions to double exponential
                "LS"            Use the Lipari Szabo model to estimate S2 values
                "extLS"         Use the extended Lipari Szabo method (method 3) from:
                                JPCB 2008, 112, 6155-6158, pubs.acs.org/doi/abs/10.1021/jp077018h

        converged : boolean, optional, default: True
            Use only converged correlation functions for averaging
        """
 
        self.method    = method
        self.converged = converged

        # select order parameter estimation method
        if self.method == "direct":
            self.estimate_direct(**kwargs)
        elif self.method == "mean":
            self.estimate_mean(converged=self.converged, **kwargs)
        elif self.method == "single exp":
            self.estimate_single_exp(**kwargs)
        elif self.method == "double exp":
            self.estimate_double_exp(**kwargs) 
        elif self.method == "LS":
            self.estimate_LS(**kwargs)   
        elif self.method == "extLS":
            self.estimate_extLS(**kwargs)  
        elif self.method == "CloreLS":
            self.estimate_CloreLS(**kwargs)   
        elif self.method == "extCloreLS":
            self.estimate_extCloreLS(**kwargs)    
        else:
            print("Order parameter estimation method unknown: {}".format(self.method))
            sys.exit(1)
 
# ==================================== #

    def estimate_direct(self):
        """
        Compute mean bond vector order parameters from direct estimates as described in:
        Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
        """
        if self.corrlist[0].S2direct is not None:
            self.S2all = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            for n, corrfun in enumerate(self.corrlist):
                self.S2all[:,n] = corrfun.S2direct
            self.S2mean  = self.S2all.mean(1)
            self.S2std   = self.S2all.std(1)
            self.S2error = self.S2all.std(1) / self.S2all.shape[1]**0.5 
        else:
            print("Direct estimate of S2 not possible. Data not present in correlation functions.")

# ==================================== #

    def estimate_mean(self, converged=True, diffThreshold=0.02, stdThreshold=0.02):
        """
        Estimate bond vector order parameters as the mean of the last quarter of the correlation function.

        Parameters
        ----------
        converged : boolean, optional
            If True, only use converged correlation functions.
        diffThreshold : float, optional
            If the quarters 3 and 4 of the correlation function differ more than this threshold,
            they are not considered converged
        stdThreshold : float, optional
            If the mean standard deviation of quarter 3 and 4 is larger than this threshold,
            they are not considered converged.
        """

        self.S2all         = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
        self.S2convergence = np.zeros_like(self.S2all, dtype=np.bool)

        # compute S2 values for each correlation function and judge convergence
        for c, corrfun in enumerate(self.corrlist):
            mean, convergence = self.estimate_mean_single(corrfun, diffThreshold=diffThreshold, stdThreshold=stdThreshold)
            self.S2all[:,c]         = mean
            self.S2convergence[:,c] = convergence

        # compute global S2 as mean of individiual S2
        if not converged:
            self.S2mean  = self.S2all.mean(1)
            self.S2std   = self.S2all.std(1)
            self.S2error = self.S2all.std(1) / self.S2all.shape[1]**0.5
        else:
            self.S2mean  = np.zeros(self.S2all.shape[0])
            self.S2std   = np.zeros(self.S2all.shape[0])
            self.S2error = np.zeros(self.S2all.shape[0])
            for n in range(self.S2all.shape[0]):
                self.S2mean[n]  = self.S2all[:,self.S2convergence[n,:]].mean()
                self.S2std[n]   = self.S2all[:,self.S2convergence[n,:]].std()
                self.S2error[n] = self.S2std[n] / self.S2convergence[n,:].sum()

        self.S2nconverged = self.S2convergence.sum(1)

        # compute global S2 from average correlation functions
        self.S2avg, self.S2avgConvergence = self.estimate_mean_single(self.avgcorr, diffThreshold=diffThreshold, stdThreshold=stdThreshold)

# ==================================== #

    def estimate_mean_single(self, corrfun, converged=True, diffThreshold=0.02, stdThreshold=0.02):
        """
        Estimate bond vector order parameters for a single set of correlation functions.
        
        Parameters
        ----------
        corrfun : (nvec, nframes) array
            Single set of correlation functions
        other parameters:
            see function "estimate_mean()"
        """

        length            = corrfun.corr.shape[1]
        quarter           = length / 4
        thirdQuarter      = corrfun.corr[:,2*quarter:3*quarter]
        fourthQuarter     = corrfun.corr[:,3*quarter:4*quarter]
        fourthQuarterMean = fourthQuarter.mean(1)

        difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
        stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
        convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)

        return fourthQuarterMean, convergence 

# ==================================== #

    def estimate_single_exp(self, guess=[1.0, 1.0], avgonly=False):
        """
        Estimate bond vector order parameters from single exponential fits to the correlation functions.

        Parameters
        ----------
        guess : [a, t1]
            Initial guess of exponential coefficient and decay constant, defaults to 1's.  
        avgonly : boolean, optional
            Compute only the fit of the average 
        """

        # compute global S2 from average correlation functions
        self.S2avg, self.tauavg, self.poptavg, self.pvaravg = self.estimate_single_exp_single(self.avgcorr, guess=guess)

        if not avgonly:
            self.S2all  = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            self.tauall = np.zeros_like(self.S2all)
            self.popt   = np.zeros([self.corrlist[0].corr.shape[0], 2, len(self.corrlist)])
            self.pvar   = np.zeros([self.corrlist[0].corr.shape[0], 2, len(self.corrlist)])

            # compute S2 values for each correlation function
            for c, corrfun in enumerate(self.corrlist):
                print("\rProgress: {:.0f}%".format(100.0*(c+1)/len(self.corrlist)), end="")
                sys.stdout.flush()
                S2, tau, popt, pvar = self.estimate_single_exp_single(corrfun, guess=guess)
                self.S2all[:,c]  = S2
                self.tauall[:,c] = tau
                self.popt[:,:,c] = popt
                self.pvar[:,:,c] = pvar
            print("\r", 100*" ", "\r", end="")
     
            # compute global S2 and tau as mean of individiual S2 and tau
            self.S2mean   = self.S2all.mean(1)
            self.S2std    = self.S2all.std(1)
            self.S2error  = self.S2all.std(1) / self.S2all.shape[1]**0.5 
            self.taumean  = self.tauall.mean(1)
            self.taustd   = self.tauall.std(1)
            self.tauerror = self.tauall.std(1) / self.tauall.shape[1]**0.5  

# ==================================== #

    def estimate_single_exp_single(self, corrfun, guess=[1.0, 1.0]):
        """
        Estimate bond vector order parameters with a single exponential for a single set of correlation functions.

        Parameters
        ----------
        corrfun : corrFunction
            Correlation function object
        guess : [a, t1]
            Initial guess of exponential coefficient and decay constant, defaults to 1's. 

        Returns
        -------
        S2 : (nvec) array
            Array of S2 values per bond vector
        tau : (nvec) array
            Decay times of bond vector correlations
        popt : (nvec, 2) array
            fitted parameters for exponential fit
        """
        dt    = corrfun.dt
        xdata = np.linspace(0, dt*corrfun.corr.shape[1], corrfun.corr.shape[1])
        S2    = np.zeros(corrfun.corr.shape[0])
        tau   = np.zeros(corrfun.corr.shape[0])
        popt  = np.zeros([corrfun.corr.shape[0], 2])
        pvar  = np.zeros([corrfun.corr.shape[0], 2])

        bounds = ((0,1), (0,None))

        for n in range(corrfun.corr.shape[0]):
            res = minimize(n_exp_obj, x0=guess, args=(xdata, self.avgcorr.corr[n,:]), bounds=bounds)
            popt[n,:] = res.x  
            #try:

            #    popt[n,:], pcov = curve_fit(n_exp, xdata, corrfun.corr[n,:], p0=guess)
            #    pvar[n,:] = np.diag(pcov)
            #except RuntimeError:
            #    popt[n,:] = np.array([0.0, 0.0])
            #    pvar[n,:] = np.array([0.0, 0.0])

            S2[n]     = 1 - popt[n,0]
            tau[n]    =     popt[n,1]

        return S2, tau, popt, pvar

# ==================================== #

    def estimate_double_exp(self, guess=[1.0, 1.0, 1.0, 1.0], avgonly=False):
        """
        Estimate bond vector order parameters from double exponential fits to the correlation functions.

        Parameters
        ----------
        guess : [a, b, t1, t2], optional
            Initial guess of exponential coefficient and decay constant, defaults to 1's.  
        avgonly : boolean, optional
            Compute only the fit of the average
        """

        # compute global S2 from average correlation functions
        self.S2avg, self.tauavg, self.poptavg, self.pvaravg = self.estimate_double_exp_single(self.avgcorr, guess=guess)

        if not avgonly:
            self.S2all  = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            self.tauall = np.zeros_like(self.S2all)
            self.popt   = np.zeros([self.corrlist[0].corr.shape[0], 4, len(self.corrlist)])
            self.pvar   = np.zeros([self.corrlist[0].corr.shape[0], 4, len(self.corrlist)])

            # compute S2 values for each correlation function
            for c, corrfun in enumerate(self.corrlist):
                print("\rProgress: {:.0f}%".format(100.0*(c+1)/len(self.corrlist)), end="")
                sys.stdout.flush()
                S2, tau, popt, pvar = self.estimate_double_exp_single(corrfun, guess=guess)
                self.S2all[:,c]  = S2
                self.tauall[:,c] = tau
                self.popt[:,:,c] = popt
                self.pvar[:,:,c] = pvar
            print("\r", 100*" ", "\r", end="")
     
            # compute global S2 and tau as mean of individiual S2 and tau
            self.S2mean   = self.S2all.mean(1)
            self.S2std    = self.S2all.std(1)
            self.S2error  = self.S2all.std(1) / self.S2all.shape[1]**0.5 
            self.taumean  = self.tauall.mean(1)
            self.taustd   = self.tauall.std(1)
            self.tauerror = self.tauall.std(1) / self.tauall.shape[1]**0.5  

# ==================================== #

    def estimate_double_exp_single(self, corrfun, guess=[1.0, 1.0, 1.0, 1.0]):
        """
        Estimate bond vector order parameters with a double exponential for a single set of correlation functions.

        Parameters
        ----------
        corrfun : corrFunction
            Correlation function object
        guess : [a, b, t1, t2]
            Initial guess of exponential coefficient and decay constant, defaults to 1's. 

        Returns
        -------
        S2 : (nvec) array
            Array of S2 values per bond vector
        tau : (nvec) array
            Decay times of bond vector correlations
        popt : (nvec, 2) array
            fitted parameters for exponential fit
        """
        dt    = corrfun.dt
        xdata = np.linspace(0, dt*corrfun.corr.shape[1], corrfun.corr.shape[1])
        S2    = np.zeros(corrfun.corr.shape[0])
        tau   = np.zeros(corrfun.corr.shape[0])
        tau2  = np.zeros(corrfun.corr.shape[0])
        popt  = np.zeros([corrfun.corr.shape[0], 4])
        pvar  = np.zeros([corrfun.corr.shape[0], 4])

        nrandomguesses = 100

        for n in range(corrfun.corr.shape[0]):
            try:
                if self.poptavg is not None and self.poptavg.shape[1] == len(guess):
                    p0 = self.poptavg[n,:]
                else:
                    p0      = guess
            except AttributeError:
                p0      = guess

            counter = 0
            while True:
                counter += 1
                try:
                    popt[n,:], pcov = curve_fit(n_exp, xdata, corrfun.corr[n,:], p0=p0)
                    pvar[n,:] = np.diag(pcov)
                    # fit was successful if it didn't raise an error and has proper error estimates
                    if not float('NaN') in pvar[n,:]:
                        break
                    else:
                        print("curve_fit failed with NaN error estimates") 
                except RuntimeError:
                    # change initial parameter guess
                    if counter == 1 and n > 0:
                        # second try with optimized values from previous iteration
                        print("    trying previously optimized parameters as initial guess")
                        p0 = popt[n-1,:]
                    elif counter == 2 and n > 0:
                        # third try with mean of optimized values from previous iterations
                        print("    trying mean of previously optimized parameters as initial guess")
                        p0 = popt[n:-1,:].mean(0)
                    elif counter < 2 + nrandomguesses:
                        # continue trying random values
                        print("    trying random parameters as initial guess", counter)
                        p0 = np.random.rand(len(p0)) + 10**np.random.randint(-2, 4, size=len(p0))
                    else:
                        # failed
                        print("    failed to converge fitting")
                        popt[n,:] = float('Nan') * np.ones([len(p0)])
                        pvar[n,:] = float('Nan') * np.ones([len(p0)])
                        break
 

            coeffs = popt[n,:2]
            tconst = popt[n,2:]
            coeffsvar = pvar[n,:2]
            tconstvar = pvar[n,2:] 
            srtidx = np.argsort(tconst)
            popt[n,:2] = coeffs[srtidx]
            popt[n,2:] = tconst[srtidx]
            pvar[n,:2] = coeffsvar[srtidx]
            pvar[n,2:] = tconstvar[srtidx] 
            S2[n]     = 1 - popt[n,0]
            tau[n]    =     popt[n,2]

        return S2, tau, popt, pvar
 
# ==================================== #

    def estimate_LS(self, guess=[0.8, 10.0, 100.0], avgonly=False):

        # compute global S2 from average correlation functions
        self.S2avg, self.poptavg, self.pvaravg = self.estimate_LS_single(self.avgcorr, guess=guess)
 
        if not avgonly:
            self.S2all  = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            self.popt   = np.zeros([self.corrlist[0].corr.shape[0], len(guess), len(self.corrlist)])
            self.pvar   = np.zeros([self.corrlist[0].corr.shape[0], len(guess), len(self.corrlist)])

            # compute S2 values for each correlation function
            for c, corrfun in enumerate(self.corrlist):
                print("\rProgress: {:.0f}%".format(100.0*(c+1)/len(self.corrlist)), end="")
                sys.stdout.flush()
                S2, popt, pvar = self.estimate_LS_single(corrfun, guess=guess)
                self.S2all[:,c]  = S2
                self.popt[:,:,c] = popt
                self.pvar[:,:,c] = pvar
            print("\r", 100*" ", "\r", end="")
     
            # compute global S2 and tau as mean of individiual S2 and tau
            self.S2mean   = self.S2all.mean(1)
            self.S2std    = self.S2all.std(1)
            self.S2error  = self.S2all.std(1) / self.S2all.shape[1]**0.5 

# ==================================== #

    def estimate_LS_single(self, corrfun, guess=[0.8, 10.0, 100.0]):

        dt    = corrfun.dt
        xdata = np.linspace(0, dt*corrfun.corr.shape[1], corrfun.corr.shape[1])
        S2    = np.zeros(corrfun.corr.shape[0])
        popt  = np.zeros([corrfun.corr.shape[0], len(guess)])
        pvar  = np.zeros([corrfun.corr.shape[0], len(guess)])

        nrandomguesses = 100

        for n in range(corrfun.corr.shape[0]):
            try:
                if self.poptavg is not None and self.poptavg.shape[1] == len(guess):
                    p0 = self.poptavg[n,:]
                else:
                    p0      = guess
            except AttributeError:
                p0      = guess

            counter = 0
            while True:
                counter += 1
                try:
                    popt[n,:], pcov = curve_fit(LipariSzabo, xdata, corrfun.corr[n,:], p0=guess)
                    pvar[n,:] = np.diag(pcov)
                    # fit was successful if it didn't raise an error and has proper error estimates
                    if not float('NaN') in pvar[n,:]:
                        break
                    else:
                        print("curve_fit failed with NaN error estimates") 
                except RuntimeError:
                    # change initial parameter guess
                    if counter == 1 and n > 0:
                        # second try with optimized values from previous iteration
                        print("    trying previously optimized parameters as initial guess")
                        p0 = popt[n-1,:]
                    elif counter == 2 and n > 0:
                        # third try with mean of optimized values from previous iterations
                        print("    trying mean of previously optimized parameters as initial guess")
                        p0 = popt[n:-1,:].mean(0)
                    elif counter < 2 + nrandomguesses:
                        # continue trying random values
                        print("    trying random parameters as initial guess", counter)
                        p0 = np.random.rand(len(p0)) + 10**np.random.randint(-2, 4, size=len(p0))
                    else:
                        # failed
                        print("    failed to converge fitting")
                        popt[n,:] = float('Nan') * np.ones([len(p0)])
                        pvar[n,:] = float('Nan') * np.ones([len(p0)])
                        break 

        S2 = popt[:,0]
        return S2, popt, pvar

# ==================================== #

    def estimate_extLS(self, guess=[0.8, 0.9, 10.0, 100.0], avgonly=False, modelSelection=False):

        # compute global S2 from average correlation functions
        self.S2avg, self.poptavg, self.pvaravg = self.estimate_extLS_single(self.avgcorr, guess=guess, modelSelection=modelSelection)
 
        if not avgonly:
            self.S2all  = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            self.popt   = np.zeros([self.corrlist[0].corr.shape[0], len(guess), len(self.corrlist)])
            self.pvar   = np.zeros([self.corrlist[0].corr.shape[0], len(guess), len(self.corrlist)])

            # compute S2 values for each correlation function
            for c, corrfun in enumerate(self.corrlist):
                print("\rProgress: {:.0f}%".format(100.0*(c+1)/len(self.corrlist)), end="")
                sys.stdout.flush()
                S2, popt, pvar = self.estimate_extLS_single(corrfun, guess=guess, modelSelection=modelSelection)
                self.S2all[:,c]  = S2
                self.popt[:,:,c] = popt
                self.pvar[:,:,c] = pvar
            print("\r", 100*" ", "\r", end="")
     
            # compute global S2 and tau as mean of individiual S2 and tau
            self.S2mean   = self.S2all.mean(1)
            self.S2std    = self.S2all.std(1)
            self.S2error  = self.S2all.std(1) / self.S2all.shape[1]**0.5 

# ==================================== #

    def estimate_extLS_single(self, corrfun, guess=[0.8, 0.9, 10.0, 100.0], modelSelection=False):

        dt    = corrfun.dt
        xdata = np.linspace(0, dt*corrfun.corr.shape[1], corrfun.corr.shape[1])
        S2    = np.zeros(corrfun.corr.shape[0])
        popt  = np.zeros([corrfun.corr.shape[0], len(guess)])
        pvar  = np.zeros([corrfun.corr.shape[0], len(guess)])

        nrandomguesses = 0

        for n in range(corrfun.corr.shape[0]):
            try:
                if self.poptavg is not None and self.poptavg.shape[1] == len(guess):
                    p0 = self.poptavg[n,:]
                else:
                    p0      = guess
            except AttributeError:
                p0      = guess

            counter = 0
            while True:
                counter += 1
                try:
                    # do the fit without the first point (t=0)
                    weights = np.linspace(1, 100, xdata.shape[0]-1)
                    popt[n,:], pcov = curve_fit(extLipariSzabo, xdata[1:], corrfun.corr[n,1:], p0=p0)
                    pvar[n,:] = np.diag(pcov)
                    # fit was successful if it didn't raise an error and has proper parameter estimates
                    if np.isnan(popt[n,:]).sum() == 0 and np.isinf(pvar[n,:]).sum() == 0 and popt[n,1] > 0 and popt[n,1] <= 1.0:
                        break
                    else:
                        #print("curve_fit failed with NaN error estimates for n = {}".format(n)) 
                        raise RuntimeError
                except RuntimeError:
                    # change initial parameter guess
                    if counter == 1:
                        # second try with all ones as initial parameter guess
                        p0 = len(guess) * [1.0]
                        print("    trying all ones as initial guess")
                    elif counter == 2 and n > 0:
                        # third try with optimized values from previous iteration
                        print("    trying previously optimized parameters as initial guess")
                        p0 = popt[n-1,:]
                    elif counter == 3 and n > 0:
                        # fourth try with mean of optimized values from previous iterations
                        print("    trying mean of previously optimized parameters as initial guess")
                        p0 = popt[n:-1,:].mean(0)
                    elif counter < 3 + nrandomguesses:
                        # continue trying random values
                        print("    trying random parameters as initial guess", counter)
                        p0 = np.random.rand(len(p0)) + 10**np.random.randint(-2, 4, size=len(p0))
                    else:
                        # failed
                        print("    failed to converge fitting for n = {}".format(n))
                        print("    trying normal LS instead")
                        popt[n,[0,2,3]], pcov = curve_fit(LipariSzabo, xdata, corrfun.corr[n,:])
                        pvar[n,[0,2,3]] = np.diag(pcov)
                        popt[n,1] = 1.0
                        pvar[n,1] = 0.0
                        popt[n,2] = 0.0 # to trigger modelSelection
                        #popt[n,:] = float('Nan') * np.ones([len(p0)])
                        #pvar[n,:] = float('Nan') * np.ones([len(p0)])
                        break 


            if modelSelection:
                if popt[n,2] < 1.0 or popt[n,2] > 40:
                    p, cov = curve_fit(n_exp, xdata, corrfun.corr[n,:], p0=[0.9, 1000.0])
                    popt[n,:] = np.array([p[0],     p[0],     1.0, p[1]    ])
                    pvar[n,:] = np.array([cov[0,0], cov[0,0], 0.0, cov[1,1]])



        S2 = popt[:,0]
        return S2, popt, pvar
 
# ==================================== #

    def estimate_CloreLS(self, guess=[0.8, 0.9, 10.0, 100.0, 1000.0]):

        dt    = self.avgcorr.dt
        xdata = np.linspace(0, dt*self.avgcorr.corr.shape[1], self.avgcorr.corr.shape[1])
        S2    = np.zeros(self.avgcorr.corr.shape[0])
        popt  = np.zeros([self.avgcorr.corr.shape[0], len(guess)])
        pvar  = np.zeros([self.avgcorr.corr.shape[0], len(guess)])

        for n in range(self.avgcorr.corr.shape[0]):
            res = minimize(CloreLS_obj, x0=guess, args=(xdata, self.avgcorr.corr[n,:]), bounds=((0,1), (0,1), (1,1e4), (1,1e4), (1,1e4)))
            popt[n,:] = res.x

        # compute global S2 from average correlation functions
        self.S2avg   = popt[:,0]
        self.poptavg = popt

# ==================================== #

    def estimate_extCloreLS(self, guess=[0.8, 0.9, 0.95, 10.0, 100.0, 1000.0]):

        dt    = self.avgcorr.dt
        xdata = np.linspace(0, dt*self.avgcorr.corr.shape[1], self.avgcorr.corr.shape[1])
        S2    = np.zeros(self.avgcorr.corr.shape[0])
        popt  = np.zeros([self.avgcorr.corr.shape[0], len(guess)])
        pvar  = np.zeros([self.avgcorr.corr.shape[0], len(guess)])


        #bounds = ((0,1), (0,1), (0,1), (1,100), (1,500), (1,1e4))
        bounds = ((0,1), (0,1), (0,1), (1,None), (1,None), (1,None))
        constraints = [{"type": "ineq", "fun":  lambda x: x[1] - x[0]},
                       {"type": "ineq", "fun":  lambda x: x[2] - x[1]}]

        for n in range(self.avgcorr.corr.shape[0]):
            res = minimize(extCloreLS_obj, x0=guess, args=(xdata[1:], self.avgcorr.corr[n,1:]), bounds=bounds, constraints=constraints)
            popt[n,:] = res.x

        # compute global S2 from average correlation functions
        self.S2avg   = popt[:,0]
        self.poptavg = popt
 
# ==================================== #

    def plot_corr(self, event=None, corrset=None, fit=None):
        """
        Plot correlation functions.
        
        Parameters
        ----------
        event : matplotlib mouse event
            Only there to update the plot on various matplotlib events, like mouse button presses
        corrset: integer, optional
            Which set of correlation functions to use for plotting.
            If None, plot the mean.
        fit: string
            What fit should be plotted? Options are:
                "double exp"
                "LS"
        """

        if fit is not None:
            self.plotfit = fit

        if self.figure is None:
            # create figure and axis
            self.figure = plt.figure()
            self.axs    = self.figure.add_subplot(111)  
            self.cids['button_press'  ] = self.figure.canvas.mpl_connect('button_press_event',   self._onclick)
            self.cids['button_release'] = self.figure.canvas.mpl_connect('button_release_event', self._onrelease)
            self.cids['motion'        ] = self.figure.canvas.mpl_connect('motion_notify_event',  self._onmove)
            self.cids['scroll'        ] = self.figure.canvas.mpl_connect('scroll_event',         self._onscroll)
            self.cids['close'         ] = self.figure.canvas.mpl_connect('close_event',          self._onclose)

            if corrset is None:
                self.corrset  = 0
                self.plotMean = True
            else:
                self.corrset  = corrset
                self.plotMean = False
            self.corridx = 0
            self.lines   = []

            if self.avgcorr is None:
                self.average_corr()

        # remove old data
        while len(self.lines) > 0:
            line = self.lines.pop()
            line.remove() 

        start_idx = 0
#        if self.plotfit in ["extLS"]:
#            start_idx = 1

        # plot data
        if self.plotMean:
            corrFun     = self.avgcorr
            xdata       = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, corrFun.corr.shape[1])    
            self.lines += self.axs.plot(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:], 'b', label="Mean correlation function over {} samples".format(len(self.corrlist)), lw=2)
            #self.lines.append(self.axs.fill_between(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:]+self.avgcorr.std[self.corridx,:],
            #                                               self.avgcorr.corr[self.corridx,start_idx:]-self.avgcorr.std[self.corridx,:],
            #                                               alpha=0.4, color='r')) 
            self.lines.append(self.axs.fill_between(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:]+1.96*self.avgcorr.error[self.corridx,start_idx:],
                                                           self.avgcorr.corr[self.corridx,start_idx:]-1.96*self.avgcorr.error[self.corridx,start_idx:],
                                                           alpha=0.4, color='r', label="95% Confidence Interval")) 

            try:
                coeffs = ",  ".join(["{:.1f}".format(n) for n in self.poptavg[self.corridx,:]])
            except AttributeError as e:
                print("Did you forget to fit the data?")
                raise e
            xdata  = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, 10*corrFun.corr.shape[1])  
            if self.plotfit in ["single exp", "double exp"]:
                self.lines += self.axs.plot(xdata[start_idx:], n_exp(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "LS":
                self.lines += self.axs.plot(xdata[start_idx:], LipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "extLS":
                self.lines += self.axs.plot(xdata[start_idx:], extLipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "CloreLS":
                coeffs = ",  ".join(["{:.1f}".format(n) for n in self.poptavg[self.corridx,:]])
                self.lines += self.axs.plot(xdata[start_idx:], CloreLS(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2) 
            elif self.plotfit == "extCloreLS":
                label = r"fit: S$^2$={:.2f},  S$_f$$^2$={:.2f},  S$_i$$^2$={:.2f},  $\tau_f$={:.1f},  $\tau_e$={:.1f},  $\tau_m$={:.1f}".format(*self.poptavg[self.corridx,:])
                self.lines += self.axs.plot(xdata[start_idx:], extCloreLS(xdata[start_idx:], *self.poptavg[self.corridx,:]), label=label, color='g', lw=2)  

        else:
            corrFun     = self.corrlist[self.corrset]
            xdata       = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, corrFun.corr.shape[1])   
            self.lines += self.axs.plot(xdata[start_idx:], corrFun.corr[self.corridx,start_idx:], 'b', lw=2)
            try:
                coeffs = " ".join(["{:.1e}".format(n) for n in self.poptavg[self.corridx,:]])
            except AttributeError as e:
                print("Did you forget to fit the data?")
                raise e 
            xdata  = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, 10*corrFun.corr.shape[1])  
            if self.plotfit in ["single exp", "double exp"]:
                self.lines += self.axs.plot(xdata[start_idx:], n_exp(xdata[start_idx:], *self.poptavg[self.corridx,start_idx:, self.corrset]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "LS":
                self.lines += self.axs.plot(xdata[start_idx:], LipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:, self.corrset]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "extLS":
                self.lines += self.axs.plot(xdata[start_idx:], extLipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2) 

        # set axis limits
        self.axs.set_ylim(min(0, corrFun.corr.min()), max(1, corrFun.corr.max()))

        # plot scrollbar
        #xmin, xmax = self.axs.get_xlim()
        #self.lines += self.axs.plot([xmin, xmax], 2*[0.98], 'k', linewidth=2)
        #self.lines += self.axs.plot(xmin + self.corridx*(xmax-xmin)/(corrFun.corr.shape[0]-1), 0.98, 'sk', markersize=15)

        # annotate plot
        self.axs.set_title("{} {} ({}-{})".format(corrFun.resname[0][self.corridx], corrFun.resid[0][self.corridx], corrFun.atomname[0][self.corridx], corrFun.atomname[1][self.corridx]))
        self.axs.set_ylabel("correlation")
        self.axs.set_xlabel("time [ps]")
        self.axs.legend(loc="upper right")
        #self.axs.legend(loc="lower left")

        self.figure.canvas.draw()
        #plt.savefig("plots/extCloreLS/corr_{:03d}.jpg".format(corrFun.resid[0][self.corridx]), dpi=100, papertype='a5', orientation='landscape', bbox_inches='tight')
        plt.show()
 
# ==================================== #

    def _onclick(self, event):
        if event.button == 1:
            self.leftButtonPressed = True
            try:
                xmin, xmax = self.axs.get_xlim()
                try:
                    ncorr = self.corrlist[self.corrset].corr.shape[0]
                except IndexError:
                    ncorr = self.avgcorr.corr.shape[0] 
                self.corridx = int(np.round((ncorr-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.corridx < 0:
                    self.corridx = 0
                elif self.corridx >= ncorr:
                    self.corridx = ncorr - 1
                self.plot_corr()
            except:
                pass 

# ==================================== #

    def _onrelease(self, event):
        if event.button == 1:
            self.leftButtonPressed = False 

# ==================================== #

    def _onmove(self, event):
        if self.leftButtonPressed:
            try:
                xmin, xmax = self.axs.get_xlim()
                try:
                    ncorr = self.corrlist[self.corrset].corr.shape[0]
                except IndexError:
                    ncorr = self.avgcorr.corr.shape[0]
                self.corridx = int(np.round((ncorr-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.corridx < 0:
                    self.corridx = 0
                elif self.corridx >= ncorr:
                    self.corridx = ncorr - 1 
                self.plot_corr()
            except:
                pass 

# ==================================== #

    def _onscroll(self, event):
        try:
            ncorr = self.corrlist[self.corrset].corr.shape[0]
        except IndexError:
            ncorr = self.avgcorr.corr.shape[0]
        self.corridx -= int(event.step)
        if self.corridx < 0:
            self.corridx = ncorr - 1
        elif self.corridx >= ncorr:
            self.corridx = 0

        self.plot_corr()
 
# ==================================== #
    
    def _onclose(self, event):
        self.figure  = None
        self.axs     = None
        self.cids    = {}
        self.lines   = []

# ============================================================================ #
 

def _fit_trj(trj, fitframe=0, ref=None, fitgroup="name CA", parallel=True):
    """
    Fit MD trajectory to reference based on the specified set of atoms.
    Fitting is done in place.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to fit.

    fitframe : int, optional
        Frame index of reference frame.

    ref : mdtraj trajectory
        Reference MDTRAJ trajectory.
        If none, the first trajectory is used as a reference.

    fitgroup : string or integer array
        mdtraj selection string or integer array of atom indices
        to select a set of atoms for RMSD fit.

    parallel : boolean, optional
        Toggle parallel fitting algorithm of mdtraj
    """
    if type(ref) == type(None):
        ref = trj

    try:
        fit_atomndx = trj.topology.select(fitgroup)
    except:
        fit_atomndx = fitgroup

    trj.superpose(ref, fitframe, atom_indices=fit_atomndx, parallel=parallel)
 

# ============================================================================ #

def _vec_atoms(trj, bondvec=None):
    """
    Compute bond atom indices.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory for topology information.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    Returns
    -------
    bondvec_ndx : (nvec, 2) array
        Integer array of all atom indices for bond vector atoms
    """

    # default option
    if type(bondvec) == type(None):
        selectiontext1   = "name N and not resname PRO and resid >= 1"
        selectiontext2   = "name H and not resname PRO and resid >= 1"
        atom1_ndx        = trj.topology.select(selectiontext1)
        atom2_ndx        = trj.topology.select(selectiontext2) 

    else:
        sel1 = bondvec[0]
        sel2 = bondvec[1]
        # if bond vectors are already specified with integer arrays
        try:
            bondvec_ndx = np.zeros([len(sel1), 2], dtype=np.int)
            bondvec_ndx[:,0] = np.array(sel1)
            bondvec_ndx[:,1] = np.array(sel2)
            return bondvec_ndx
        # if bond vectors are specified with selection texts
        except:
            atom1_ndx        = trj.topology.select(sel1)
            atom2_ndx        = trj.topology.select(sel2)  

    bondvec_ndx      = np.zeros([atom1_ndx.shape[0], 2], dtype=np.int)
    bondvec_ndx[:,0] = atom1_ndx
    bondvec_ndx[:,1] = atom2_ndx

    return bondvec_ndx


# ============================================================================ #


def _bond_vec(trj, bondvec_ndx):
    """
    Compute bond vectors from trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory

    bondvec_ndx : (nvec, 2) array
        Bond vector indices

    Returns
    -------
    bondvec : (nframes, nvec, 3) array
        Array of normalized bond vectors
    info : dict
        Dictionary with information on the bondvectors.
        Keys are:
        'dt'            time step in pico seconds
        'bondlength'    one bondlength per bond vector
        'resnames'      two sequences of resnames, one for each atom in bondvector
        'resids'                         resids
        'resindex'                       resindices
        'atomnames'                      atomnames
        'atomindex'                      atomindices
        'element'                        elements
        'chain'                          chains ids
    """

    # extract bond vector subtrajectories
    atom1_xyz = trj.xyz[:,bondvec_ndx[:,0],:]
    atom2_xyz = trj.xyz[:,bondvec_ndx[:,1],:]
    #atom1_trj = trj.atom_slice(bondvec_ndx[:,0])
    #atom2_trj = trj.atom_slice(bondvec_ndx[:,1])

    # compute bond vectors
    bondvec = atom1_xyz - atom2_xyz
    #bondvec = atom1_trj.xyz - atom2_trj.xyz

    # normalize bond vectors
    bondlength = (((bondvec**2).sum(2))**0.5)
    bondvec /= bondlength.reshape([bondvec.shape[0], bondvec.shape[1], 1])

    # estimate S2 with ensemble average formula from:
    # Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
    S2 = np.zeros_like(bondvec[0,:,0])
    for i in range(3):
        for j in range(3):
            S2 += (bondvec[:,:,i] * bondvec[:,:,j]).mean(0)**2
    S2 = list(0.5 * (3*S2 - 1))

    info = {}
    info['dt'        ] = trj.timestep
    info['S2'        ] = S2
    info['bondlength'] = bondlength
    info['resnames'  ] = [[], []]
    info['resid'     ] = [[], []]
    info['resindex'  ] = [[], []]
    info['atomnames' ] = [[], []]
    info['atomindex' ] = [[], []]
    info['element'   ] = [[], []]
    info['chain'     ] = [[], []]

    # get info on atoms and residues
    #for atom1, atom2 in zip(atom1_trj.top.atoms, atom2_trj.top.atoms):
    for atom1_ndx, atom2_ndx in bondvec_ndx:
        atom1 = trj.topology.atom(atom1_ndx)
        atom2 = trj.topology.atom(atom2_ndx)
        info['resnames'  ][0].append(atom1.residue.name)
        info['resnames'  ][1].append(atom2.residue.name)
        info['resid'     ][0].append(atom1.residue.resSeq)
        info['resid'     ][1].append(atom2.residue.resSeq)
        info['resindex'  ][0].append(atom1.residue.index)
        info['resindex'  ][1].append(atom2.residue.index)
        info['atomnames' ][0].append(atom1.name)
        info['atomnames' ][1].append(atom2.name)
        info['atomindex' ][0].append(atom1.index) 
        info['atomindex' ][1].append(atom2.index) 
        info['element'   ][0].append(atom1.element) 
        info['element'   ][1].append(atom2.element)  
        info['chain'     ][0].append(atom1.residue.chain.index) 
        info['chain'     ][1].append(atom2.residue.chain.index)   

    return bondvec, info
 

# ============================================================================ #

 
def _angular_correlation_function(bondvec, verbose=True):
    """
    Compute 2nd order Legandre polynomial correlation functions of bond vectors.

    Parameters
    ----------
    bondvec : (nframes, nvectors, ndim) array
        normalized bond vectors
    verbose : boolean, optional
        Print progress messages or not

    Returns
    -------
    corr : (nvectors, nframes/2) array
        correlation functions
    corr_std : (nvectors, nframes/2) array
        standard deviations of correlation values
    corr_stdmean : (nvectors, nframes/2) array
        standard error of the mean
    """

    # compute bond vector correlation functions
    nframes, nvectors, ndim  = bondvec.shape
    corrnframes   = int(round(0.5*nframes))
    corr          = np.zeros([nvectors, corrnframes]) # correlation functions
    corr_std      = np.zeros([nvectors, corrnframes]) # standard deviation
    corr_stdmean  = np.zeros([nvectors, corrnframes]) # standard error of the mean
 
    for vector in range(nvectors):

        if verbose:
            print("\rProgress: {:3.0f}%".format(100.0*(vector+1)/nvectors), end="")
            sys.stdout.flush()

        # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
        P2 = np.polynomial.legendre.legval(np.dot(bondvec[:,vector,:], bondvec[:,vector,:].T), [0,0,1])

        # compute the correlation function for each lag time
        for frame in range(corrnframes):
            d = np.diagonal(P2, frame)
            corr        [vector, frame] = d.mean()
            corr_std    [vector, frame] = d.std()
            corr_stdmean[vector, frame] = corr_std[vector, frame] / d.shape[0]**0.5 

    if verbose:
        print()

    return corr, corr_std, corr_stdmean

 
# ============================================================================ #


def save_corr(savefilename, corr, corrstd, corrstdmean, bondvecinfo, topfilename, trjfilename, frames):
    """
    Save correlation functions to disk including all additional information.

    Parameters
    ----------
    savefilename : string
        Filename of the resulting zipfile containing all the data

    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    bondvecinfo : dict
        Dictionary with information on bondlength, resids and resnames

    trjfilename : string
        Filename of the trajectory the original MD data was taken from.
        Stored as is, so better provide a complete path.

    topfilename : string
        Filename of the topology the original MD data was taken from.
        Stored as is, so better provide a complete path.

    frames : (startframe, endframe) tuple
        Range of frames used for correlation function computation
    """

    tmpdirnamelen     = 8
    tmpdirname        = "tmpdir_" + ''.join(np.random.choice([i for i in string.letters + string.digits], size=tmpdirnamelen))
    os.mkdir(tmpdirname)
    savefilepath      = os.path.dirname(savefilename)
    npzfilename       = "corr.npz"
    infofilename      = "info.dat"

    if not os.path.isdir(savefilepath):
        os.mkdir(savefilepath)

    info = {}
    info['npzfilename'] = npzfilename
    info['bondvecinfo'] = bondvecinfo
    info['trjfilename'] = trjfilename
    info['topfilename'] = topfilename
    info['frames']      = frames
 
    # save arrays
    with open(tmpdirname + '/' + npzfilename, 'w') as outfile:
        #np.savez(outfile, corr, corrstd, corrstdmean)
        np.savez_compressed(outfile, corr=corr, corrstd=corrstd, corrstdmean=corrstdmean)

    # save info
    with open(tmpdirname + '/' + infofilename, 'w') as outfile:
        outfile.write(bz2.compress(pickle.dumps(info)))

    # pack both into zipfile
    with zipfile.ZipFile(savefilename, 'w') as outfile:
        outfile.write(tmpdirname + '/' + infofilename, arcname=infofilename)
        outfile.write(tmpdirname + '/' +  npzfilename, arcname=npzfilename )

    # remove npz and pickle files
    os.remove(tmpdirname + '/' + infofilename)
    os.remove(tmpdirname + '/' +  npzfilename)
    os.rmdir (tmpdirname)


# ============================================================================ #


def load_corr(filename):
    """
    Load correlation functions including all additional information.

    Parameters
    ----------
    filename : string
        Filename of the zipfile containing all the data

    Returns
    -------
    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    info
        Dictionary with information on the loaded correlation functions
    """

    # create and change to working directory
    cwd           = os.getcwd()
    absfilename   = os.path.abspath(filename)
    tmpdirnamelen = 8
    tmpdirname    = "/tmp/" + ''.join(np.random.choice([i for i in string.letters + string.digits], size=tmpdirnamelen))
    os.mkdir(tmpdirname)
    os.chdir(tmpdirname)

    # extract files
    with zipfile.ZipFile(absfilename, 'r') as infile:
        zipfilenames = infile.namelist()
        for zipfilename in zipfilenames:
            infile.extract(zipfilename)

    for zipfilename in zipfilenames:
        with open(zipfilename, 'r') as infile:
            # load info dictionary
            if zipfilename.find('.dat') >= 0:
                info = pickle.loads(bz2.decompress(infile.read()))

            # load corr arrays
            if zipfilename.find('.npz') >= 0:
                arrays      = np.load(infile)
                corr        = arrays['corr']
                corrstd     = arrays['corrstd']
                corrstdmean = arrays['corrstdmean']

        # remove extracted files
        os.remove(zipfilename)

    # change back to original cwd
    os.chdir(cwd)
    os.rmdir(tmpdirname)

    return corr, corrstd, corrstdmean, info


# ============================================================================ #
 

def bondvec_corr(trj, bondvec=None, fitgroup=None, parallel=True, saveinfo=None, verbose=True):
    """
    Compute bond vector correlation functions from mdtraj trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to compute the order parameters of.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.

    parallel : boolean, optional
        Do the fitting in parallel on as many processors as available on the machine.

    saveinfo : dict, optional
        Information needed to save the correlation functions.
        Required keys are:
            'topfilename'   the topology   file the data comes from
            'trjfilename'   the trajectory file the data comes from
            'frames'        tuple with initial and final frame
            'zipfilename'   name of the zipfile to store correlation functions in

    verbose : boolean, optional
        Toggle verbosity level.


    Returns
    -------
    corr : (nvec, t) array
        Bond vector correlation functions.

    bondvecinfo :
        Dictionary with information on the bond vectors.
    """

    # fit trajectory
    if fitgroup is not None:
        _fit_trj(trj, fitgroup=fitgroup, parallel=parallel)

    bondvec_ndx                = _vec_atoms(trj, bondvec)
    bondvectors, bondvecinfo   = _bond_vec(trj, bondvec_ndx)
    corr, corrstd, corrstdmean = _angular_correlation_function(bondvectors, verbose=verbose)

    # store additional bond vector information
    if type(fitgroup) != type(None):
        bondvecinfo['fit' ] = True
    else:
        bondvecinfo['fit' ] = False
    bondvecinfo['fitgroup'] = fitgroup
    bondvecinfo['bondvec' ] = bondvec

    if type(saveinfo) == type(dict()):
        save_corr(saveinfo['zipfilename'], corr, corrstd, corrstdmean, bondvecinfo, saveinfo['topfilename'], saveinfo['trjfilename'], saveinfo['frames'])

    return corr, bondvecinfo

    
# ============================================================================ #


def bondvec_corr_batch_mpi(topfilename, trjfilenames, savepath, subtrjlength=None, bondvec=None, fitgroup=None):
    """
    Compute bond vector correlation functions from mdtraj trajectory.
    MPI parallelized batch version that stores the results on disk.

    Parameters
    ----------
    topfilename : string
        topology filename, e.g. PDB

    trjfilenames : sequence of strings
        trajectory filenames, e.g. XTC

    savepath : string
        path to store the results at

    subtrjlength: float, optional
        Subtrajetory length in ps.
        If none, whole trajectories are used.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.
    """

    # set up some timer and counter
    tc = {}
    tc['runtimer']  = time.time()
    tc['loadtimer'] = 0.0
    tc['corrtimer'] = 0.0
    tc['nsubtrjs' ] = 0
    tc['nframes'  ] = 0
 

    # set up MPI
    comm   = MPI.COMM_WORLD
    nranks = comm.Get_size()
    myrank = comm.Get_rank()
    root   = 0

    # to fit or not to fit?
    fit = "nofit"
    if type(fitgroup) != type(None):
        fit = "fit"

    # distribute tasks over ranks
    if myrank == root:
        ntasks        = len(trjfilenames)
        ntasksperrank = np.array(nranks * [ntasks / nranks], dtype=np.int)
        ntasksperrank[:ntasks%nranks] += 1

        for rank in range(0, nranks):
            taskstaken = sum(ntasksperrank[:rank])
            task = {}
            task['topfilename']  = topfilename
            task['trjindices']   = range(rank, len(trjfilenames), nranks)
            task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
            task['savepath']     = savepath

            comm.send(task, dest=rank, tag=rank)

    task = comm.recv(source=root, tag=myrank)
#    print ("rank {}: ".format(myrank), task['trjindices'], task['trjfilenames'])
#    sys.exit(0)

    # do the assigned piece of work
    for nf, trjfilename in enumerate(task['trjfilenames']):
#        if nf > 3:
#            break
        # determinde dt and chunksize
        trjs      = mdtraj.iterload(trjfilename, top=task['topfilename'], chunk=2)
        trj       = trjs.next()
        dt        = trj.timestep
        chunksize = int(subtrjlength / dt)
        trjindex  = task['trjindices'][nf]

        if myrank == root:
            print("\r", 100*" ", end="")
            print("\rProgress rank root: {:3.0f}% ".format(100.0*nf/len(task['trjfilenames']), nf), end="")
            sys.stdout.flush()
 
        loadstarttime = time.time()
        for ntrj, trj in enumerate(mdtraj.iterload(trjfilename, top=task['topfilename'], chunk=chunksize)):
#            if ntrj > 3:
#                break

            tc['loadtimer'] += time.time() - loadstarttime

            if trj.n_frames == chunksize:
                tc['nsubtrjs' ] += 1
                tc['nframes']   += trj.n_frames

                saveinfo = {}
                saveinfo['topfilename'] = task['topfilename']
                saveinfo['trjfilename'] = trjfilename
                saveinfo['frames'     ] = (ntrj * chunksize, (ntrj+1)*chunksize-1)
                saveinfo['zipfilename'] = savepath + '/' + ''.join(os.path.basename(trjfilename).split('.')[:-1]) + '_{}.zip'.format(fit)
                saveinfo['zipfilename'] = "{}/{}_no{}_frames{}to{}_{}.zip".format(savepath,
                                                                                  ''.join(os.path.basename(trjfilename).split('.')[:-1]),
                                                                                  trjindex,
                                                                                  saveinfo['frames'][0],
                                                                                  saveinfo['frames'][1],
                                                                                  fit)
            
                corrstarttime = time.time()
                bondvec_corr(trj, bondvec=bondvec, fitgroup=fitgroup, parallel=False, saveinfo=saveinfo, verbose=False)
                tc['corrtimer'] += time.time() - corrstarttime
            loadstarttime = time.time()

            if myrank == root:
                print(".", end="")
                sys.stdout.flush()

    if myrank == root:
        print("")


    # report timers and counters to root
    if myrank == root:
        for rank in range(1, nranks):
            tci = comm.recv(source=rank, tag=rank)
            tc['loadtimer'] += tci['loadtimer'] 
            tc['corrtimer'] += tci['corrtimer'] 
            tc['nsubtrjs' ] += tci['nsubtrjs' ] 
            tc['nframes'  ] += tci['nframes'  ] 
    else:
        comm.send(tc, dest=root, tag=myrank)

    comm.barrier()

    tc['runtimer'] = time.time() - tc['runtimer']
    if myrank == root:
        print("Summary")
        print("-------\n")
        print("Loaded {} subtrajectories with a total of {} frames".format(tc['nsubtrjs'], tc['nframes']))
        print("Total runtime:                        {:8.0f} sec.".format(tc['runtimer' ]))
        print("Aggregated time for loading:          {:8.0f} sec.".format(tc['loadtimer']))
        print("Aggregated time for corr computation: {:8.0f} sec.".format(tc['corrtimer']))
        

# ============================================================================ #


def _check_corr_convergence(corr, diffThreshold=0.02, stdThreshold=0.02):
    """
    Check the convergence of the bond vector correlation functions

    Parameters
    ----------
    corr : (nbonds, nframes) array
        Correlation functions

    diffThreshold: float, optional
        Maximum mean difference for convergence check.

    stdThreshold: float, optional
        Maximum stdev difference for convergence check.

    Returns
    -------
    

        (convergence values, boolean array of converged correlation functions)
    """
    length            = corr.shape[1]
    quarter           = length/4
    thirdQuarter      = corr[:,2*quarter:3*quarter]
    fourthQuarter     = corr[:,3*quarter:4*quarter]
    fourthQuarterMean = fourthQuarter.mean(1)
    difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
    stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
    convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)
    return fourthQuarterMean, convergence
 


# ============================================================================ #


def order_parameter_mean(corr, converged=True, diffThreshold=0.02, stdThreshold=0.02):

    length            = corr.shape[1]
    quarter           = length/4
    thirdQuarter      = corr[:,2*quarter:3*quarter]
    fourthQuarter     = corr[:,3*quarter:4*quarter]
    fourthQuarterMean = fourthQuarter.mean(1)

    difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
    stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
    convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)

    return fourthQuarterMean, convergence 


# ============================================================================ #


def order_parameter(corrfilenames, method="mean", converged=True, verbose=True, **kwargs):
    """
    Compute bond vector order parameters from correlation functions.

    Parameters
    ----------
    corrfilenames : list
        List with *.zip filenames containing correlation functions.
        
    method: string, optional
        The method to use for order parameter computation.
        Options are:
            "mean"          Use the mean of the final quarter as order parameter
            "single exp"    Fit correlation functions to single exponential
            "double exp"    Fit correlation functions to double exponential
            "extLS"         Use the extended least squares method (method 3) from:
                            JPCB 2008, 112, 6155-6158, pubs.acs.org/doi/abs/10.1021/jp077018h
            "iRED"          Use the iRED method from:
                            Brüschweiler, JACS 2002, 124, 4522-4534, pubs.acs.org/doi/abs/10.1021/ja012750u

    converged : boolean, optional, default: True
        Use only converged correlation functions for averaging

    verbose : boolean, optional
        Report progress and other verbose output.

    **kwargs : optional keyword arguments
        All remaining keyword arguments are passed to the order parameter method.
    """

    starttime = time.time()

    # load correlation functions
    corrlist = []
    print("Loading {} set{} of correlation functions:".format(len(corrfilenames), *['s' if i>1 else '' for i in [len(corrlist)]]))
    for nf, filename in enumerate(corrfilenames):
        corr, corrstd, corrstdmean, info = load_corr(filename)
        corrlist.append(corrFunction(corr, corrstd, corrstdmean, info))
        if verbose:
            print("\rProgress: {:3.0f}%".format(100.0*nf/len(corrfilenames)), end="")
            sys.stdout.flush()

    # select order parameter estimation method
    if method == "mean":
        S2 = order_parameter_mean(corrlist, **kwargs)
    else:
        print("Order parameter estimation method unknown: {}".format(method))
        sys.exit(1)

    # report runtime
    print("\rRuntime: {:.2f} sec.".format(time.time() - starttime))

    return S2
 

# ============================================================================ #


def n_exp(x, *args):
    """
    Evaluate the sum of n decaying exponentials depending on the length of *args.
    n coefficients and n exponential decay constants are expected.
    The coefficients are constraint to sum up to 1:
        y = sum_i(ai * exp(x/-ti)) + 1 - sum_i(ai)
    Example call:
        n_exp(x, a1, a2, a3, t1, t2, t3)
    """
    coefficients = args[:len(args)/2]
    decayconst   = args[len(args)/2:]

    if (np.array(decayconst) < 0).sum() > 0:
        return x * float('NaN')

    y = 0
    for i, c, t in zip(range(len(decayconst)), coefficients, decayconst):
        y += c * np.exp(x/(-1*t))
    y += 1 - sum(coefficients)

    return y


# ============================================================================ #


def n_exp_obj(p, x, y):
    """
    Objective function for least squares minimization for n exponentials.

    Parameters
    ----------
    p: Parameters for n_exp
    x: n_exp(x)
    y: y = n_exp(x)
    """

    sumsquares = ((y - n_exp(x, *p))**2).sum()
    return sumsquares


# ============================================================================ #


def LipariSzabo(t, S2, te, tm):
    """
    Evaluate the Liapri Szabo model:
    C(t) = exp(-t/tm) * (S2 + (1-S2)*exp(-t/te)))
    """
    C = np.exp(-t/tm) * (S2 + (1-S2)*np.exp(-t/te))
    return C


# ============================================================================ #


def extLipariSzabo(t, S2, Sf, te, tm):
    """
    Evaluate the extended Liapri Szabo model:
    C(t) = exp(-t/tm) * (S2 + (Sf-S2)*exp(-t/te)))
    """
    C = np.exp(-t/tm) * (S2 + (Sf-S2)*np.exp(-t/te))
    return C
 

# ============================================================================ #

def CloreLS(t, S2, Sf, tf, ts, tm):
    """
    CloreLS model from:
    Clore, M.; Szabo, A.; et al.
    Deviations from the Simple Two-Parameter Model-Free Approach to the Interpretation of Nitrogen-15 Nuclear Magnetic Relaxation of Proteins.
    J. Am. Chem. Soc. 1990.
    """

    C_I = S2 + (1-Sf)*np.exp(-t/tf) + (Sf-S2)*np.exp(-t/ts)
    C   = np.exp(-t/tm) * C_I 

    return C

# ============================================================================ #

def CloreLS_obj(p, t, C):
    """
    Objective function for CloreLS least squares fit
    
    Parameters
    ----------
    p : sequence
        Parameters for CloreLS model
    t : (N,) array
        Time sequence
    C : (N,) array
        Correlation function
    """
    sumsquares = ((C - CloreLS(t, *p))**2).sum()
    return sumsquares

# ============================================================================ #
 
def extCloreLS(t, S2, Sf, Si, tf, ts, tm, motionAveraged=False):
    """
    extends the CloreLS model by a immediate component, Si, that is quicker than the fast component:
    Clore, M.; Szabo, A.; et al.
    Deviations from the Simple Two-Parameter Model-Free Approach to the Interpretation of Nitrogen-15 Nuclear Magnetic Relaxation of Proteins.
    J. Am. Chem. Soc. 1990.

    motionAveraged: Only use the motionally averaged model
    """

    C_I = S2 + (Si-Sf)*np.exp(-t/tf) + (Sf-S2)*np.exp(-t/ts)

    if motionAveraged:
        return C_I
    else:
        C   = np.exp(-t/tm) * C_I 
        return C
 
# ============================================================================ #

def extCloreLS_obj(p, t, C, motionAveraged=False):
    """
    Objective function for extCloreLS least squares fit
    
    Parameters
    ----------
    p : sequence
        Parameters for extCloreLS model
    t : (N,) array
        Time sequence
    C : (N,) array
        Correlation function
    """
    sumsquares = ((C - extCloreLS(t, *p))**2).sum()
    return sumsquares

# ============================================================================ #

def generalLS(t, S, tau, Sf=1.0, tm=float('inf')):
    """
    Generalization of the Lipari Szabo model for Bond vector correlation function fits.

    Parameters
    ----------
    t : (N,) array
        timeseries array
    S : (M,) array
        S2 and all other decay coefficients to use in the model.
    tau: (M,) array]
        decay time constants
    Sf : float, default 1.0
        initial decay before first data point
    tm : float, default 'inf'
        global rotation decay time constant.
        Set to float('inf') for internal motion only (default).

    Retruns
    -------
    C : (N,) array
        Correlation function estimates for fitting
    """

#    print('S:  ', S)
#    print('tau:', tau)
#    print('Sf: ', Sf)
#    print('tm: ', tm)
#    print(' ')

    S_ = list(S) + [Sf]

    #C_I = S0 + (S1-S0)*np.exp(-t/t0) + (S2-S1)*np.exp(-t/t1) + ...
    #    = S0 + sum_i [(Si+1 - Si)*np.exp(-t/ti)]

    C_I = S_[0]
    for i, taui in enumerate(tau):
        C_I += (S_[i+1] - S_[i]) * np.exp(-t/taui)

    C = C_I * np.exp(-t/tm)

    return C

# ============================================================================ #
 
#def generalLS_obj(p, t, C, sigma=None, fast=False, internal=False):
#    """
#    Objective function for least squares minimization of 
#    generalization of the Lipari Szabo model for Bond vector correlation function fits.
#
#    Parameters
#    ----------
#    p : (N,) array
#        generalLS parameters: [S0, .., Sn, (Sf), t1, .. tn, (tm)]
#        shape: (N,), with N = 2*n + 1 (+1) (+1)
#        Sf is only present if fast     == True
#        tm is only present if internal == True
#    t : (M,) array
#        time series
#    C : (M,) array
#        True correlation function values
#    sigma : (M,) array, optional
#        Optionally variances of the correlation function for weighting
#    fast : float, optional
#        If True, p contains an Sf parameter without decay constant describing decay before the first datapoint
#        Minimization should then be done without the first datapoint C(t=0.0) = 1.0
#    internal : float, optional
#        If False, p contains a tm decay constant accounting for global rotations.
#
#    Returns
#    -------
#    sumsqurares : float
#        Sum of squares of correlation function deviations from the model with the current parameters.
#        (optinally weighted by variances)
#    """
#
#    if fast and not internal:
#        # len(p) = 2*n + 1 + 2
#        n   = (len(p) - 3) / 2
#        tau = p[n+2:-1]
#        Sf  = p[n+1]
#        tm  = p[-1]
#    elif fast and internal:
#        # len(p) = 2*n + 1 + 1
#        n   = (len(p) - 2) / 2
#        tau = p[n+2:] 
#        Sf  = p[n+1]
#        tm  = float('inf')
#    elif not fast and not internal:
#        n   = (len(p) - 2) / 2
#        tau = p[n+1:-1]  
#        Sf  = 1.0
#        tm  = p[-1]
#    elif not fast and internal:
#        n   = (len(p) - 1) / 2
#        tau = p[n+1:]  
#        Sf  = 1.0
#        tm  = float('inf')
#
#    S = p[:n+1]
#
#    if sigma is None:
#        sigma = np.ones_like(C)
#
#    print('S:  ', S)
#    print('tau:', tau)
#    print('Sf: ', Sf)
#    print('tm: ', tm)
#    print('')
#
#    sumsquares = ((generalLS(t, S, tau, Sf=Sf, tm=tm) - C)**2 / sigma).sum()
#    return sumsquares
#
## ============================================================================ #
#
#def generalLS_fit(t, C, p, n, sigma=None, fast=False, internal=False):
#    """
#    Least squares fit of the
#    generalization of the Lipari Szabo model for Bond vector correlation function fits.
#    """
#
##    # initial parameter values
##    p = list(np.linspace(0.1,0.9,n+1))
##    if fast:
##        p += [0.95]
##    p += n*[1.0]
##    if not internal:
##        p += [1.0]
#
#    # bounds
#    bounds = (n+1)*[(0,1)]
#    if fast:
#        bounds += [(0,1)]
#    bounds += n*[(0,None)]
#    if not internal:
#        bounds += [(0,None)]
#
#    # constraints
#    constraints = []
#    for i in range(n):
#        constraints += [{"type": "ineq", "fun":  lambda x: x[i+1] - x[i]}]
#    if fast:
#        constraints += [{"type": "ineq", "fun":  lambda x: x[n+1] - x[n]}]
#    
##    print('p:', p)
##    print('bounds:', bounds)
##    print('const.:', constraints)
##    print('')
##    generalLS_obj(p, t, C, sigma=sigma, fast=fast, internal=internal)
#    res = minimize(generalLS_obj, x0=p, args=(t, C, sigma, fast, internal), bounds=bounds, constraints=constraints)
#
#    S = res.x[:n+1]
#    Sf = 1.0
#    tm = float('inf')
#    if fast:
#        Sf = res.x[n+1]
#        if not internal:
#            tau = res.x[n+2:-1]
#            tm  = res.x[-1]
#        else:
#            tau = res.x[n+2:]
#    else:
#        if not internal:
#            tau = res.x[n+1:-1]
#            tm  = res.x[-1]
#        else:
#            tau = res.x[n+1:]
#    return (S, tau, Sf, tm)


# ============================================================================ #

def generalLS_I_obj(p, t, C, sigma=None):
    """
    Objective function for general Lipari Szabo model least squares fit.
    Internal, no Sf

    Parameters
    ----------
    p : (N,) array
        Lipari Szabo parameters
    t : (M,) array
        Timeseries
    C : (M,) array
        Correlation function values
    sigma : (M,) array
        variances for weighted least squares fit
    """

    n   = len(p) / 2
    S   = p[:n]
    tau = p[n:]

    variances = sigma
    if sigma is None:
        variances = np.ones_like(C)

    sumsquares  = ((generalLS(t, S, tau, Sf=1.0, tm=float('inf')) - C)**2 / variances).sum()
    return sumsquares

# ============================================================================ #

def generalLS_I_Sf_obj(p, t, C, sigma=None):
    """
    Objective function for general Lipari Szabo model least squares fit.
    Internal, Sf

    Parameters
    ----------
    p : (N,) array
        Lipari Szabo parameters
    t : (M,) array
        Timeseries
    C : (M,) array
        Correlation function values
    sigma : (M,) array
        variances for weighted least squares fit
    """

    n   = len(p) / 2
    S   = p[:n]
    tau = p[n:2*n]
    Sf  = p[-1]

    variances = sigma
    if sigma is None:
        variances = np.ones_like(C)

    sumsquares  = ((generalLS(t, S, tau, Sf=Sf, tm=float('inf')) - C)**2 / variances).sum()
    return sumsquares
 
# ============================================================================ #

def generalLS_Sf_obj(p, t, C, sigma=None):
    """
    Objective function for general Lipari Szabo model least squares fit.
    global, Sf

    Parameters
    ----------
    p : (N,) array
        Lipari Szabo parameters
    t : (M,) array
        Timeseries
    C : (M,) array
        Correlation function values
    sigma : (M,) array
        variances for weighted least squares fit
    """

    print(p)

    n   = (len(p)-2) / 2
    S   = p[:n]
    tau = p[n:2*n]
    Sf  = p[-1]
    tm  = p[-2]

    variances = sigma
    if sigma is None:
        variances = np.ones_like(C)

    sumsquares  = ((generalLS(t, S, tau, Sf=Sf, tm=tm) - C)**2 / variances).sum()
    return sumsquares
 
# ============================================================================ #

def generalLS_obj(p, t, C, sigma=None):
    """
    Objective function for general Lipari Szabo model least squares fit.
    global, no Sf

    Parameters
    ----------
    p : (N,) array
        Lipari Szabo parameters
    t : (M,) array
        Timeseries
    C : (M,) array
        Correlation function values
    sigma : (M,) array
        variances for weighted least squares fit
    """

    n   = len(p) / 2
    S   = p[:n]
    tau = p[n:2*n]
    tm  = p[-1]

    variances = sigma
    if sigma is None:
        variances = np.ones_like(C)

    sumsquares  = ((generalLS(t, S, tau, Sf=1.0, tm=tm) - C)**2 / variances).sum()
    return sumsquares
 

# ============================================================================ #

def generalLS_fit(t, C, n, sigma=None, fast=False, internal=False):
    """
    Fit the general Lipari Szabo model with a sum of n exponentials.

    Parameters
    ----------
    t : (N, ) array
        timeseries
    C : (N, ) array
        correlation function values
    n : int
        Number of exponentials for fit
    sigma : (N, ) array, optional
        Variances for least squares fit weighting
    fast : bool, optional
        Include a parameter for fast decay before the first datapoint.
        Omit first correlation function datapoint equalint 1.0
    internal : bool, optional
        Only consider internal bond vector motions

    Returns
    -------
    p : dict
        Dictionary containing all optimized parameters for generalLS()
    """

    # initialize parameter guesses
    p0  = list(np.linspace(0.1, 0.9, n))     # S
    p0 += list(np.logspace(1, n, n, base=10)) # tau
    if fast:
        p0 += [0.95]                         # Sf
    if not internal:                         # tm
        p0 += [10.0**n]

    # set bounds
    bounds  = n * [(0, 1)]    # S
    bounds += n * [(0, None)] # tau
    if fast:
        bounds += [(0, 1)]    # Sf
    if not internal:
        bounds += [(0, None)] # tm
    bounds = tuple(bounds)

    # set constraints
    constraints = []
    for i in range(n-1):
        constraints += [{"type": "ineq", "fun":  lambda x: x[i+1] - x[i]}]    # Si+1 > Si
        print (i+1, ">", i)
    if fast:
        constraints += [{"type": "ineq", "fun":  lambda x: x[2*n] - x[n-1]}]  # Sf > Sn-1
        print (2*n, ">", n-1)
    for i in range(n,2*n-1):
        constraints += [{"type": "ineq", "fun":  lambda x: x[i+1] - x[i]}]    # tau_{i+1} > tau_{i}
        print (i+1, ">", i)
    if not internal:
        constraints += [{"type": "ineq", "fun":  lambda x: x[-1] - x[2*n-1]}] # tau_{m} > tau_{2*n-1}
        print (-1, ">", 2*n-1)

    print(p0)
    print(bounds)

    # minimize objective function
    if not fast and internal:
        res = minimize(generalLS_I_obj,    x0=p0, args=(t, C, sigma), bounds=bounds, constraints=constraints)
    elif not fast and not internal:
        res = minimize(generalLS_obj,      x0=p0, args=(t, C, sigma), bounds=bounds, constraints=constraints)
    elif fast and internal:
        res = minimize(generalLS_I_Sf_obj, x0=p0, args=(t, C, sigma), bounds=bounds, constraints=constraints)
    elif fast and not internal:
        res = minimize(generalLS_Sf_obj,   x0=p0, args=(t, C, sigma), bounds=bounds, constraints=constraints)

    p =  {'S': res.x[:n], 'tau': res.x[n:2*n], 'Sf': 1.0, 'tm': float('inf')}
    if fast:
        p['Sf'] = res.x[-2]
    if not internal:
        p['tm'] = res.x[-1]
    return p

# ============================================================================ #

