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
import LS

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

    def __init__(self, corrfilenames=None, converged=True, verbose=True, sort=True, **kwargs):
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
            self.load(sort=sort)
 
# ==================================== #

    def load(self, corrfilenames=None, sort=True):
        """
        Load sets of correlation functions.

        Parameters
        ----------
        corrfilenames : list
            List with *.zip filenames containing correlation functions.
        sort: boolean
            if True, sort corrfilenames before loading
 
        """

        if corrfilenames is not None:
            self.corrfilenames = corrfilenames

        # The order in which correlationfunctions are loaded can affect the computation
        # of the mean and standard deviations later on.
        # To ensure reproducability, sort them here!
        if sort:
            self.corrfilenames.sort()

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

    def average_corr(self, offset=0):
        """
        Compute the average correlation function and its standard deviation

        Parameters:
        -----------
        offset: int
            determines the order in which the correlation functions are averaged
        """
        nvec, nframes = self.corrlist[0].corr.shape
        ncorr         = len(self.corrlist)
        allcorr = np.zeros([ncorr, nvec, nframes], dtype=self.corrlist[0].corr.dtype)

        for i, c in enumerate(range(offset, ncorr+offset)):
            allcorr[i,:,:] = self.corrlist[c % ncorr].corr
        #for c, corr in enumerate(self.corrlist):
        #    allcorr[c,:,:] = corr.corr

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

        return allcorr

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
        elif self.method == "generalLS":
            self.estimate_generalLS(**kwargs)
        elif self.method == "generalLSselection":
            self.estimate_generalLS_modelSelection(**kwargs)
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
                self.S2mean[n]  = self.S2all[n,self.S2convergence[n,:]].mean()
                self.S2std[n]   = self.S2all[n,self.S2convergence[n,:]].std()
                self.S2error[n] = self.S2std[n] / self.S2convergence[n,:].sum()**0.5

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

    def estimate_generalLS(self, n, fast=True, internal=False, weighted=False, **kwargs):

        dt      = self.avgcorr.dt
        ncorr   = self.avgcorr.corr.shape[0]
        nframes = self.avgcorr.corr.shape[1]
        t       = np.linspace(0, dt*nframes, nframes)
        firstf  = 0
        if fast:
            firstf = 1

        self.para = []
        self.S2   = np.zeros([ncorr, n])
        self.tau  = np.zeros([ncorr, n])

#        nS   = n
#        ntau = n 
#        startvals  = []
#        startvals.append(2*n*[1.0])
#        startvals.append(list(np.linspace(0.1, 0.9, nS)) + ntau*[1.0])
#        startvals.append(nS*[1.0] + list(np.logspace(0, ntau-1, ntau))[-1::-1])
#        startvals.append(list(np.linspace(0.1, 0.9, nS)) + list(np.logspace(0, ntau-1, ntau))[-1::-1])
#        self.lsq = np.zeros([ncorr, len(startvals)])

        for nc in range(ncorr):
            if weighted:
                self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:], sigma=self.avgcorr.std[nc,firstf:])
            else:
                self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:])
            p = self.ls.fit(n, fast=fast, internal=internal, **kwargs)

#            for i, p0 in enumerate(startvals):
#                p = self.ls.fit(n, fast=fast, internal=internal, p0=p0, **kwargs)
#                self.lsq[nc,i] = p["lsq"]
#
#            print("{:3d} done".format(nc))

            self.para.append(p)
            self.S2[nc,:]  = p["S"]
            self.tau[nc,:] = p["tau"]

# ==================================== #

    def check_overfitting(self, parameters, mintauratio=2.0, minS2diff=0.01):
        """
        Check generalLS parameters for overfitting
        """
        overfitted = False

        if not parameters['success']:
            return True

        if len(parameters['tau']) > 1:
            minratio = min(parameters['tau'][:-1] / parameters['tau'][1:])
            if minratio < mintauratio:
                overfitted = True
#                print("tau overfitted:", parameters['tau'])

        if len(parameters['S']) > 1:
            mindiff = min(abs(parameters['S'][:-1] - parameters['S'][1:]))
            if mindiff < minS2diff:
                overfitted = True
#                print("S overfitted:", parameters['S'])

        if np.isnan(parameters['p']).sum() > 0:
            overfitted = True

        return overfitted


# ==================================== #

    def estimate_generalLS_modelSelection(self, fast=True, internal=False, weighted=False, maxdecays=int(1e3), nfits=1, **kwargs):
        """
        nfits
        The number of fits to perform per amino acid.
        The order of the correlation function for the mean correlation function is randomized for each fit.
        As fitting is an ill-posed problem, we get better estimates in this way
        """

        randomstate = np.random.get_state()
        np.random.seed(23)

        dt      = self.avgcorr.dt
        ncorr   = self.avgcorr.corr.shape[0]
#        ncorr   = 10
        nframes = self.avgcorr.corr.shape[1]
        t       = np.linspace(0, dt*nframes, nframes)
        firstf  = 0
        if fast:
            firstf = 1

        self._AIC      = np.zeros([ncorr, nfits, maxdecays])
        self._lsq      = np.zeros([ncorr, nfits, maxdecays])
        self._para     = np.zeros([ncorr, nfits, maxdecays, 2*maxdecays])
        self._S2       = np.zeros([ncorr, nfits, maxdecays,   maxdecays])
        self._tau      = np.zeros([ncorr, nfits, maxdecays,   maxdecays])
        self._success  = np.zeros([ncorr, nfits, maxdecays], dtype=np.bool)
        self._paralist = [[[] for nfit in range(nfits)] for nc in range(ncorr)]

        self.para = []
#        totaldecays = 1 # maximum number of decays with successful fit

        def make_progress_msg(progress_percent, ETA):
            progress_msg = "Progress: {:3.0f}% ".format(progress_percent)
            if ETA > 60*60*24:
                ETA_msg = "(ETA: {:5.1f} days)".format(ETA/(60*60*24))
            elif ETA > 60*60:
                ETA_msg = "(ETA: {:5.1f} hours)".format(ETA/(60*60))
            elif ETA > 60:
                ETA_msg = "(ETA: {:5.1f} minutes)".format(ETA/(60))
            else:
                ETA_msg = "(ETA: {:5.0f} seconds)".format(ETA)
            return progress_msg + ETA_msg


        ETA = 0.0
        progress_msg = make_progress_msg(0.0, ETA)
        print(progress_msg, end="")
        sys.stdout.flush()

        starttime = time.time()
        for nc in range(ncorr):

            for nfit in range(nfits):
                print(len(progress_msg)*'\b' + len(progress_msg)*' ' + len(progress_msg)*'\b', end="")
                progress_percent = 100.0*(nc*nfits+nfit)/(nfits*ncorr)
                runtime          = time.time() - starttime
                if progress_percent > 0:
                    ETA          = runtime * (100-progress_percent) / progress_percent
                progress_msg = make_progress_msg(progress_percent, ETA)
                print(progress_msg, end="")
                sys.stdout.flush()

                # compute new average correlation function
                np.random.shuffle(self.corrlist)
                self.average_corr()

                # set up Lipari Szabo fitter
                if weighted:
                    self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:], sigma=self.avgcorr.std[nc,firstf:])
                else:
                    self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:]) 

                # fit for all correlation functions, correlation function shuffles and number of decays
                for ndecays in range(1,maxdecays+1):
                    decay_ndx = ndecays - 1
                    p = self.ls.fit(ndecays, fast=fast, internal=internal, **kwargs)
                    self._paralist[nc][nfit].append(p)
                    self._AIC    [nc, nfit, decay_ndx]             = p['AIC']
                    self._lsq    [nc, nfit, decay_ndx]             = p['lsq']
                    self._para   [nc, nfit, decay_ndx, :2*ndecays] = p['p']
                    self._S2     [nc, nfit, decay_ndx,   :ndecays] = p['S']
                    self._tau    [nc, nfit, decay_ndx,   :ndecays] = p['tau']
                    self._success[nc, nfit, decay_ndx]             = p['success']

        print(len(progress_msg)*'\b' + len(progress_msg)*' ' + len(progress_msg)*'\b', end="")
        progress_msg = make_progress_msg(100.0, 0.0)
        print(progress_msg)

        # reset random number generator state
        np.random.set_state(randomstate)

        # compute probabilities for each model to be the best
        self._probability  = np.exp((self._AIC.min(2, keepdims=True) - self._AIC)/2)
        self._probability /= self._probability.sum(2, keepdims=True)
        self._meanprob     = self._probability.mean(1)

        # select best model based on AIC
        self._bestmodel = self._meanprob.argmax(1)
        self._bestprob  = self._meanprob.max(1)
        self._goodfits  = self._probability.argmax(2) - self._bestmodel.reshape([self._bestmodel.shape[0],1]) == 0

        # store S2 and tau of the best model for each residue
        self.S2      = np.array([self._S2 [nc,self._goodfits[nc,:],self._bestmodel[nc],:].mean(0) for nc in range(ncorr)])
        self.S2std   = np.array([self._S2 [nc,self._goodfits[nc,:],self._bestmodel[nc],:].std (0) for nc in range(ncorr)])
        self.tau     = np.array([self._tau[nc,self._goodfits[nc,:],self._bestmodel[nc],:].mean(0) for nc in range(ncorr)])
        self.taustd  = np.array([self._tau[nc,self._goodfits[nc,:],self._bestmodel[nc],:].std (0) for nc in range(ncorr)]) 
        self.ndecays = self._bestmodel + 1
        self.prob    = self._bestprob

        # store list of one original fitting result for each residue
        self.para = []
        for nc in range(ncorr):
            clearestfit = self._probability[nc,:,self._bestmodel[nc]].argmax()
            self.para.append(self._paralist[nc][clearestfit][self._bestmodel[nc]])
            self.para[-1]['S'  ] = self.S2 [nc,:self.ndecays[nc]]
            self.para[-1]['tau'] = self.tau[nc,:self.ndecays[nc]]
            self.para[-1]['p'  ] = np.concatenate((self.para[-1]['S'], self.para[-1]['tau']))


# ==================================== #

    def plot_all(self, savepath, internal=False):

        fig = plt.figure()
        axs = fig.add_subplot(111)

        nc, npoints  = self.avgcorr.corr.shape
        dt      = self.avgcorr.dt
        t       = np.linspace(0, dt*npoints, npoints)

        if not os.path.isdir(savepath):
            os.mkdir(savepath)

        for c in range(nc):
            corr    = self.avgcorr.corr[c,:]
            corrstd = self.avgcorr.std[c,:]
            correrr = self.avgcorr.error[c,:]
            para    = self.para[c]
            ls      = LS.LS(t, corr)
            ls.internal = internal
            corrfit = ls.generalLS(para['p'])

            coeffs  = ",  ".join([r"$S_{{{}}}^{{2}}={:.2f}$".format(i,S) for i,S in enumerate(para['S'])])
            coeffs += ",  "
            coeffs += ",  ".join([r"$\tau_{{{}}}={:.1f}$".format(i,tau) for i,tau in enumerate(para['tau'])]) 

            lines = []
            lines += axs.plot(t, corr   , 'b', label="Mean correlation function over {} samples".format(len(self.corrlist)), lw=2)
            lines += axs.plot(t, corrfit, 'g', label="fit: {}".format(coeffs), lw=2)
            #lines += axs.fill_between(t, corr+1.96*correrr, corr-1.96*correrr, alpha=0.4, color='r', label="95% Confidence Interval")

            resid   = self.avgcorr.resid[0][c]
            resname = self.avgcorr.resname[0][c]
            atom1   = self.avgcorr.atomname[0][c]
            atom2   = self.avgcorr.atomname[1][c]

            axs.set_ylim(min(0, corr.min()), max(1, corr.max()))
            axs.set_title("{} {} ({}-{})".format(resname, resid, atom1, atom2))
            axs.set_ylabel("correlation")
            axs.set_xlabel("time [ps]")
            lgd = axs.legend(loc="upper right", bbox_to_anchor=(1.2,1)) 

            filename = savepath + "/corr_{:03d}.jpg".format(resid)
            print(filename)
            plt.savefig(filename, dpi=100, papertype='a4', orientation='landscape', bbox_inches='tight', bbox_extra_artists=(lgd,))

            while len(lines) > 0:
                line = lines.pop()
                line.remove()  



        plt.close()

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
                #coeffs = ",  ".join(["{:.1f}".format(n) for n in self.poptavg[self.corridx,:]])
                coeffs  = ",  ".join([r"$S_{{{}}}^{{2}}={:.2f}$".format(i,S) for i,S in enumerate(self.para[self.corridx]['S'])])
                coeffs += ",  "
                coeffs += ",  ".join([r"$\tau_{{{}}}={:.1f}$".format(i,S) for i,S in enumerate(self.para[self.corridx]['tau'])])
            except AttributeError as e:
                print("Did you forget to fit the data?")
                raise e
            xdata  = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, 10*corrFun.corr.shape[1])  
            if self.plotfit == "generalLS":
                self.lines += self.axs.plot(self.ls.t, self.ls.generalLS(self.para[self.corridx]['p']), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit in ["single exp", "double exp"]:
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
        xmin, xmax = self.axs.get_xlim()
        self.lines += self.axs.plot([xmin, xmax], 2*[0.98], 'k', linewidth=2)
        self.lines += self.axs.plot(xmin + self.corridx*(xmax-xmin)/(corrFun.corr.shape[0]-1), 0.98, 'sk', markersize=15)

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

    if not os.path.isdir(savefilepath) and not savefilepath == '':
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

        for rank in range(1, nranks):
            taskstaken = sum(ntasksperrank[:rank])
            task = {}
            task['topfilename']  = topfilename
            task['trjindices']   = range(rank, len(trjfilenames), nranks)
            task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
            task['savepath']     = savepath

#            print("Sending task to rank ", rank)
            comm.send(task, dest=rank, tag=rank)

#        print("Done with sending, receiving")

    if myrank != root:
        task = comm.recv(source=root, tag=myrank)
    else:
        rank = myrank
        taskstaken = sum(ntasksperrank[:rank])
        task = {}
        task['topfilename']  = topfilename
        task['trjindices']   = range(rank, len(trjfilenames), nranks)
        task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
        task['savepath']     = savepath 

#    print ("rank {}: ".format(myrank), task['trjindices'], task['trjfilenames'])
#    sys.stdout.flush()
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










