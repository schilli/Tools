#!/usr/bin/env python

from __future__ import print_function, division

import sys, os, time, itertools, argparse
import numpy as np
import cPickle as pickle
import mdtraj as md
import scipy.constants as cn
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from saltbridges import *

usage = "usage: {} <pdbfile> <gmxtopfile> <trjfile> [<moretrjfiles>]".format(sys.argv[0])

positive_resnames = ["ARG", "LYS"]
negative_resnames = ["ASP", "GLU"]

# Boltzmann's electrostatic constant
ke    = 1/(4*cn.pi*cn.epsilon_0)                  #  J    m / C^2
ke_gu = ke / cn.kilo / cn.nano * cn.e**2 * cn.N_A # kJ * nm / e**2 / mol (gromacs units)

closed_saltbridge_threshold = -320 # kJ/mol


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Compute saltbridge energies and derived properties')
    parser.add_argument('-pdb', action='store', type=str, required=True, dest='pdbfilename', metavar='pdbfilename', help="PDB filename (GRO also allowed)")
    parser.add_argument('-top', action='store', type=str, required=True, dest='gmxtopfilename', metavar='gmxtopfilename', help="Gromacs topology filename")
    parser.add_argument('-trj', action='store', nargs='+', type=str, required=True, dest='trjfilenames', metavar='trjfilenames', help="trajectory filenames")
    parser.add_argument('-sfx', action='store', type=str, required=False, default="", dest='suffix', metavar='suffix', help="Output suffix")
    parser.add_argument('-stride', action='store', type=int, required=False, default=1, dest='stride', metavar='stride', help="Stride for trajectory reader")
    parser.add_argument('-rb', action='store', type=int, required=False, default=None, dest='resSeq_border', metavar='resSeq_border',
            help="Indicate separate protein domains, separated at this residue sequence number")
    args = parser.parse_args()

    main(args.pdbfilename, args.gmxtopfilename, args.trjfilenames, args.stride, args.suffix, args.resSeq_border)

def main(topfilename, gmxtopfilename, trjfilenames, stride=1, suffix="", resSeq_border=None):

    if gmxtopfilename.split('.')[-1] != "top":
        print("Not a gromacs topology file: {}".format(gmxtopfilename))
        print(usage)
        sys.exit(1)

    if len(suffix) > 0 and suffix[0] != "_":
        suffix = "_" + suffix

    Epot = None

    for trjfilename in trjfilenames:
        st = time.time()
        trj = md.load(trjfilename, top=topfilename, stride=stride)
        #trj = md.iterload(trjfilename, top=topfilename, chunk=10).next()
        loadtime = time.time() - st
        filesize = os.stat(trjfilename).st_size/1024**2 # in Mb
        print("Loaded {} frames from {} in {:.1f} sec. (IO rate: {:.2f} Mb/s)".format(trj.n_frames, trjfilename, loadtime, filesize/loadtime))

        protein_residues = [res for res in trj.top.residues if res.is_protein]

        positive_residues = []
        negative_residues = []

#        # all residues
#        for res in protein_residues:
#            positive_residues.append(Residue(res, gmxtop=gmxtopfilename))
#            negative_residues.append(Residue(res, gmxtop=gmxtopfilename))  

#        # positive and negative only
#        positive_residues.append(Residue(list(protein_residues)[ 0], gmxtop=gmxtopfilename)) # N-terminus
#        negative_residues.append(Residue(list(protein_residues)[ 0], gmxtop=gmxtopfilename)) # N-terminus
#        for res in protein_residues:
#            if res.name in positive_resnames or res.name in negative_resnames:
#                positive_residues.append(Residue(res, gmxtop=gmxtopfilename))
#                negative_residues.append(Residue(res, gmxtop=gmxtopfilename)) 
#        positive_residues.append(Residue(list(protein_residues)[-1], gmxtop=gmxtopfilename)) # C-terminus
#        negative_residues.append(Residue(list(protein_residues)[-1], gmxtop=gmxtopfilename)) # C-terminus
 

        # positive and negative separately
        positive_residues.append(Residue(list(protein_residues)[ 0], gmxtop=gmxtopfilename)) # N-terminus
        for res in protein_residues:
            if res.name in positive_resnames:
                positive_residues.append(Residue(res, gmxtop=gmxtopfilename))
            if res.name in negative_resnames:
                negative_residues.append(Residue(res, gmxtop=gmxtopfilename)) 
        negative_residues.append(Residue(list(protein_residues)[-1], gmxtop=gmxtopfilename)) # C-terminus


        #distance_min_, distance_mean_, distance_sc_min_, distance_sc_mean_ = compute_residue_distance(trj, positive_residues, negative_residues, verbose=True)
        Epot_, ideal_dist_ = compute_Epot_saltbridges(trj, positive_residues, negative_residues, verbose=True)

        if Epot is None:
            Epot             = Epot_
            ideal_dist       = ideal_dist_
            #distance_min     = distance_min_
            #distance_mean    = distance_mean_
            #distance_sc_min  = distance_sc_min_
            #distance_sc_mean = distance_sc_mean_
        else:
            Epot       = np.concatenate((Epot, Epot_), axis=2)
            ideal_dist = np.concatenate((ideal_dist, ideal_dist_), axis=2)
            #distance_min     = np.concatenate((distance_min    , distance_min_    ), axis=2)
            #distance_mean    = np.concatenate((distance_mean   , distance_mean_   ), axis=2)
            #distance_sc_min  = np.concatenate((distance_sc_min , distance_sc_min_ ), axis=2)
            #distance_sc_mean = np.concatenate((distance_sc_mean, distance_sc_mean_), axis=2) 

    save_data(Epot, ideal_dist, positive_residues, negative_residues, "data/saltbridges{}.dat".format(suffix))
    plot_Epot_saltbridges(Epot, ideal_dist, positive_residues, negative_residues, "plots", suffix=suffix, resSeq_border=resSeq_border)
    plot_Epot_saltbridges(Epot[:,:,:1], ideal_dist[:,:,:1], positive_residues, negative_residues, "plots", suffix=suffix+"_crystal", resSeq_border=resSeq_border)
    plot_saltbridge_correlations(Epot, positive_residues, negative_residues, "plots", suffix=suffix, resSeq_border=resSeq_border)

    return Epot, ideal_dist, positive_residues, negative_residues #, distance_min, distance_mean, distance_sc_min, distance_sc_mean
     



def save_data(Epot, ideal_dist, residues1, residues2, datafilename):
    """
    Pickle computed saltbridge quantities to disk.
    """
    datadir = os.path.dirname(datafilename)
    if not os.path.isdir(datadir):
        os.mkdir(datadir) 

    with open(datafilename, 'w') as datafile:
        pickler = pickle.Pickler(datafile, protocol=pickle.HIGHEST_PROTOCOL)
        pickler.dump(Epot)
        pickler.dump(ideal_dist)
        pickler.dump(residues1)
        pickler.dump(residues2)



def load_data(datafilename):
    """
    Unpickle computed saltbridge quantities from disk.
    """ 
    with open(datafilename, 'r') as datafile:
        pickler = pickle.Unpickler(datafile)
        Epot       = pickler.load()
        ideal_dist = pickler.load()
        residues1  = pickler.load()
        residues2  = pickler.load()
    return Epot, ideal_dist, residues1, residues2
 



def plot_Epot_saltbridges(Epot, ideal_dist, positive_residues, negative_residues, plotdir, suffix="", resSeq_border=None):
    """
    Plot saltbridge energies in 3 different ways:
        * mean Epot per residue pair
        * minimum Epot per residue pair
        * fraction of time the saltbridge is closed (Epot < closed_saltbridge_threshold) per residue pair

    resSeq_border can be provided as an integer to indicate separate domains separated by this residue sequence number in the plots
    """

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    if len(suffix) > 0 and suffix[0] != '_':
        suffix = "_" + suffix

    maxlabels = 40

    # set labels and ticks
    if len(negative_residues) < maxlabels:
        xlabels = ["{}-{}".format(r.name, r.resSeq) for r in negative_residues]
        xlabels[-1] = "C-Term"
        xticks = range(Epot.shape[1])
    else:
        xtickres = np.round(np.linspace(0, len(negative_residues)-1, 20)).astype(int)
        xlabels = []
        for r_idx in xtickres:
            xlabels.append("{}-{}".format(negative_residues[r_idx].name, negative_residues[r_idx].resSeq))
        xticks = xtickres
    if len(positive_residues) < maxlabels:
        ylabels = ["{}-{}".format(r.name, r.resSeq) for r in positive_residues]
        ylabels[0] = "N-Term"
        yticks = range(Epot.shape[0])
    else:
        ytickres = np.round(np.linspace(0, len(positive_residues)-1, 20)).astype(int)
        ylabels = []
        for r_idx in ytickres:
            ylabels.append("{}-{}".format(positive_residues[r_idx].name, positive_residues[r_idx].resSeq))
        yticks = ytickres

    if resSeq_border is not None:
        for r, residue in enumerate(positive_residues):
            if residue.resSeq >= resSeq_border:
                borderx = r - 0.5
                break
        for r, residue in enumerate(negative_residues):
            if residue.resSeq >= resSeq_border:
                bordery = r - 0.5
                break

    # Epot_mean colormap and vmin/vmax
    Epot_mean = Epot.mean(2)
    vmin, vmax = (min(0, Epot_mean.min()), max(0, Epot_mean.max()))
    if 0 in (vmin, vmax):
        cmap = plt.get_cmap('Blues_r')
        #cmap = plt.get_cmap('Blues')
        #cmap = plt.get_cmap('terrain')
    else:
        cmap = plt.get_cmap('bwr')
        vmin = -max(np.abs(vmin), np.abs(vmax))
        vmax = -vmin

    # plot Epot_mean
    fig = plt.figure()
    axs = fig.add_subplot(111) 
    image = axs.imshow(Epot_mean, interpolation='none', cmap=cmap, vmin=vmin, vmax=vmax)
    axs.set_xticks(xticks)
    axs.set_xticklabels(xlabels, rotation="vertical")   
    axs.set_yticks(yticks)
    axs.set_yticklabels(ylabels, rotation="horizontal") 
    if resSeq_border is not None:
        axs.plot(axs.get_xlim(), 2*[borderx], '-k')
        axs.plot(2*[bordery], axs.get_ylim(), '-k')
    fig.colorbar(image, label='kJ/mol')
    plt.savefig(os.path.join(plotdir, "saltbridges_energy_mean{}.png".format(suffix)), dpi=600, orientation="landscape", papertype="a5", bbox_inches="tight")
    plt.close()
 

    # Epot_min colormap and vmin/vmax
    Epot_min  = Epot.min(2)
    vmin, vmax = (min(0, Epot_min.min()), max(0, Epot_min.max()))
    if 0 in (vmin, vmax):
        cmap = plt.get_cmap('Blues_r')
        #cmap = plt.get_cmap('terrain')
    else:
        cmap = plt.get_cmap('bwr')
        vmin = -max(np.abs(vmin), np.abs(vmax))
        vmax = -vmin

    # plot Epot_min
    fig = plt.figure()
    axs = fig.add_subplot(111) 
    image = axs.imshow(Epot_min, interpolation='none', cmap=cmap, vmin=vmin, vmax=vmax)
    axs.set_xticks(xticks)
    axs.set_xticklabels(xlabels, rotation="vertical")   
    axs.set_yticks(yticks)
    axs.set_yticklabels(ylabels, rotation="horizontal") 
    if resSeq_border is not None:
        axs.plot(axs.get_xlim(), 2*[borderx], '-k')
        axs.plot(2*[bordery], axs.get_ylim(), '-k')  
    fig.colorbar(image, label='kJ/mol')
    plt.savefig(os.path.join(plotdir, "saltbridges_energy_min{}.png".format(suffix)), dpi=600, orientation="landscape", papertype="a5", bbox_inches="tight")
    plt.close()


    # fraction_closed colormap and vmin/vmax
    fraction_closed = 1.0 * (Epot < closed_saltbridge_threshold).sum(2) / Epot.shape[2]
    vmin, vmax = (0, 1)
    cmap = plt.get_cmap('Blues')

    # plot fraction_closed
    fig = plt.figure()
    axs = fig.add_subplot(111) 
    image = axs.imshow(fraction_closed, interpolation='none', cmap=cmap, vmin=vmin, vmax=vmax)
    axs.set_xticks(xticks)
    axs.set_xticklabels(xlabels, rotation="vertical")   
    axs.set_yticks(yticks)
    axs.set_yticklabels(ylabels, rotation="horizontal") 
    if resSeq_border is not None:
        axs.plot(axs.get_xlim(), 2*[borderx], '-k')
        axs.plot(2*[bordery], axs.get_ylim(), '-k') 
    fig.colorbar(image, label='fraction of time closed')
    plt.savefig(os.path.join(plotdir, "saltbridges_fraction_closed{}.png".format(suffix)), dpi=600, orientation="landscape", papertype="a5", bbox_inches="tight")
    plt.close()
 
 


def plot_saltbridge_correlations(Epot, residues1, residues2, plotdir, suffix="", resSeq_border=None):
    """
    Plot all correlations between saltbridge energies

    resSeq_border can be provided as an integer to only plot inter-domain saltbridges of domains
    separated by this residue sequence number
    """

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    if len(suffix) > 0 and suffix[0] != '_':
        suffix = "_" + suffix
 
    fraction_closed = 1.0 * (Epot < closed_saltbridge_threshold).sum(2) / Epot.shape[2]
    res1_indices, res2_indices = np.where(fraction_closed > 0.05)

    # Hacky code for GABARAP: filter only saltbridges between N-terminus (resSeq < 28) and protein body (resSeq >= 28)
    if resSeq_border is not None:
        res1_indices_tmp = []
        res2_indices_tmp = []
        for r1i, r2i in zip(res1_indices, res2_indices):
            resSeq1 = residues1[r1i].resSeq
            resSeq2 = residues2[r2i].resSeq
            if (resSeq1 < resSeq_border and resSeq2 >= resSeq_border) or (resSeq1 >= resSeq_border and resSeq2 < resSeq_border):
                res1_indices_tmp.append(r1i)
                res2_indices_tmp.append(r2i)
        res1_indices = np.array(res1_indices_tmp, dtype=np.int)
        res2_indices = np.array(res2_indices_tmp, dtype=np.int)

    saltbridges = np.zeros([len(res1_indices), Epot.shape[2]])
    xlabels     = []

    for sb_ndx, res1_ndx, res2_ndx in zip(range(len(res1_indices)), res1_indices, res2_indices):
        saltbridges[sb_ndx,:] = Epot[res1_ndx, res2_ndx, :]
        xlabels.append("{}{} - {}{}".format(residues1[res1_ndx].name, residues1[res1_ndx].resSeq, residues2[res2_ndx].name, residues2[res2_ndx].resSeq))
        #print(residues1[res1_ndx]._residue, residues2[res2_ndx]._residue)

    R = np.corrcoef(saltbridges)

    # plot
    cmap = plt.get_cmap('bwr_r')
    fig = plt.figure()
    axs = fig.add_subplot(111) 
    image = axs.imshow(R, interpolation='none', cmap=cmap, vmin=-1, vmax=1)
    axs.set_xticks(range(saltbridges.shape[0]))
    axs.set_xticklabels(xlabels, rotation="vertical")   
    axs.set_yticks(range(saltbridges.shape[0]))
    axs.set_yticklabels(xlabels, rotation="horizontal") 
    fig.colorbar(image, label='correlation coefficient')
    plt.savefig(os.path.join(plotdir, "saltbridges_corrcoef{}.png".format(suffix)), dpi=600, orientation="landscape", papertype="a5", bbox_inches="tight")
    plt.close() 



def compute_Epot_saltbridges(trj, positive_residues, negative_residues, periodic=True, opt=True, verbose=False):
    """
    Compute the electrostatic energy of each pair of positive and negative residues during the trajectory
    Return:
        Epot:       array with shape [len(positive_residues), len(negative_residues), trj.n_frames]
        ideal_dist: idealized distance of two charged particles with the residue charges and Epot
                    array with shape [len(positive_residues), len(negative_residues), trj.n_frames]
    Energies are reported in kJ/mol
    """

    Epot       = np.zeros([len(positive_residues), len(negative_residues), trj.n_frames])
    ideal_dist = np.zeros_like(Epot)

    st = time.time()
    for i, posres in enumerate(positive_residues):
        for j, negres in enumerate(negative_residues):
            Epot_residues, ideal_dist_residues = compute_Epot_residues(trj, posres, negres, periodic=periodic, opt=opt)
            Epot[i,j,:]       = Epot_residues
            ideal_dist[i,j,:] = ideal_dist_residues
    et = time.time()

    return Epot, ideal_dist




def compute_Epot_residues(trj, res1, res2, periodic=True, opt=True):
    """
    Compute the electrostatic energy of each pair of atoms in the two residues during the trajectory
    Return:
        Epot:       array that contains the energy per frame in kJ/mol
        ideal_dist: idealized distance of two charged particles with the residue charges and Epot
    """

    # no energy for interaction with itself
    if res1 == res2:
        Epot       = np.zeros(trj.n_frames)
        ideal_dist = np.zeros(trj.n_frames)
        return Epot, ideal_dist

    qq        = np.array([q1*q2 for (q1, q2) in itertools.product(res1.atom_charges, res2.atom_charges)])
    pairs     = itertools.product(res1.atom_indices, res2.atom_indices)
    distances = md.compute_distances(trj, pairs, periodic=periodic, opt=opt)

    Epot = (qq / distances).sum(1) # kJ / mol

    # Idealised distance of two charges with this energy
    qq_ideal = res1.charge * res2.charge
    if np.abs(qq_ideal) > 0:
        ideal_dist = (qq_ideal) / Epot 
    else:
        ideal_dist = np.zeros(trj.n_frames)

    # convert to gromacs units
    Epot *= ke_gu

    return Epot, ideal_dist



def compute_residue_distance(trj, residues1, residues2, periodic=True, opt=True, verbose=False):

    distance_min     = np.zeros([len(residues1), len(residues2), trj.n_frames])
    distance_mean    = np.zeros([len(residues1), len(residues2), trj.n_frames])
    distance_sc_min  = np.zeros([len(residues1), len(residues2), trj.n_frames])
    distance_sc_mean = np.zeros([len(residues1), len(residues2), trj.n_frames])

    for i, res1 in enumerate(residues1):
        for j, res2 in enumerate(residues2):
            distance    = md.compute_distances(trj, itertools.product(res1.atom_indices,      res2.atom_indices),      periodic=periodic, opt=opt)
            distance_sc = md.compute_distances(trj, itertools.product(res1.sidechain_indices, res2.sidechain_indices), periodic=periodic, opt=opt)

            distance_min    [i,j,:] = distance.min(1)
            distance_mean   [i,j,:] = distance.mean(1)
            distance_sc_min [i,j,:] = distance_sc.min(1)
            distance_sc_mean[i,j,:] = distance_sc.mean(1) 

    return distance_min, distance_mean, distance_sc_min, distance_sc_mean


class Residue(object):

    def __init__(self, residue, gmxtop=None):
        """
        Residue: mdtraj residue object
        """
        self._residue      = residue
        self._atoms        = [a for a in residue.atoms]
        self.atom_charges  = np.zeros(len(self._atoms))
        self.atom_masses   = np.zeros(len(self._atoms))

        if gmxtop is not None:
            self.read_charges(gmxtop)

    @property
    def name(self):
        return self._residue.name

    @property
    def resid(self):
        return self._residue.index

    @property
    def resSeq(self):
        return self._residue.resSeq

    @property
    def segment_id(self):
        return self._residu.segment_id

    @property
    def charge(self):
        return self.atom_charges.sum()

    @property
    def atom_indices(self):
        return [a.index for a in self._atoms]

    @property
    def sidechain_indices(self):
        return [a.index for a in self._atoms if a.is_sidechain]

    @property
    def backbone_indices(self):
        return [a.index for a in self._atoms if a.is_backbone]

    @property
    def atom_names(self):
        return [a.name for a in self._atoms] 

    @property
    def is_protein(self):
        return self._residue.is_protein

    def __eq__(self, other):
        return self._residue == other._residue
    def __ne__(self, other):
        return self._residue != other._residue

    def read_charges(self, gmxtop):
        """
        Read charges of atoms from gromacs topology
        """

        chargesAssigned = np.zeros_like(self.atom_charges, dtype=np.bool)

        # read topology from file
        with open(gmxtop) as f:
            lines = f.readlines()

        inAtomSection = False
        for ln, line in enumerate(lines):

            # start reading atoms
            if not inAtomSection and line.strip().find("[ atoms ]") == 0:
                inAtomSection = True

            # stop reading atoms
            elif inAtomSection and line.strip().find("[") == 0:
                break

            # read atoms
            elif inAtomSection and line.strip().find(";") != 0:
                fields = line.split()
                # if this line contains an atom of this residue
                if len(fields) >= 8 and fields[2] == str(self.resSeq) and fields[3] == self.name:
                    charge   = float(fields[6])
                    mass     = float(fields[7])
                    atomName = fields[4]
                    resName  = fields[3]

                    # some atom renaming needs to be done:
                    if resName in ["GLY"]:
                        if   atomName == "HA1": atomName = "HA2"
                        elif atomName == "HA2": atomName = "HA3" 
                    if resName in ["ARG", "HIS", "LYS", "ASP", "GLU", "SER", "ASN", "GLN", "CYS", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"]:
                        if   atomName == "HB1": atomName = "HB2"
                        elif atomName == "HB2": atomName = "HB3"
                    if resName in ["ARG", "LYS", "GLU", "GLN", "MET", "PRO"]:
                        if   atomName == "HG1": atomName = "HG2"
                        elif atomName == "HG2": atomName = "HG3" 
                    if resName in ["LYS", "PRO", "ARG"]:
                        if   atomName == "HD1": atomName = "HD2"
                        elif atomName == "HD2": atomName = "HD3"  
                    if resName in ["LYS"]:
                        if   atomName == "HE1": atomName = "HE2"
                        elif atomName == "HE2": atomName = "HE3"  
                    if resName in ["ILE"]:
                        if   atomName == "HG11": atomName = "HG12"
                        elif atomName == "HG12": atomName = "HG13"  
                        elif atomName == "CD":   atomName = "CD1"  
                        elif atomName == "HD1":  atomName = "HD11"  
                        elif atomName == "HD2":  atomName = "HD12"  
                        elif atomName == "HD3":  atomName = "HD13"  

                    # Termini:
                    if   atomName == "H1":  atomName = "H"  
                    elif atomName in ["OC1", "OT1"]:  atomName = "O"
                    elif atomName in ["OC2", "OT2"]:  atomName = "OXT"

                    # store charge and mass
                    try:
                        idx = self.atom_names.index(atomName)
                        if chargesAssigned[idx] == False:
                            self.atom_charges[idx]    = charge
                            self.atom_masses[idx]     = mass
                            chargesAssigned[idx] = True
                        else:
                            print("Reading atom for the second time: {}-{} {} (topology line nr. {})".format(self.name, self.resSeq, atomName, ln))
                            sys.exit(1)

                    except ValueError as e:
                        print("Cannot find atom {} from topology line number {} in".format(atomName, ln))
                        print("residue {}-{} with atoms:".format(self.name, self.resSeq))
                        for a in self.atom_names:
                            print("\t", a)
                        sys.exit(1)

        # check that all charges have been read
        if not chargesAssigned.sum() == len(chargesAssigned):
            print("Could not find these atoms of {}-{} in topology:".format(self.name, self.resSeq))
            for idx, assigned in enumerate(chargesAssigned):
                if not assigned:
                    print("\t{}".format(self.atom_names[idx]))
            sys.exit(1)

