#!/usr/bin/env python

import nmr

nmr.bondvec_corr_batch_mpi("bla.pdb", [str(i) for i in range(26)], subtrjlength=20000.0, savepath="bla")
