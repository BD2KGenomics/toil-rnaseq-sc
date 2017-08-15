# This is a modified version of a source file from the repository "scRNA-Seq-tcc-prep" by the Pachter Lab which can be found here: https://github.com/pachterlab/scRNA-Seq-TCC-prep/blob/0469873bdadcc48e34782882dbd24c3939c0542a/source/prep_TCC_matrix.py
# The citation for the paper with which this repository is associated is Ntranos, V., Kamath, G. M., Zhang, J. M., Pachter, L. & Tse, D. N. Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts. Genome Biology 17, 112 (2016).
# The entire source of "scRNA-Seq-tcc prep" is also used in Dockerized form in this pipeline.
import os
import sys, gc

import numpy as np
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize

from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import *
from scipy.stats import entropy

import pickle

def prep_tcc_matrix(job, threads, tcc_output_dir, save_dir):
    """
    For some reason, changing the number of threads to more than one results in a crash.
    """
    print "Setting threads to 1... threads value was ignored."
    threads = 1
    # matrix.ec file
    ecfile_dir = os.path.join(tcc_output_dir, "matrix.ec")
    tsvfile_dir = os.path.join(tcc_output_dir, "matrix.tsv")

    print "Loading TCCs.."

    COOinput = np.loadtxt( tsvfile_dir, delimiter='\t' , dtype=float)
    rows,cols,data = COOinput.T
    nonzero_ec = np.unique(rows)
    map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
    map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
    TCCmatrix   = coo_matrix( (data.astype(float),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) )

    NUM_OF_CELLS = TCCmatrix.shape[1]
    print "NUM_OF_CELLS =", NUM_OF_CELLS
    
    T = TCCmatrix.tocsr()
    T_norm = normalize(T, norm='l1', axis=0)
    T_normT = T_norm.transpose()
    del TCCmatrix;
    _ = gc.collect()
    
    # Pairwise_distances
    def L1_distance(p,q):
        return cityblock(p,q).sum()

    # def jensen_shannon(p, q):
    #     m=0.5*p+0.5*q
    #     p = np.transpose(p[p > 0])
    #     q = np.transpose(q[q > 0])
    #     m = np.transpose(m[m > 0])
    #     return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))

    num_of_threads = threads
    print "Calculating pairwise L1 distances... ( num_threads =",num_of_threads,")"

    # D_js = pairwise_distances(T_normT,metric=jensen_shannon,n_jobs=num_of_threads)
    D_l1 = pairwise_distances(T_normT,metric=L1_distance,n_jobs=num_of_threads)

    print "writing data..."

    # Save data
    with open(os.path.join(save_dir, "TCC_matrix.dat"), 'wb') as f:
        pickle.dump(T,f)
    with open(os.path.join(save_dir, "pwise_dist_L1.dat"), 'wb') as f:
        pickle.dump(D_l1,f)
    with open(os.path.join(save_dir, "nonzero_ec.dat"), 'wb') as f:
        pickle.dump(nonzero_ec,f)

    print "DONE."
