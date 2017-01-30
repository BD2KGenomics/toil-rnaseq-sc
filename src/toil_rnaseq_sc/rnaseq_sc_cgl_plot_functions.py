#!/usr/bin/env python2.7
from __future__ import print_function

import os
import pickle
import sys
from urlparse import urlparse

import numpy as np
from bd2k.util.files import mkdir_p
from sklearn import cluster,manifold
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from toil_lib.files import tarball_files, copy_files
from toil_lib.urls import s3am_upload

# Matplotlib backend nonsense
import matplotlib
if sys.platform == 'darwin':
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# source: https://github.com/pachterlab/scRNA-Seq-TCC-prep (/blob/master/notebooks/10xResults.ipynb)
def run_data_analysis(job, config, tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id, kallisto_matrix_id):
    """
    Generates graphs and plots of results.  Uploads images to savedir location.
    :param job: toil job
    :param config: toil job configuration
    :param tcc_matrix_id: jobstore location of TCC matrix (.dat)
    :param pwise_dist_l1_id: jobstore location of L1 pairwise distance (.dat)
    :param nonzero_ec_id: jobstore loation of nonzero ec (.dat)
    :param kallisto_matrix_id: id of kallisto output matrix (.ec)
    """
    # source: https://github.com/pachterlab/scRNA-Seq-TCC-prep (/blob/master/notebooks/10xResults.ipynb)
    # extract output
    job.fileStore.logToMaster('Performing data analysis')
    # read files
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(tcc_matrix_id, os.path.join(work_dir, "TCC_matrix.dat"))
    job.fileStore.readGlobalFile(pwise_dist_l1_id, os.path.join(work_dir, "pwise_dist_L1.dat"))
    job.fileStore.readGlobalFile(nonzero_ec_id, os.path.join(work_dir, "nonzero_ec.dat"))
    job.fileStore.readGlobalFile(kallisto_matrix_id, os.path.join(work_dir, 'kallisto_matrix.ec'))

    ##############################################################
    # load dataset
    with open(os.path.join(work_dir, "TCC_matrix.dat"), 'rb') as f:
        tcc_matrix = pickle.load(f)
    with open(os.path.join(work_dir, "pwise_dist_L1.dat"), 'rb') as f:
        pwise_dist_l1 = pickle.load(f)
    with open(os.path.join(work_dir, "nonzero_ec.dat"), 'rb') as f:
        nonzero_ec = pickle.load(f)

    ecfile_dir = os.path.join(work_dir, 'kallisto_matrix.ec')
    eclist = np.loadtxt(ecfile_dir, dtype=str)

    tcc = tcc_matrix.T
    T_norm = normalize(tcc_matrix, norm='l1', axis=0)
    t_normt = T_norm.transpose()

    num_of_cells = np.shape(tcc_matrix)[1]
    print("NUM_OF_CELLS =", num_of_cells)
    print("NUM_OF_nonzero_EC =", np.shape(tcc_matrix)[0])

    #################################

    EC_dict = {}
    for i in range(np.shape(eclist)[0]):
        EC_dict[i] = [int(x) for x in eclist[i, 1].split(',')]

    union = set()
    for i in nonzero_ec:
        new = [tx for tx in EC_dict[i] if tx not in union]  # filter out previously seen transcripts
        union.update(new)
    NUM_OF_TX_inTCC = len(union)
    print("NUM_OF_Transcripts =", NUM_OF_TX_inTCC)  # number of distinct transcripts in nonzero eq. classes

    ##############################################################
    # inspect

    # sort eq. classes based on size
    size_of_ec = [len(EC_dict[i]) for i in nonzero_ec]
    ec_idx = [i[0] for i in sorted(enumerate(size_of_ec), key=lambda x: x[1])]
    index_ec = np.array(ec_idx)

    ec_sort_map = {}
    nonzero_ec_srt = []  # init
    for i in range(len(nonzero_ec)):
        nonzero_ec_srt += [nonzero_ec[index_ec[i]]]
        ec_sort_map[nonzero_ec[index_ec[i]]] = i

    sumi = np.array(tcc_matrix.sum(axis=1))
    sumi_sorted = sumi[index_ec]
    total_num_of_umis = int(sumi_sorted.sum())
    total_num_of_umis_per_cell = np.array(tcc_matrix.sum(axis=0))[0, :]

    print("Total number of UMIs =", total_num_of_umis)

    #################################

    fig, ax1 = plt.subplots()
    ax1.plot(sorted(total_num_of_umis_per_cell)[::-1], 'b-', linewidth=2.0)
    ax1.set_title('UMI counts per cell')
    ax1.set_xlabel('cells (sorted by UMI counts)')
    ax1.set_ylabel('UMI counts')
    ax1.set_yscale("log", nonposy='clip')
    ax1.grid(True)
    ax1.grid(True, 'minor')
    umi_counts_per_cell = os.path.join(work_dir, "UMI_counts_per_cell.png")
    plt.savefig(umi_counts_per_cell, format='png')

    fig, ax1 = plt.subplots()
    ax1.plot(sorted(sumi.reshape(np.shape(sumi)[0]))[::-1], 'r-', linewidth=2.0)
    ax1.set_title('UMI counts per eq. class')
    ax1.set_xlabel('ECs (sorted by UMI counts)')
    ax1.set_ylabel('UMI counts')
    ax1.set_yscale("log", nonposy='clip')
    ax1.grid(True)
    ax1.grid(True, 'minor')
    umi_counts_per_class = os.path.join(work_dir, "UMI_counts_per_class.png")
    plt.savefig(umi_counts_per_class, format='png')

    cell_nonzeros = np.array(((T_norm != 0)).sum(axis=0))[0]

    fig, ax1 = plt.subplots()
    ax1.plot(total_num_of_umis_per_cell, cell_nonzeros, '.g', linewidth=2.0)
    ax1.set_title('UMI counts vs nonzero ECs')
    ax1.set_xlabel('total num of umis per cell')
    ax1.set_ylabel('total num of nonzero ecs per cell')
    ax1.set_yscale("log", nonposy='clip')
    ax1.set_xscale("log", nonposy='clip')
    ax1.grid(True)
    ax1.grid(True, 'minor')
    umi_counts_vs_nonzero_ecs = os.path.join(work_dir, "UMI_counts_vs_nonzero_ECs.png")
    plt.savefig(umi_counts_vs_nonzero_ecs, format='png')

    # TCC MEAN-VARIANCE
    #todo verify this works
    TCC_var=np.var(tcc.todense(),axis=0)
    TCC_mean=np.mean(tcc.todense(),axis=0)
    TCC_mean=np.array(TCC_mean)[0]
    TCC_var=np.array(TCC_var)[0]
    fig = plt.figure()
    N=tcc.sum()
    C=tcc.shape[0]
    ax = plt.gca()
    ax.plot(TCC_mean ,TCC_var,'.', c='blue', alpha=0.5, markeredgecolor='none')
    xlims=[0.0001,10*TCC_mean.max()]
    ax.set_xlim(xlims)
    ax.set_ylim([0.0001,10*TCC_var.max()])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(xlims, [(C-1)*(xlims[0])**2, (C-1)*(xlims[1])**2], color='g', linestyle='-', linewidth=2)
    ax.plot(xlims, [(xlims[0]), (xlims[1])], color='k', linestyle='--', linewidth=1)
    ax.set_title("TCC Mean-Variance ["+str(tcc.shape[1])+" TCCs in "+str(C)+" Cells]")
    ax.set_xlabel("mean(TCC)")
    ax.set_ylabel("var(TCC)")
    tcc_mean_variance = os.path.join(work_dir, "TCC_mean_variance.png")
    plt.savefig(tcc_mean_variance, format='png')

    ##############################################################
    # clustering

    #################################
    # t-SNE
    x_tsne = tSNE_pairwise(pwise_dist_l1)

    #################################
    # spectral clustering
    num_of_clusters = 2
    similarity_mat = pwise_dist_l1.max() - pwise_dist_l1
    labels_spectral = spectral(num_of_clusters, similarity_mat)

    spectral_clustering = stain_plot(x_tsne, labels_spectral, [], "TCC -- tSNE, spectral clustering", work_dir=work_dir,
                                     filename="spectral_clustering_tSNE")

    #################################
    # affinity propagation
    pref = -np.median(pwise_dist_l1) * np.ones(num_of_cells)
    labels_aff = AffinityProp(-pwise_dist_l1, pref, 0.5)
    np.unique(labels_aff)

    affinity_propagation_tsne = stain_plot(x_tsne, labels_aff, [], "TCC -- tSNE, affinity propagation", work_dir,
                                           "affinity_propagation_tSNE")

    #################################
    # pca
    pca = PCA(n_components=2)
    x_pca = pca.fit_transform(t_normt.todense())

    affinity_propagation_pca = stain_plot(x_pca, labels_aff, [], "TCC -- PCA, affinity propagation", work_dir,
                                          "affinity_propagation_PCA")

    # build tarfile of output plots
    output_files = [umi_counts_per_cell, umi_counts_per_class, umi_counts_vs_nonzero_ecs, tcc_mean_variance,
                    spectral_clustering, affinity_propagation_tsne, affinity_propagation_pca]
    tarball_files(tar_name='single_cell_plots.tar.gz', file_paths=output_files, output_dir=work_dir)
    # return file id for consolidation
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'single_cell_plots.tar.gz'))


def AffinityProp(D, pref, damp):
    aff = cluster.AffinityPropagation(affinity='precomputed',
                                      preference=pref, damping=damp, verbose=True)
    labels = aff.fit_predict(D)
    return labels


def spectral(k, D):
    spectral = cluster.SpectralClustering(n_clusters=k, affinity='precomputed')
    spectral.fit(D)
    labels = spectral.labels_
    return labels


def tSNE_pairwise(D):
    tsne = manifold.TSNE(n_components=2, random_state=213, metric='precomputed', n_iter=2000, verbose=1);
    X_tsne = tsne.fit_transform(D);
    return X_tsne


def stain_plot(X, labels, stain, title, work_dir, filename, filetype='png', nc=2, ax_lim=0, marksize=46):
    file_location = os.path.join(work_dir, filename + "." + filetype)
    unique_labels = np.unique(labels)
    N = len(unique_labels)
    max_value = 16581375  # 255**3
    interval = int(max_value / N)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    color = [(int(i[:2], 16) / float(255), int(i[2:4], 16) / float(255),
              int(i[4:], 16) / float(255)) for i in colors]
    i = 0;
    plt.figure(figsize=(15, 10))
    for label in unique_labels:
        ind = np.squeeze(labels == label)
        if label in stain:
            plt.scatter(X[ind, 0], X[ind, 1], c='red', s=146, edgecolor='black',
                        lw=0.5, alpha=1, marker='*', label=label)
        else:
            plt.scatter(X[ind, 0], X[ind, 1], c=color[i], s=marksize, edgecolor='lightgray',
                        lw=0.5, label=label)
        i += 1
    plt.title(title)
    plt.gray()
    plt.legend(loc='upper right', bbox_to_anchor=(1.18, 1.01), ncol=nc)
    if ax_lim > 0:
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])
    plt.axis('off')
    plt.savefig(file_location, format=filetype)
    return file_location

