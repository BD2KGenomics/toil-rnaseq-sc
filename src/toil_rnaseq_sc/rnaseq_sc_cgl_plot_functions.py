#!/usr/bin/env python2.7

# This is a modified version of a source file from the repository "scRNA-Seq-tcc-prep" by the Pachter Lab which can be found here: https://github.com/pachterlab/scRNA-Seq-TCC-prep/blob/201469940e138c2f09bcd058a6291b17794f7c88/notebooks/10xResults.ipynb
# The citation for the paper with which this repository is associated is Ntranos, V., Kamath, G. M., Zhang, J. M., Pachter, L. & Tse, D. N. Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts. Genome Biology 17, 112 (2016).
# The entire source of "scRNA-Seq-tcc prep" is also used in Dockerized form in this pipeline.
# The original "scRNA-Seq-TCC-prep" repository was released under GPLv3, as is this repository (and thus this source file). For more details, see the 'README.md' of this repository which contains the full text of the GPL.

from __future__ import print_function

import os
import pickle
import sys
from subprocess import CalledProcessError
from urlparse import urlparse

import numpy as np
from bd2k.util.files import mkdir_p
from sklearn import cluster,manifold
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from toil.lib.docker import dockerCall
from toil_lib.files import tarball_files, copy_files
from toil_lib.urls import s3am_upload

from string import lstrip

# Matplotlib backend nonsense
import matplotlib
if sys.platform == 'darwin':
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

SC3_OUTPUT_DIRECTORY = "SC3"
MATRIX_TSV_FILENAME = "matrix.tsv"
MATRIX_CELLS_FILENAME = "matrix.cells"
DOCKER_WORK_DIR = "/data"

# TODO: Refactor to use ids
def run_data_analysis(job, config, tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id, kallisto_matrix_id, matrix_tsv_id, matrix_cells_id):
    """
    Generates graphs and plots of results.  Uploads images to savedir location.
    :param job: toil job
    :param config: toil job configuration
    :param tcc_matrix_id: jobstore location of TCC matrix (.dat)
    :param pwise_dist_l1_id: jobstore location of L1 pairwise distance (.dat)
    :param nonzero_ec_id: jobstore loation of nonzero ec (.dat)
    :param kallisto_matrix_id: id of kallisto output matrix (.ec)
    :param matrix_tsv_id: id of kallisto output matrix (.tsv)
    :param matrix_cells_id: id of kallisto output matrix (.cells)
    """
    # source: https://github.com/pachterlab/scRNA-Seq-TCC-prep (/blob/master/notebooks/10xResults.ipynb)
    # extract output
    job.fileStore.logToMaster('Performing data analysis')
    # read files
    work_dir = job.fileStore.getLocalTempDir()
    tcc_matrix = job.fileStore.readGlobalFile(tcc_matrix_id, os.path.join(work_dir, "TCC_matrix.dat"))
    pwise_dist_l1 = job.fileStore.readGlobalFile(pwise_dist_l1_id, os.path.join(work_dir, "pwise_dist_L1.dat"))
    nonzero_ec = job.fileStore.readGlobalFile(nonzero_ec_id, os.path.join(work_dir, "nonzero_ec.dat"))
    kallisto_matrix = job.fileStore.readGlobalFile(kallisto_matrix_id, os.path.join(work_dir, 'kallisto_matrix.ec'))
    matrix_tsv = job.fileStore.readGlobalFile(matrix_tsv_id, os.path.join(work_dir, MATRIX_TSV_FILENAME))
    matrix_cells = job.fileStore.readGlobalFile(matrix_cells_id, os.path.join(work_dir, MATRIX_CELLS_FILENAME))
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
    ax.set_yscale('symlog')
    ax.set_xscale('symlog')
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
    x_tsne = tSNE_pairwise(2, pwise_dist_l1)

    #################################
    # spectral clustering
    n_clusters = config.n_clusters
    similarity_mat = pwise_dist_l1.max() - pwise_dist_l1
    labels_spectral = spectral(n_clusters, similarity_mat)

    spectral_clustering = stain_plot(x_tsne, labels_spectral, [], "TCC -- tSNE, spectral clustering with " + str(n_clusters) + " n_clusters", work_dir=work_dir,
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

    # SC3
    outfilePath = job.fileStore.getLocalTempFile()
    SC3OutputPath = os.path.join(work_dir, SC3_OUTPUT_DIRECTORY)
    os.mkdir(SC3OutputPath)
    shouldUseSC3Output = True
    with open(outfilePath, "r+") as outfile:
        def dockerPathTo(resource): return os.path.join(DOCKER_WORK_DIR, resource)
        def boolForR(aBool): return "TRUE" if aBool else "FALSE"
        try:
            dockerCall(job, tool='rscript', workDir=work_dir, parameters=map(str, [config.min_k, config.max_k, dockerPathTo(MATRIX_TSV_FILENAME), dockerPathTo(MATRIX_CELLS_FILENAME), dockerPathTo(SC3_OUTPUT_DIRECTORY), boolForR(config.use_estimated_k), boolForR(config.debug)]), outfile=outfile)
            pass
        except CalledProcessError:
            outfile.seek(0, 0)
            job.fileStore.logToMaster("Docker failed with the following log:  " + str(outfile.read()))
            shouldUseSC3Output = False
    # build tarfile of output plots
    output_files = [umi_counts_per_cell, umi_counts_per_class, umi_counts_vs_nonzero_ecs, tcc_mean_variance,
                    spectral_clustering, affinity_propagation_tsne, affinity_propagation_pca, outfilePath] + ([os.path.join(work_dir, SC3_OUTPUT_DIRECTORY, x) for x in os.listdir(SC3OutputPath)] if shouldUseSC3Output else [])
    tarball_files(tar_name='single_cell_plots.tar.gz', file_paths=output_files, output_dir=work_dir)
    # return file id for consolidation
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'single_cell_plots.tar.gz'))
    
def AffinityProp(D, pref, damp):
    """
    Perform SKLearn affinity propagation (clustering) with specified data and parameters, returning labels.
    :param pref: preference parameter for the affinity propagation
    :param damp: damping parameter for the affinity propagation
    :return: labels
    """
    aff = cluster.AffinityPropagation(affinity='precomputed',
                                      preference=pref, damping=damp, verbose=True)
    labels = aff.fit_predict(D)
    return labels
    
def spectral(n, D):
    """
    Perform spectral clustering on the distance matrix.
    :param n: Number of clusters (for some reason, this may not equal the number displayed on the stain plot?)
    :param D: Distance matrix to analyze
    :return: labels from the spectral clustering
    """
    spectral = cluster.SpectralClustering(n_clusters=n, affinity='precomputed')
    spectral.fit(D)
    labels = spectral.labels_
    return labels
    
def tSNE_pairwise(n, D):
    """
    Perform t-SNE dimensionality reduction on the distance matrix D, using n components.
    :param n: the number of components to use (passed as n_components to sklearn.manifold.TSNE.__init__)
    :param D: Distance matrix to be processed
    :return: t-SNE reduced version of D
    """
    tsne = manifold.TSNE(n_components=n, random_state=213, metric='precomputed', n_iter=2000, verbose=1);
    X_tsne = tsne.fit_transform(D);
    return X_tsne
    
def stain_plot(X, labels, stain, title, work_dir, filename, filetype='png', nc=2, ax_lim=0, marksize=46):
    """
    Create a matplotlib plot from the specified parameters, including cluster labels and dimensionally reduced points to plot
    :param X: the reduced matrix returned by a dimensionality reduction routine e.g. tSNE or PCA
    :param labels: the labels to use to group the points into clusters
    :param stain: labels to stain
    :param title: plot title
    :param work_dir: working directory to create the file
    :param filename: name of the file to be saved to work_dir
    :param filetype: extension of the created file
    :param nc: number of columns in the legend
    :param ax_lim: limits of x- and y- axes (e.g. ax_lim = 3 -> [-3, 3] x [-3, 3] bounding box)
    :param marksize: size of the scatter-plot points that are NOT stained (stained are always 146
    """
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

