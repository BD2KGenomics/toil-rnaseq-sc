#!/usr/bin/env python2.7
from __future__ import print_function

import os
import numpy as np
import pickle
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn import cluster,manifold
import matplotlib.pyplot as plt


# source: https://github.com/pachterlab/scRNA-Seq-TCC-prep (/blob/master/notebooks/10xResults.ipynb)

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


def stain_plot(X, labels, stain, title, work_dir, filename=None, filetype='png', nc=2, ax_lim=0, marksize=46):
    #todo do this better
    if filename is None:
        filename = title
    file_location = os.path.join(work_dir, filename + "." + filetype)
    #/todo
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


# Plot function with colors(heatmap) corresponding to the total number of umis per cell  (log10)
def umi_counts_plot(X_tsne, total_num_of_umis_per_cell, title, ax_lim=0, marksize=26):
    plt.figure(figsize=(15, 10))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=np.log10(total_num_of_umis_per_cell), s=marksize, edgecolor='black',
                lw=0.25, cmap='OrRd')
    plt.colorbar()
    plt.title(title)
    if ax_lim > 0:
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])
    plt.axis('off')


def color_plot(X, colorvec, title, ax_lim=0, marksize=26):
    plt.figure(figsize=(10, 6))
    plt.scatter(X[:, 0], X[:, 1], c=np.log10(1 + colorvec), s=marksize, edgecolor='black', lw=0.25, cmap='OrRd')
    plt.colorbar()
    plt.title(title)
    if ax_lim > 0:
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])
    plt.axis('off')
    plt.show()

## TCC MEAN-VARIANCE
def meanvar_plot(TCC_, alph=0.05):
    TCC_var = np.var(TCC_.todense(), axis=0)
    TCC_mean = np.mean(TCC_.todense(), axis=0)
    TCC_mean = np.array(TCC_mean)[0]
    TCC_var = np.array(TCC_var)[0]

    fig = plt.figure()
    N = TCC_.sum()
    C = TCC_.shape[0]
    ax = plt.gca()

    ax.plot(TCC_mean, TCC_var, '.', c='blue', alpha=alph, markeredgecolor='none')

    xlims = [0.0001, 10 * TCC_mean.max()]
    ax.set_xlim(xlims)
    ax.set_ylim([0.0001, 10 * TCC_var.max()])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(xlims, [(C - 1) * (xlims[0]) ** 2, (C - 1) * (xlims[1]) ** 2], color='g', linestyle='-', linewidth=2)
    ax.plot(xlims, [(xlims[0]), (xlims[1])], color='k', linestyle='--', linewidth=1)
    ax.set_title("TCC Mean-Variance [" + str(TCC_.shape[1]) + " TCCs in " + str(C) + " Cells]")
    ax.set_xlabel("mean(TCC)")
    ax.set_ylabel("var(TCC)")
