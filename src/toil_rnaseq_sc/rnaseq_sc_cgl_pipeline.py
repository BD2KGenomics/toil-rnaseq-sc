#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import multiprocessing
import os
import re
import subprocess
import sys
import tarfile
import textwrap
from contextlib import closing
from subprocess import PIPE
from urlparse import urlparse

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job, PromisedRequirement
from toil_lib import flatten
from toil_lib import require, UserError
from toil_lib.files import copy_files
from toil_lib.jobs import map_job
from toil_lib.tools.QC import run_fastqc
from toil_lib.tools.aligners import run_star
from toil_lib.tools.preprocessing import run_cutadapt
from toil_lib.tools.quantifiers import run_kallisto, run_rsem, run_rsem_postprocess
from toil_lib.urls import download_url_job, s3am_upload

from toil_lib.programs import docker_call
from toil_lib.urls import download_url
from toil_lib.files import tarball_files

#analysis imports
import numpy as np
import pickle
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn import cluster,manifold
import matplotlib.pyplot as plt

from rnaseq_sc_cgl_plot_functions import *


SCHEMES = ('http', 'file', 's3', 'ftp')

DEFAULT_CONFIG_NAME = 'config-toil-rnaseqsc'
DEFAULT_CONFIG_YAML = DEFAULT_CONFIG_NAME + ".yaml"
DEFAULT_CONFIG_MANIFEST = DEFAULT_CONFIG_NAME + ".tsv"

# Pipeline specific functions
def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :param list sample_urls: Sample URLs
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')
            require(len(sample) == 4, 'Bad manifest format! '
                                      'Expected 4 tab separated columns, got: {}'.format(sample))
            file_type, paired, uuid, url = sample
            require(file_type == 'tar' or file_type == 'fq',
                    '1st column must be "tar" or "fq": {}'.format(sample[0]))
            require(paired == 'paired' or paired == 'single',
                    '2nd column must be "paired" or "single": {}'.format(sample[1]))
            if file_type == 'fq' and paired == 'paired':
                require(len(url.split(',')) == 2, 'Fastq pair requires two URLs separated'
                                                  ' by a comma: {}'.format(url))
            samples.append(sample)
    return samples


def run_single_cell(job, config, samples):
    work_dir = job.fileStore.getLocalTempDir()
    # get config for patcherlab
    with open(os.path.join(work_dir, "config.json"), 'w') as config_file:
        config_file.write(build_patcherlab_config(config))
    # get kallisto index
    download_url(url=config.kallisto_index, name='kallisto_index.idx', work_dir=work_dir)
    # get input files
    input_location = os.path.join(work_dir, "fastq_input")
    os.mkdir(input_location)
    for sample in samples:
        _, _, _, sample_location = sample
        download_url(sample_location, name=os.path.basename(sample_location), work_dir=input_location)
    # create other locations for patcherlab stuff
    os.mkdir(os.path.join(work_dir, "tcc"))
    os.mkdir(os.path.join(work_dir, "output"))

    get_data_directory_contents(work_dir)
    get_patcherlab_config_contents(work_dir)
    docker_call(tool='quay.io/ucsc_cgl/kallisto_sc:latest', work_dir=work_dir, parameters=["/data/config.json"])

    # save output for next step
    tcc_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'TCC_matrix.dat'))
    pwise_dist_l1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'pwise_dist_L1.dat'))
    nonzero_ec_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'nonzero_ec.dat'))
    kallisto_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'tcc', 'matrix.ec'))
    job.addFollowOnJobFn(run_data_analysis, config, tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id, kallisto_matrix_id)

    # build tarfile of output
    output_files = [os.path.join(work_dir, "tcc", x) for x in ['run_info.json', 'matrix.tsv', 'matrix.ec', 'matrix.cells']]
    output_files.extend([os.path.join(work_dir, "save", x) for x in ['TCC_matrix.dat', 'pwise_dist_L1.dat', 'nonzero_ec.dat']])
    tarball_files(tar_name='single_cell.tar.gz', file_paths=output_files, output_dir=work_dir)
    output_file_location = os.path.join(work_dir, 'single_cell.tar.gz')
    # save tarfile to output directory (intermediate step)
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(output_file_location, config.output_dir))
        s3am_upload(fpath=output_file_location, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_file_location, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[output_file_location], output_dir=config.output_dir)


def build_patcherlab_config(config):
    #todo: this needs to be parameterized I think
    return textwrap.dedent("""
        {
            "NUM_THREADS": %(num_threads)d,
            "WINDOW": [500, 5000],
            "SOURCE_DIR": "/opt/single_cell/source/",
            "BASE_DIR": "/data/fastq_input/",
            "sample_idx": ["ATCGCTCC","CCGTACAG","GATAGGTA","TGACTAGT"],
            "SAVE_DIR": "/data/save/",
            "dmin": 5,
            "BARCODE_LENGTH": 14,
            "OUTPUT_DIR": "/data/output/",
            "kallisto":{
                "binary": "/usr/local/bin/kallisto",
                "index": "/data/kallisto_index.idx",
                "TCC_output" : "/data/tcc/"
            }
        }
    """) % {'num_threads':1}

def run_data_analysis(job, config, tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id, kallisto_matrix_id):
    # source: https://github.com/pachterlab/scRNA-Seq-TCC-prep (/blob/master/notebooks/10xResults.ipynb)
    # extract output
    print("running data analysis")
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

    TCC = tcc_matrix.T
    T_norm = normalize(tcc_matrix, norm='l1', axis=0)
    T_normT = T_norm.transpose()

    NUM_OF_CELLS = np.shape(tcc_matrix)[1]
    print("NUM_OF_CELLS =", NUM_OF_CELLS)
    print("NUM_OF_nonzero_EC =", np.shape(tcc_matrix)[0])

    #################################

    EC_dict = {}
    for i in range(np.shape(eclist)[0]):
        EC_dict[i] = [int(x) for x in eclist[i, 1].split(',')]

    union = set()
    for i in nonzero_ec:
        new = [tx for tx in EC_dict[i] if tx not in union]  # filter out previously seen transcripts
        union.update(new)
    union_list = list(union)  # union of all transctipt ids seen in nonzero eq.classes
    NUM_OF_TX_inTCC = len(union)
    print("NUM_OF_Transcripts =", NUM_OF_TX_inTCC)  # number of distinct transcripts in nonzero eq. classes

    ##############################################################
    # inspect

    # sort eq. classes based on size
    size_of_ec = [len(EC_dict[i]) for i in nonzero_ec]
    ec_idx = [i[0] for i in sorted(enumerate(size_of_ec), key=lambda x: x[1])]
    index_ec = np.array(ec_idx)

    ec_sort_map = {}
    nonzero_ec_srt = [];  # init
    for i in range(len(nonzero_ec)):
        nonzero_ec_srt += [nonzero_ec[index_ec[i]]]
        ec_sort_map[nonzero_ec[index_ec[i]]] = i
    nonzero_ec_srt = np.array(nonzero_ec_srt)

    ec_size_sort = np.array(size_of_ec)[index_ec]
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
    # ax1.set_xscale("log", nonposy='clip')
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
    # ax1.set_xscale("log", nonposy='clip')
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
    umi_counts_vs_nonzero_ECs = os.path.join(work_dir, "UMI_counts_vs_nonzero_ECs.png")
    plt.savefig(umi_counts_vs_nonzero_ECs, format='png')

    #################################

    # TCC MEAN-VARIANCE

    meanvar_plot(TCC, alph=0.5)

    ##############################################################
    # clustering

    #################################
    # t-SNE
    X_tsne = tSNE_pairwise(pwise_dist_l1)

    #################################
    # spectral clustering
    num_of_clusters = 2
    similarity_mat = pwise_dist_l1.max() - pwise_dist_l1
    labels_spectral = spectral(num_of_clusters, similarity_mat)

    spectral_clustering = stain_plot(X_tsne, labels_spectral, [], "TCC -- tSNE, spectral clustering", work_dir=work_dir, filename="spectral_clustering_tSNE")

    #################################
    # affinity propagation
    pref = -np.median(pwise_dist_l1) * np.ones(NUM_OF_CELLS)
    labels_aff = AffinityProp(-pwise_dist_l1, pref, 0.5)
    np.unique(labels_aff)

    affinity_propagation_tsne = stain_plot(X_tsne, labels_aff, [], "TCC -- tSNE, affinity propagation", work_dir, "affinity_propagation_tSNE")

    #################################
    # pca
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(T_normT.todense())

    affinity_propagation_pca = stain_plot(X_pca, labels_aff, [], "TCC -- PCA, affinity propagation", work_dir, "affinity_propagation_PCA")


    # build tarfile of output plots
    output_files = [umi_counts_per_cell, umi_counts_per_class, umi_counts_vs_nonzero_ECs,
                    spectral_clustering, affinity_propagation_tsne, affinity_propagation_pca]
    tarball_files(tar_name='single_cell_plots.tar.gz', file_paths=output_files, output_dir=work_dir)
    output_file_location = os.path.join(work_dir, 'single_cell_plots.tar.gz')
    # save tarfile to output directory (intermediate step)
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(output_file_location, config.output_dir))
        s3am_upload(fpath=output_file_location, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(output_file_location, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[output_file_location], output_dir=config.output_dir)

def get_kallisto_version(job, sample, config):
    # verifying kallisto version
    docker_call(tool='quay.io/ucsc_cgl/kallisto_sc:latest',
                parameters=['version'], docker_parameters=["--entrypoint=ls"])

def get_data_directory_contents(work_dir):
    docker_call(tool='quay.io/ucsc_cgl/kallisto_sc:latest',
                work_dir=work_dir, parameters=["-laR", "/data"], docker_parameters=["--entrypoint=ls"])

def get_patcherlab_config_contents(work_dir):
    docker_call(tool='quay.io/ucsc_cgl/kallisto_sc:latest',
                work_dir=work_dir, parameters=["/data/config.json"], docker_parameters=["--entrypoint=cat"])

def generate_config():
    #todo TPESOUT: verify this is right
    return textwrap.dedent("""
        # RNA-seq CGL Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-rnaseq-sc run"
        # Just Kallisto or STAR/RSEM can be run by supplying only the inputs to those tools
        #
        # URLs can take the form: http://, ftp://, file://, s3://, gnos://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################
        # Required: URL {scheme} to kallisto index file.
        kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
        ssec:

        # Optional: If true, uses resource requirements appropriate for continuous integration
        ci-test:
    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    #todo I think this format isn't quite right
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 4 tab-separated columns: filetype, paired/unpaired, UUID, URL(s) to sample
        #
        #   filetype    Filetype of the sample. Options: "tar" or "fq", for tarball/tarfile or fastq/fastq.gz
        #   paired      Indicates whether the data is paired or single-ended. Options:  "paired" or "single"
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL {scheme} pointing to the sample
        #
        #   If sample is being submitted as a fastq pair, provide two URLs separated by a comma.
        #   Samples must have the same extension - do not mix and match gzip and non-gzipped sample pairs.
        #
        #   Samples consisting of tarballs with fastq files inside must follow the file name convention of
        #   ending in an R1/R2 or _1/_2 followed by one of the 4 extensions: .fastq.gz, .fastq, .fq.gz, .fq
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   tar paired  UUID_1  file:///path/to/sample.tar
        #   fq  paired  UUID_2  file:///path/to/R1.fq.gz,file:///path/to/R2.fq.gz
        #   tar single  UUID_3  http://sample-depot.com/single-end-sample.tar
        #   tar paired  UUID_4  s3://my-bucket-name/directory/paired-sample.tar.gz
        #   fq  single  UUID_5  s3://my-bucket-name/directory/single-end-file.fq
        #
        #   Place your samples below, one per line.
        """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))

def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil RNA-seq single cell pipeline

    Imagine this is filled out

    =======================================
    Dependencies
    Curl:       apt-get install curl
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL)
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the RNA-seq single cell pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_YAML, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_CONFIG_MANIFEST, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')
    group.add_argument('--samples', default=None, nargs='+', type=str,
                       help='Space delimited sample URLs (any number). Samples must be tarfiles/tarballs that contain '
                            'fastq files. URLs follow the format: http://foo.com/sample.tar, '
                            'file:///full/path/to/file.tar. The UUID for the sample will be derived from the file.'
                            'Samples passed in this way will be assumed to be paired end, if using single-end data, '
                            'please use the manifest option.')
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_YAML), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_MANIFEST), generate_manifest)
    # Pipeline execution
    elif args.command == 'run':
        # sanity check
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-rnaseq generate-config"'.format(args.config))
        require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                               'run "toil-rnaseq generate-manifest"'.format(args.manifest))
        # get samples
        samples = parse_samples(path_to_manifest=args.manifest)
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint
        # Config sanity checks
        require(config.kallisto_index,
                'URLs not provided for Kallisto index, so there is nothing to do!')
        require(config.output_dir, 'No output location specified: {}'.format(config.output_dir))
        for input in [x for x in [config.kallisto_index] if x]:
            require(urlparse(input).scheme in SCHEMES,
                    'Input in config must have the appropriate URL prefix: {}'.format(SCHEMES))
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow by using map_job() to run the pipeline for each sample
        # Job.Runner.startToil(Job.wrapJobFn(map_job, get_kallisto_version, samples, config), args)
        # tpesout: start the job without map_job
        Job.Runner.startToil(Job.wrapJobFn(run_single_cell, config, samples), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
