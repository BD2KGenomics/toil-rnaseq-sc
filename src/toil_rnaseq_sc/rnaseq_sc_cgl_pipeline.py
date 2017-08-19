#!/usr/bin/env python2.7

from __future__ import print_function

import argparse
import os
import shutil
import multiprocessing
import subprocess
import sys
import textwrap
import tarfile
from urlparse import urlparse
from contextlib import closing
from string import rstrip
from glob import glob

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

import rnaseq_sc_cgl_plot_functions
from rnaseq_sc_cgl_plot_functions import run_data_analysis
from quant_to_pseudo import quant_to_pseudo
from pachterlab_post import prep_tcc_matrix

SCHEMES = ('http', 'file', 's3', 'ftp')

DEFAULT_CONFIG_NAME = 'config-toil-rnaseqsc.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-rnaseqsc.tsv'

# Kallisto output file names (incomplete, currently only for post-processing output)
TCC_MATRIX_FILENAME = 'TCC_matrix.dat'
PWISE_DIST_FILENAME = 'pwise_dist_L1.dat'
NONZERO_EC_FILENAME = 'nonzero_ec.dat'
KALLISTO_MATRIX_FILENAME = 'matrix.ec'

def formattedSchemes():
    return [scheme + '://' for scheme in SCHEMES]

# Pipeline specific functions
# TODO: Validate type, keep a list of types or something
def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    
    def validate(url):
        require(urlparse(url).scheme in SCHEMES, 'URL "{}" not valid. Schemes:{}'.format(url, SCHEMES))
    
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')
            if len(sample) != 3:
                raise UserError('Bad manifest format! Expected 3 tab separated columns, got: {}'.format(sample))

            # If a directory is passed in, use all samples in that directory
            uuid, type, url = sample
            if urlparse(url).scheme == '':
                urls = ['file://' + os.path.join(url, x) for x in os.listdir(url)]
            # If url is a tarball
            elif url.endswith('tar.gz') or url.endswith('tar'):
                validate(url)
                urls = [url]
            # If URL is a fastq or series of fastqs
            elif url.endswith('fastq.gz') or url.endswith('fastq') or url.endswith('fq.gz') or url.endswith('fq'):
                urls = url.split(',')
                [validate(x) for x in urls]
            else:
                raise UserError('URL does not have approved extension: .tar.gz, .tar, .fastq.gz, .fastq, .fq.gz, .fq')
                
            # use urls instead of url s.t. if url is not processed, error is thrown *here* instead of having a non-array passed forward through the program
            sample = (uuid, type, urls)
            samples.append(sample)
    return samples


def run_single_cell(job, sample, config):
    """
    Performs single cell analysis through the quay.io/ucsc_cgl/kallisto_sc image (which uses code from the repo:
    https://github.com/pachterlab/scRNA-Seq-TCC-prep).  Output includes TCC matrix from kallisto process.

    :param job: toil job
    :param config: configuration for toil job
    :param sample: a [UUID, url(s)] pair as constructed by parse_samples
    """
    # Common logic (for handling pre- and post- Kallisto data)
    config = argparse.Namespace(**vars(config)) # why?
    config.cores = min(config.maxCores, multiprocessing.cpu_count())
    work_dir = job.fileStore.getLocalTempDir()
    # Get input files
    uuid, type, urls = sample
    config.uuid = uuid
    # Handle kallisto output file (only works w/ one file for now)
    if type == "plot":
        filename = os.path.basename(urls[0])
        download_url(job, url=urls[0], name=filename, work_dir=work_dir)
        tar = tarfile.open(name=os.path.join(work_dir, filename))
        root_dir=rstrip(os.path.basename(urls[0]), ".tar.gz") # post, kallisto, plots folders are in this root folder, with same name as the archive
        kallisto_output = None # could just forward the kallisto output
        post_processing_output = None # same with this
        # method that, given the location of the file in the tar, writes it to the global job store
        def tarToGlobal(folder, path):
            with closing(tar.extractfile(os.path.join(root_dir, folder, path))) as file:
                data = file.read()
                with job.fileStore.writeGlobalFileStream() as (stream, id):
                    stream.write(data)
                    return id
        tcc_matrix_id = tarToGlobal("post", TCC_MATRIX_FILENAME)
        pwise_dist_l1_id = tarToGlobal("post", PWISE_DIST_FILENAME)
        nonzero_ec_id = tarToGlobal("post", NONZERO_EC_FILENAME)
        kallisto_matrix_id = tarToGlobal("post", KALLISTO_MATRIX_FILENAME)
        matrix_tsv_id = tarToGlobal("kallisto", "matrix.tsv")
        matrix_cells_id = tarToGlobal("kallisto", "matrix.cells")
    # Handle fastq file(s)
    else:
        input_location = os.path.join(work_dir, "fastq_input")
        os.mkdir(input_location)
        for url in urls:
            if url.endswith('.tar') or url.endswith('.tar.gz'):
                tar_path = os.path.join(work_dir, os.path.basename(url))
                download_url(job, url=url, work_dir=work_dir)
                subprocess.check_call(['tar', '-xvf', tar_path, '-C', input_location])
                os.remove(tar_path)
            elif url.endswith('.gz'):
                download_url(job, url=url, work_dir=input_location)
                subprocess.check_call(['gunzip', os.path.join(input_location, os.path.basename(url))])
            else:
                job.fileStore.logToMaster("Download url " + str(url))
                download_url(job, url=url, work_dir=input_location)
        # Generate configuration JSON
        with open(os.path.join(work_dir, "config.json"), 'w') as config_file:
            config_file.write(build_patcherlab_config(config))
        # Get Kallisto index
        download_url(job, url=config.kallisto_index, name='kallisto_index.idx', work_dir=work_dir)
        # Create other locations for patcherlab stuff
        os.mkdir(os.path.join(work_dir, "tcc"))
        os.mkdir(os.path.join(work_dir, "output"))
        if type == "pseudo":
            # Call docker image
            dockerCall(job, tool='quay.io/ucsc_cgl/kallisto_sc:latest', workDir=work_dir, parameters=["/data/config.json"])
        else: # quantification of quake brain-style paired end fastqs, each for a different cell
            require(type == "quant", "invalid type " + type + " found in manifest ")
            os.mkdir(os.path.join(work_dir, "quant_output"))
            # Call docker image
            dockerCall(job, tool='kallisto_sc_quant', workDir = work_dir, parameters=["/data/kallisto_index.idx", "/data/quant_output", str(config.cores), "/data/fastq_input"])
            # Consolidate abundances for the various cells
            quant_output = os.path.join(work_dir, "quant_output")
            consolidated = os.path.join(work_dir, "quant_consolidated")
            os.mkdir(consolidated)
            for output_folder in os.listdir(quant_output):
                shutil.copy(os.path.join(quant_output, output_folder, "abundance.tsv"), os.path.join(consolidated, output_folder+".tsv"))
            # quant to pseudo
            quant_to_pseudo(None, consolidated, os.path.join(work_dir, "tcc"))
            # run post-processing
            save_dir = os.path.join(work_dir, "save")
            os.mkdir(save_dir)
            prep_tcc_matrix(job, threads = config.cores, tcc_output_dir = os.path.join(work_dir, "tcc"), save_dir = save_dir) # this should be the same as specified in build_pachterlab_config. It may be worth refactoring so that these don't have to be manually synced, although there's no reason for these values to ever change and thus become desynced.
        # Irrespective of whether quant or pseudo, because of quant-to-pseudo conversion
        # Build tarfile of output
        output_files = glob(os.path.join(work_dir, "tcc", "*"))
        tarball_files(tar_name='kallisto_output.tar.gz', file_paths=output_files, output_dir=work_dir)
        kallisto_output = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'kallisto_output.tar.gz'))
        # Consolidate post-processing output
        tcc_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', TCC_MATRIX_FILENAME))
        pwise_dist_l1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', PWISE_DIST_FILENAME))
        nonzero_ec_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', NONZERO_EC_FILENAME))
        kallisto_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'tcc', 'matrix.ec'))
        post_processing_output = {
                                 TCC_MATRIX_FILENAME : tcc_matrix_id,
                                 PWISE_DIST_FILENAME : pwise_dist_l1_id,
                                 NONZERO_EC_FILENAME : nonzero_ec_id,
                                 KALLISTO_MATRIX_FILENAME : kallisto_matrix_id # technically redundant
                                 }
        # Prepare files to send to plots for SC3
        matrix_tsv_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, "tcc", "matrix.tsv"))
        matrix_cells_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, "tcc", "matrix.cells"))
    # Graphing step
    if config.generate_graphs:
        graphical_output = job.addChildJobFn(run_data_analysis, config, tcc_matrix_id, pwise_dist_l1_id,
                          nonzero_ec_id, kallisto_matrix_id, matrix_tsv_id, matrix_cells_id).rv()
        job.addFollowOnJobFn(consolidate_output, config, kallisto_output, graphical_output, post_processing_output)
    else:
        # converts to UUID name scheme and transfers to output location
        consolidate_output(job, config, kallisto_output=kallisto_output, graphical_output=None, post_processing_output=post_processing_output)


def build_patcherlab_config(config):
    """
    Builds configuration file for patcherlab code.  Parameters include:
        config.maxCores
        config.sample_idx
        config.window_min
        config.window_max
        config.barcode_length

    :param config: toil job config
    :return: JSON string
    :rtype: str
    """
    return textwrap.dedent("""
        {{
            "NUM_THREADS": {cores},
            "WINDOW": [{wmin}, {wmax}],
            "SOURCE_DIR": "/opt/single_cell/source/",
            "BASE_DIR": "/data/fastq_input/",
            "sample_idx": {idx},
            "SAVE_DIR": "/data/save/",
            "dmin": 5,
            "BARCODE_LENGTH": {barcode},
            "OUTPUT_DIR": "/data/output/",
            "kallisto":{{
                "binary": "/usr/local/bin/kallisto",
                "index": "/data/kallisto_index.idx",
                "TCC_output" : "/data/tcc/"
            }}
        }}
        """).format(cores=config.cores, wmin=config.window_min, wmax=config.window_max,
                    barcode=config.barcode_length, idx=str(config.sample_idx).replace("'", "\""))

# TODO: Take os.path.join(config.uuid, ~) and consider function-ify it
def consolidate_output(job, config, kallisto_output, graphical_output, post_processing_output):
    """
    Combines the contents of the outputs into one tarball and places in output directory or s3

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str kallisto_output: FileStoreID for Kallisto output
    :param str graphical_output: FileStoreID for output of graphing step
    :param list[str] post_processing_output: 4 files of post-processing (tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id, kallisto_matrix_id)
    """
    job.fileStore.logToMaster('Consolidating output: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    graphical_tar, kallisto_tar = None, None
    # Retrieve output file paths to consolidate
    if kallisto_output is not None:
        kallisto_tar = job.fileStore.readGlobalFile(kallisto_output, os.path.join(work_dir, 'kallisto_output.tar.gz'))
    if graphical_output is not None:
        graphical_tar = job.fileStore.readGlobalFile(graphical_output, os.path.join(work_dir, 'single_cell_plots.tar.gz'))
    # Retrive post-processing output
    job.fileStore.logToMaster("This is the value of post_processing_output: " + str(post_processing_output))
    # I/O
    out_tar = os.path.join(work_dir, config.uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    tar_list = [x for x in [graphical_tar, kallisto_tar] if x is not None]
    with tarfile.open(out_tar, 'w:gz') as f_out:
        for tar in tar_list:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == kallisto_tar:
                            tarinfo.name = os.path.join(config.uuid, 'kallisto', os.path.basename(tarinfo.name))
                        elif tar == graphical_tar:
                            tarinfo.name = os.path.join(config.uuid, 'plots', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
        # Add post processing output
        if post_processing_output is not None:
            for filename, fileID in post_processing_output.items():
                f_out.add(job.fileStore.readGlobalFile(fileID), arcname=os.path.join(config.uuid, 'post', filename))
    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(config.uuid, config.output_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(config.uuid, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[os.path.join(work_dir, config.uuid + '.tar.gz')], output_dir=config.output_dir)


def generate_config():
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

        # Generates graphs of output after completion
        generate-graphs: true
        
        # Number of clusters for spectral clustering / tSNE dimensionality reduction, stain plot
        n_clusters: 3
        
        # Parameters for SC3 (k = # of clusters, debug should be FALSE)
        min_k: 2
        max_k: 3
        use_estimated_k: TRUE
        debug: FALSE

        # The length of the barcodes
        barcode-length: 14

        # The lower threshold for expected number of cells in the experiment
        window-min: 500

        # The upper threshold for expected number of cells in the experiment
        window-max: 5000

        # The sample index used for the run as a comma-separated list of genomic strings
        sample-idx: [ATCGCTCC, CCGTACAG, GATAGGTA, TGACTAGT]

        # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
        ssec:
    """.format(scheme=formattedSchemes())[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 3 tab-separated columns: UUID, TYPE, URL(s) to sample
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   TYPE        One of the following:
        #                   pseudo (run kallisto pseudo on 10xChromium data using the Pachter Lab's code from https://github.com/pachterlab/scRNA-Seq-TCC-prep, then run clustering/plots)
        #                   quant (run kallisto quant on paired-end fastqs of the form ID_1.fastq and ID_2.fastq)
        #                   plot (run clustering/plots on the tarball that is output by this program after pseudo or quant finish)
        #   URL         A URL {scheme} pointing to the sample or a full path to a directory
        #
        #   If sample is being submitted as fastqs, provide URLs separated by a comma.
        #   Samples must have the same extension - do not mix and match gzip and non-gzipped sample pairs.
        #
        #   Sample tarballs with fastq files inside must follow one of the 4 extensions: fastq.gz, fastq, fq.gz, fq
        #
        #   Sample tarballs with pipeline output should be .tar.gz just as the pipeline automatically outputs them
        #
        #   If a full path to a directory is provided for a sample, every file inside needs to be a fastq(.gz).
        #   Do not have other superfluous files / directories inside or the pipeline will complain.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  pseudo	file:///path/to/sample.tar
        #   UUID_2  pseudo	file:///path/to/first.fq.gz,file:///path/to/second.fq.gz,file:///path/to/third.fq.gz
        #   UUID_3  pseudo	http://sample-depot.com/sample.tar
        #   UUID_4  pseudo	s3://my-bucket-name/directory/paired-sample.tar.gz
        #   UUID_5  quant	/full/path/to/directory/of/fastqs/
        #   UUID_6  plot	file:///path/to/sample.tar.gz
        #
        #   Place your samples below, one per line.
        """.format(scheme=formattedSchemes()))


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
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_NAME, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_MANIFEST_NAME, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')
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
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_NAME), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_MANIFEST_NAME), generate_manifest)

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
        require(urlparse(config.kallisto_index).scheme in SCHEMES,
                'Kallisto index in config must have the appropriate URL prefix: {}'.format(SCHEMES))
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, run_single_cell, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
