#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import sys
import textwrap

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload
from urlparse import urlparse
import os

from rnaseq_sc_cgl_plot_functions import run_data_analysis

SCHEMES = ('http', 'file', 's3', 'ftp')

DEFAULT_CONFIG_NAME = 'config-toil-rnaseqsc.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-rnaseqsc.tsv'


# Pipeline specific functions
def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')
            if len(sample) != 2:
                raise UserError('Bad manifest format! Expected 2 tab separated columns, got: {}'.format(sample))

            # If a directory is passed in, use all samples in that directory
            uuid, url = sample
            if urlparse(url).scheme == '':
                url = ['file://' + os.path.join(url, x) for x in os.listdir(url)]
            # If url is a tarball
            elif url.endswith('tar.gz') or url.endswith('tar'):
                require(urlparse(url).scheme in SCHEMES, 'URL "{}" not valid. Schemes:{}'.format(url, SCHEMES))
                url = [url]
            # If URL is a fastq or series of fastqs
            elif url.endswith('fastq.gz') or url.endswith('fastq') or url.endswith('fq.gz') or url.endswith('fq'):
                url = url.split(',')
                [require(urlparse(x).scheme in SCHEMES,
                         'URL "{}" not valid. Schemes:{}'.format(url, SCHEMES)) for x in url]
            else:
                raise UserError('URL is ')

            sample = [uuid, url]
            samples.append(sample)
    return samples


def run_single_cell(job, config, samples):
    """
    Performs single cell analysis through the quay.io/ucsc_cgl/kallisto_sc image (which uses code from the repo:
    https://github.com/pachterlab/scRNA-Seq-TCC-prep).  Output includes TCC matrix from kallisto process.
    :param job: toil job
    :param config: configuration for toil job
    :param samples: list of samples as constucted by 'parse_samples' method
    """
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
    # call docker image
    docker_call(tool='quay.io/ucsc_cgl/kallisto_sc:latest', work_dir=work_dir, parameters=["/data/config.json"])
    # save output for graphing step
    if config.generateGraphs:
        tcc_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'TCC_matrix.dat'))
        pwise_dist_l1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'pwise_dist_L1.dat'))
        nonzero_ec_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'save', 'nonzero_ec.dat'))
        kallisto_matrix_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'tcc', 'matrix.ec'))
        job.addFollowOnJobFn(run_data_analysis, config, tcc_matrix_id, pwise_dist_l1_id, nonzero_ec_id,
                             kallisto_matrix_id)
    # build tarfile of output
    output_files = [os.path.join(work_dir, "tcc", x) for x in ['run_info.json', 'matrix.tsv', 'matrix.ec',
                                                               'matrix.cells']]
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
    """
    Builds configuration file for patcherlab code.  Parameters include:
        config.maxCores
        config.sampleIdx
        config.windowMin
        config.windowMax
        config.barcodeLength
    :param config: toil job config
    :return: JSON string to send to patcherlab code
    """
    return textwrap.dedent("""
        {
            "NUM_THREADS": %(num_threads)d,
            "WINDOW": [%(window_min)d, %(window_max)d],
            "SOURCE_DIR": "/opt/single_cell/source/",
            "BASE_DIR": "/data/fastq_input/",
            "sample_idx": %(sample_idx)s,
            "SAVE_DIR": "/data/save/",
            "dmin": 5,
            "BARCODE_LENGTH": %(barcode_length)d,
            "OUTPUT_DIR": "/data/output/",
            "kallisto":{
                "binary": "/usr/local/bin/kallisto",
                "index": "/data/kallisto_index.idx",
                "TCC_output" : "/data/tcc/"
            }
        }
    """) % {'num_threads':config.maxCores,'window_min':config.windowMin,'window_max':config.windowMax,'barcode_length':
            config.barcodeLength, 'sample_idx':'["' + '","'.join(config.sampleIdx) + '"]'}



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
        kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/transcript.idx

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Generates graphs of output after completion
        generate-graphs: true

        # The length of the barcodes
        barcode-length: 14

        # The lower threshold for expected number of cells in the experiment
        window-min: 500

        # The upper threshold for expected number of cells in the experiment
        window-max: 5000

        # The sample index used for the run as a comma-separated list of genomic strings
        sample-idx: [ACAGCAAC, CGCAATTT, GAGTTGCG, TTTCGCGA]

        # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
        ssec:

        # Optional: If true, uses resource requirements appropriate for continuous integration
        ci-test:
    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2 tab-separated columns: UUID, URL(s) to sample
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL {scheme} pointing to the sample or a full path to a directory
        #
        #   If sample is being submitted as a fastqs, provide URLs separated by a comma.
        #   Samples must have the same extension - do not mix and match gzip and non-gzipped sample pairs.
        #
        #   Samples consisting of tarballs with fastq files inside must follow the file name convention of
        #   ending in an R1/R2 or _1/_2 followed by one of the 4 extensions: .fastq.gz, .fastq, .fq.gz, .fq
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/sample.tar
        #   UUID_2  file:///path/to/first.fq.gz,file:///path/to/second.fq.gz,file:///path/to/third.fq.gz
        #   UUID_3  http://sample-depot.com/sample.tar
        #   UUID_4  s3://my-bucket-name/directory/paired-sample.tar.gz
        #   UUID_5  /full/path/to/directory/of/fastqs/
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
                'Kallisto in config must have the appropriate URL prefix: {}'.format(SCHEMES))
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, run_single_cell, config, samples), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
