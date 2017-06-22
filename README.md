# **Warning: Some of this documentation is inaccurate**

Test commit 1

## University of California, Santa Cruz Genomics Institute
### Guide: Running the Single Cell RNA-seq Pipeline using Toil

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request. Feedback of any kind is appreciated.

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Inputs](#inputs)
- [Usage](#general-usage)
- [Methods](#methods)


## Overview

RNA-seq fastqs generated from 10x Chromium single-cell experiments are quantified to produce a gene by cell matrix.
Additional QC plots are generated 

This pipeline produces a tarball (tar.gz) file for a given sample that contains n subdirectories:

- 
- 
- 

The output tarball is prepended with the UUID for the sample (e.g. UUID.tar.gz). 

# Dependencies

This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have `sudo` 
privileges you will need to build these tools from source, or bug a sysadmin about how to get them (they don't mind). 

#### General Dependencies

    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       http://docs.docker.com/engine/installation/

#### Python Dependencies

    1. Toil         pip install toil
    2. S3AM         pip install --pre s3am (optional, needed for uploading output to S3)
    
    
#### System Dependencies


# Installation

 
# Inputs

The CGL RNA-seq pipeline requires an index file in order to run. This file is hosted on Synapse and can 
be downloaded after creating an account which takes about 1 minute and is free. 

* Register for a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)
* Either download the samples from the [website GUI](https://www.synapse.org/#!Synapse:syn5886029) or use the Python API
* `pip install synapseclient`
* `python`
    * `import synapseclient`
    * `syn = synapseclient.Synapse()`
    * `syn.login('foo@bar.com', 'password')`
    * Get the Kallisto index reference
        * `syn.get('syn5889216', downloadLocation='.')`
        
 
All samples and inputs must be submitted as URLs with support for the following schemas: 
`http://`, `file://`, `s3://`, `ftp://`.

Samples consisting of tarballs with fastq files inside _must_ follow the file name convention of ending in an 
R1/R2 or \_1/\_2 followed by `.fastq.gz`, `.fastq`, `.fq.gz` or `.fq.`.

# General Usage

Type `toil-rnaseq` to get basic help menu and instructions
 
1. Type `toil-rnaseq-sc generate` to create an editable manifest and config in the current working directory.
2. Parameterize the pipeline by editing the config.
3. Fill in the manifest with information pertaining to your samples.
4. Type `toil-rnaseq-sc run [jobStore]` to execute the pipeline.

### Example Commands

Run sample(s) locally using the manifest

1. `toil-rnaseq-sc generate`
2. Fill in config and manifest
3. `toil-rnaseq-sc run ./example-jobstore`

Toil options can be appended to `toil-rnaseq run`, for example:
`toil-rnaseq-sc run ./example-jobstore --retryCount=1 --workDir=/data`

For a complete list of Toil options, just type `toil-rnaseq run -h`

Run a variety of samples locally

1. `toil-rnaseq-sc generate-config`
2. Fill in config
3. `toil-rnaseq-sc run ./example-jobstore --retryCount=1 --workDir=/data --samples \
    s3://example-bucket/sample_1.tar file:///full/path/to/sample_2.tar https://sample-depot.com/sample_3.tar`

### Example Config

```
kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx
output-dir: /data/my-toil-run
ssec: 
ci-test:
```


## Distributed Run

To run on a distributed AWS cluster, see [CGCloud](https://github.com/BD2KGenomics/cgcloud) for instance provisioning, 
then run `toil-rnaseq-sc run aws:us-west-2:example-jobstore-bucket --batchSystem=mesos --mesosMaster mesos-master:5050`
to use the AWS job store and mesos batch system. 

# Methods

## Tools
    
All tool containers can be found on our [quay.io account](quay.io/organization/ucsc_cgl).

## Reference Data


## Tool Options
