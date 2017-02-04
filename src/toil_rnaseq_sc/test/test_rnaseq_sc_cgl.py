import logging
import os
import textwrap
import subprocess
from toil_lib import require, UserError
from toil_lib.files import tarball_files
from toil_rnaseq_sc.rnaseq_sc_cgl_pipeline import parse_samples

log = logging.getLogger(__name__)


def test_parse_samples(tmpdir):
    ### helpers ###
    # specify location
    # creates manifest as tsv

    ### verifying cases work as expected ###
    # specific files - all file types and schemes

    SCHEMES = ('http', 'file', 's3', 'ftp')
    in1 = [["test1", "http://www.com.org/file.fastq.gz,file:///example/file.fastq,s3:///example/file.fq.gz,ftp://host/file.fq"]]
    manifest_location = _generate_manifest(tmpdir, in1)
    out1 = parse_samples(manifest_location)
    require(len(out1) == 1, "expected to have single output for test1")
    uuid1, url1 = out1[0]
    require(uuid1 == "test1", "expected uuid to be 'test1'")
    require(len(url1) == 4, "expected to have four files specified in url1")
    
    # tarball
    in2 = [["test2", "file:///example/tarball.tar"]]
    _generate_manifest(tmpdir, in2)
    out2 = parse_samples(manifest_location)
    require(len(out2) == 1, "expected to have single output for test2")
    uuid2, url2 = out2[0]
    require(uuid2 == "test2", "expected uuid to be 'test2'")
    require(len(url2) == 1, "expected to have two files specified in url2")

    # s3 zipped tarball
    in3 = [["test3", "file:///example/tarball.tar.gz"]]
    _generate_manifest(tmpdir, in3)
    out3 = parse_samples(manifest_location)
    require(len(out3) == 1, "expected to have single output for test3")
    uuid3, url3 = out3[0]
    require(uuid3 == "test3", "expected uuid to be 'test3'")
    require(len(url3) == 1, "expected to have two files specified in url3")

    # directory on filesystem
    test4_location = os.path.join(str(tmpdir), "test4")
    os.mkdir(test4_location)
    open(os.path.join(test4_location, "test4.1.fastq.gz"), 'a').close()
    open(os.path.join(test4_location, "test4.2.fastq"), 'a').close()
    open(os.path.join(test4_location, "test4.3.fq.gz"), 'a').close()
    open(os.path.join(test4_location, "test4.4.fq"), 'a').close()
    in4 = [["test4", test4_location]]
    _generate_manifest(tmpdir, in4)
    out4 = parse_samples(manifest_location)
    require(len(out4) == 1, "expected to have single output for test4")
    uuid4, url4 = out4[0]
    require(uuid4 == "test4", "expected uuid to be 'test4'")
    require(len(url4) == 4, "expected to have four files specified in url4")

    # all inputs
    in5 = [in1[0], in2[0], in3[0], in4[0]]
    _generate_manifest(tmpdir, in5)
    out5 = parse_samples(manifest_location)
    require(len(out5) == 4, "expected to have four outputs for test5")

    ### verifying cases which fail as expected ###
    # only uuid specified
    _generate_manifest(tmpdir, [["only_uuid"]])
    try:
        parse_samples(manifest_location)
        require(False, "parse_samples should have failed due to malformed manifest")
    except UserError:
        pass

    # multiple samples specified
    _generate_manifest(tmpdir, [["uuid", "file:///file1.fq", "file:///file2.fq"]])
    try:
        parse_samples(manifest_location)
        require(False, "parse_samples should have failed due to malformed manifest")
    except UserError:
        pass

    _generate_manifest(tmpdir, [["uuid", "badscheme:///file.fq"]])
    try:
        parse_samples(manifest_location)
        require(False, "parse_samples should have failed due to bad scheme")
    except UserError:
        pass

    require(False, "how to get assertRaises working?")

def test_pipeline_output_with_graphs(tmpdir):

    uuid = "test_rnaseqsc_g"
    output_dir = os.path.join(str(tmpdir), uuid)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, uuid + ".tar.gz")

    jobstore = os.path.join(str(tmpdir), uuid + "_jobstore")

    input = "file://" + _get_test_fastq_files(tmpdir, True)
    config = _generate_config(tmpdir, True)
    manifest = _generate_manifest(tmpdir, [[uuid, input]])

    subprocess.check_call(['toil-rnaseq-sc', 'run',
                                '--config', config,
                                '--manifest', manifest,
                                '--retryCount', '1',
                                jobstore])
    # ensure file and directories exist
    require(os.path.isfile(output_file), "expected outputfile to exist: " + output_file)
    subprocess.check_call(['tar', '-xvf', output_file, '-C', output_dir])
    require(os.path.isfile(os.path.join(output_dir, "matrix.tsv")))
    require(os.path.isdir(os.path.join(output_dir, "plots")))




    # get pachterlab fastq files
    # give directory location for manifest?
    # generate config?
    # run the pipeline
    pass


def _generate_manifest(tmpdir, list_of_lines):
    manifest_location = os.path.join(str(tmpdir), "manifest-toil-rnaseqsc-test.tsv")
    # clear old manifest if it exists
    if os.path.isfile(manifest_location):
        os.remove(manifest_location)
    with open(manifest_location, 'w') as manifest:
        for line in list_of_lines:
            manifest.write("\t".join(line) + "\n")
    return manifest_location



def _get_test_fastq_files(tmpdir, tarball=True):
    # path stuff
    tmpdir_str = str(tmpdir)
    tmpdir_repo = os.path.join(tmpdir_str, "single_cell")
    # get from pachterlab
    subprocess.check_call(
        "git clone --no-checkout https://github.com/pachterlab/scRNA-Seq-TCC-prep.git single_cell".split(),
        cwd=tmpdir_str)
    subprocess.check_call("git config core.sparseCheckout true".split(),cwd=tmpdir_repo)
    subprocess.check_call('echo "example_dataset/fastq_files/*ATCGCTCC*" > .git/info/sparse-checkout'.split(),cwd=tmpdir_repo)
    subprocess.check_call("git checkout 0469873bdadcc48e34782882dbd24c3939c0542a".split(),cwd=tmpdir_repo)
    # return location if not tarballed
    fastqs_location = os.path.join(tmpdir_str, "single_cell", "example_dataset", "fastq_files")
    if not tarball:
        return fastqs_location
    # esle, tarball and return that location
    # todo: need this? fastq_files = os.path.listdir(fastqs_location)
    tarball_files(output_dir=tmpdir_str, tar_name='test_fastq.tar.gz', file_paths=[fastqs_location])
    return os.path.join(tmpdir_str, 'test_fastq.tar.gz')



def _generate_config(tmpdir, output_dir, generate_graphs):
    tmpdir_str = str(tmpdir)
    config_location = os.path.join(tmpdir_str, 'config-toil-rnaseqsc-test.yaml')
    if os.path.isfile(config_location):
        os.remove(config_location)
    #todo: the real index location
    with open(config_location, 'w') as f:
        f.write(textwrap.dedent("""
        kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/transcript.idx
        output-dir: {output_dir}
        generate-graphs: {generate_graphs}
        barcode-length: 14
        window-min: 500
        window-max: 5000
        sample-idx: [ATCGCTCC]
        ci-test: true
                """.format(output_dir=output_dir, generate_graphs="true" if generate_graphs else "false")))
    return config_location
