import logging
import os
import textwrap
import subprocess
from toil_lib import require, UserError
from toil_lib.files import tarball_files
from toil_rnaseq_sc.rnaseq_sc_cgl_pipeline import parse_samples

log = logging.getLogger(__name__)


def test_parse_samples(tmpdir):

    ### verifying parse_samples succeed as expected ###

    # specific files - all file types and schemes
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
    require(len(url2) == 1, "expected to have one file specified in url2")

    # s3 zipped tarball
    in3 = [["test3", "file:///example/tarball.tar.gz"]]
    _generate_manifest(tmpdir, in3)
    out3 = parse_samples(manifest_location)
    require(len(out3) == 1, "expected to have single output for test3")
    uuid3, url3 = out3[0]
    require(uuid3 == "test3", "expected uuid to be 'test3'")
    require(len(url3) == 1, "expected to have one file specified in url3")

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

    ### verifying parse_samples fail as expected ###
    # only uuid specified
    _generate_manifest(tmpdir, [["only_uuid"]])
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, "parse_samples should have failed due to malformed manifest")

    # multiple samples specified
    _generate_manifest(tmpdir, [["uuid", "file:///file1.fq", "file:///file2.fq"]])
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, "parse_samples should have failed due to malformed manifest")

    # tests bad scheme fails
    _generate_manifest(tmpdir, [["uuid", "badscheme:///file.fq"]])
    try:
        parse_samples(manifest_location)
    except UserError:
        pass
    else:
        require(False, "parse_samples should have failed due to bad scheme")


def test_pipeline_output_with_graphs(tmpdir):
    uuid = "test_rnaseqsc_g"
    output_dir = os.path.join(str(tmpdir), uuid)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, uuid + ".tar.gz")
    jobstore = os.path.join(str(tmpdir), uuid + "_jobstore")

    input = "file://" + _get_test_fastq_files(tmpdir, tarball=True)
    config = _generate_config(tmpdir, output_dir, generate_graphs=True)
    manifest = _generate_manifest(tmpdir, [[uuid, input]])

    subprocess.check_call(['toil-rnaseq-sc', 'run',
                                '--config', config,
                                '--manifest', manifest,
                                '--maxCores', "1",
                                jobstore])
    # ensure file and directories exist
    require(os.path.isfile(output_file), "expected outputfile to exist: " + output_file)
    subprocess.check_call(['tar', '-xvf', output_file, '-C', output_dir])
    require(os.path.isfile(os.path.join(output_dir, uuid, "matrix.tsv")),
            "matrix.tsv file should exist in output tarball")
    require(os.path.isdir(os.path.join(output_dir, uuid, "plots")),
            "plots directory should exist in output tarball")
    require(len(os.listdir(os.path.join(output_dir, uuid, "plots"))) > 0,
            "plots directory should not be empty in output tarball")


def test_pipeline_output_without_graphs(tmpdir):
    uuid = "test_rnaseqsc_ng"
    output_dir = os.path.join(str(tmpdir), uuid)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, uuid + ".tar.gz")

    jobstore = os.path.join(str(tmpdir), uuid + "_jobstore")

    input = _get_test_fastq_files(tmpdir, tarball=False)
    config = _generate_config(tmpdir, output_dir, generate_graphs=False)
    manifest = _generate_manifest(tmpdir, [[uuid, input]])

    subprocess.check_call(['toil-rnaseq-sc', 'run',
                                '--config', config,
                                '--manifest', manifest,
                                '--maxCores', "1",
                                jobstore])
    # ensure file and directories exist (or don't)
    require(os.path.isfile(output_file), "expected outputfile to exist: " + output_file)
    subprocess.check_call(['tar', '-xvf', output_file, '-C', output_dir])
    require(os.path.isfile(os.path.join(output_dir, uuid, "matrix.tsv")),
            "matrix.tsv file should exist in output tarball")
    require(not os.path.isdir(os.path.join(output_dir, uuid, "plots")),
            "plots directory should not exist in output tarball")


def _generate_manifest(tmpdir, list_of_lines):
    manifest_location = os.path.join(str(tmpdir), "manifest-toil-rnaseqsc-test.tsv")
    # clear old manifest if it exists
    if os.path.isfile(manifest_location):
        os.remove(manifest_location)
    with open(manifest_location, 'w') as manifest:
        for line in list_of_lines:
            manifest.write("\t".join(line) + "\n")
    return manifest_location

def _generate_config(tmpdir, output_dir, generate_graphs):
    tmpdir_str = str(tmpdir)
    config_location = os.path.join(tmpdir_str, 'config-toil-rnaseqsc-test.yaml')
    if os.path.isfile(config_location):
        os.remove(config_location)
    with open(config_location, 'w') as f: #todo get rid of my kallisto idx!!
        f.write(textwrap.dedent("""
        kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx
        output-dir: {output_dir}
        generate-graphs: {generate_graphs}
        barcode-length: 14
        window-min: 500
        window-max: 5000
        sample-idx: [ATCGCTCC]
                """.format(output_dir=output_dir, generate_graphs="true" if generate_graphs else "false")))
    return config_location



def _get_test_fastq_files(tmpdir, tarball=True):
    # path stuff
    tmpdir_str = str(tmpdir)
    tmpdir_repo = os.path.join(tmpdir_str, "single_cell")
    # get from pachterlab
    subprocess.check_call(
        "git clone --no-checkout https://github.com/pachterlab/scRNA-Seq-TCC-prep.git single_cell".split(),
        cwd=tmpdir_str)
    subprocess.check_call("git config core.sparseCheckout true".split(),cwd=tmpdir_repo)
    with open(os.path.join(tmpdir_repo, ".git", "info", "sparse-checkout"), 'w') as sparse:
        sparse.write("example_dataset/fastq_files/*ATCGCTCC*")
    # subprocess.check_call('echo "example_dataset/fastq_files/*ATCGCTCC*" > .git/info/sparse-checkout'.split(),cwd=tmpdir_repo)
    subprocess.check_call("git checkout 0469873bdadcc48e34782882dbd24c3939c0542a".split(),cwd=tmpdir_repo)
    # return location if not tarballed
    fastqs_location = os.path.join(tmpdir_str, "single_cell", "example_dataset", "fastq_files")
    if not tarball:
        return fastqs_location
    # else, tarball and return that location
    tarball_files(output_dir=tmpdir_str, tar_name='test_fastq.tar.gz', file_paths=
            [os.path.join(fastqs_location, x) for x in os.listdir(fastqs_location)])
    return os.path.join(tmpdir_str, 'test_fastq.tar.gz')

