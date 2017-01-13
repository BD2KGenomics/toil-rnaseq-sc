import logging
import shlex
import shutil
import subprocess
import tempfile
import textwrap
from contextlib import closing
from unittest import TestCase
from urlparse import urlparse
from uuid import uuid4

import os
import posixpath
from bd2k.util.iterables import concat
from boto.s3.connection import S3Connection, Bucket

log = logging.getLogger(__name__)


class RNASeqSCCGLTest(TestCase):
    """
    These tests *can* be parameterized with the following optional environment variables:

    TOIL_SCRIPTS_TEST_TOIL_OPTIONS - a space-separated list of additional command line arguments to pass to Toil via
    the script entry point. Default is the empty string.

    TOIL_SCRIPTS_TEST_JOBSTORE - the job store locator to use for the tests. The default is a file: locator pointing
    at a local temporary directory.

    TOIL_SCRIPTS_TEST_NUM_SAMPLES - the number of sample lines to generate in the input manifest
    """

    @classmethod
    def setUpClass(cls):
        super(RNASeqSCCGLTest, cls).setUpClass()
        # FIXME: pull up into common base class
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        log.info("In Setup")

    def test_the_test(self):
        log.info("Looks like the test worked")

    def tearDown(self):
        log.info("In TearDown")
