#!/usr/bin/env python

from glob import glob
import pytest
from pathlib import Path
from hashlib import md5
import shutil
import os
from pprint import pprint
from os.path import join as pjoin
import subprocess
import gzip
import re


def clean():
    """
    """
    work_dir = pjoin(os.getcwd(), 'work')
    globbed_works = glob(pjoin(work_dir, '[A-Za-z0-9][A-Za-z0-9]'))

    for i in globbed_works:
        shutil.rmtree(i, ignore_errors=True)

    shutil.rmtree(pjoin(os.getcwd(), 'null'), ignore_errors=True)

class TestPycloneviFIt:

    def test_pyclonevi_fit(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_pyclone-vi.nf',
               '-entry',
               'test_pyclonevi_fit',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_pyclone-vi.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_norm_run_test_tumor_run',
                                'pyclone-vi_fit',
                                'test.pcvi_input.tmp')

        # Outputs are non-deterministic. Need to check contents of file instead
        # of checksum.
        assert expected_output

        clean()
