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

class TestDeepVariant:

    def test_deepvariant(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_deepvariant.nf',
               '-entry',
               'test_deepvariant',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_deepvariant.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'deepvariant',
                                'test_dataset-test_patient-test_run.deepv.germline.vcf')

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 in ('fc8ef361518be0cac3913fe63538012c', '130f56666b06ba53fb4aca3db96c6356', '260cfe5d3d26523381b46fa049fae8be')

        clean()
