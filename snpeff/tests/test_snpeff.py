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

class TestSnpeffAnn:

    def test_snpeff_ann(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_snpeff.nf',
               '-entry',
               'test_snpeff_ann',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_snpeff.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'norm_test_run_tumor_test_run',
                                'snpeff_ann',
                                'test.annot.vcf')


        eo_md5 = md5()
        with open(pjoin(expected_output) , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == 'da86cb8190f9e94b8f17a27c2ee13f2f'

        clean()

class TestSnpeffAnnGermline:

    def test_snpeff_ann_germline(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_snpeff.nf',
               '-entry',
               'test_snpeff_ann_germline',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_snpeff.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'snpeff_ann',
                                'test.annot.vcf')


        eo_md5 = md5()
        with open(pjoin(expected_output) , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == 'da86cb8190f9e94b8f17a27c2ee13f2f'

        clean()
