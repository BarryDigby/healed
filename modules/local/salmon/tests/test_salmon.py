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

class TestSalmonIndex:

    def test_salmon_index(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_salmon.nf',
               '-entry',
               'test_salmon_index',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_salmon.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test.tx.fa',
                                'salmon_index',
                                'test.tx.fa.index')


        # Need to do checksums on individual files, likely.
        assert expected_output

        # Should not clean here. These reference files are needed for the
        # test_salmon_map_quant test.
#        clean()

class TestSalmonMapQuant:

    def test_salmon_map_quant(self):
        """
        """

        cmd = ['nextflow',
               'test_salmon.nf',
               '-entry',
               'test_salmon_map_quant',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_salmon.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'salmon_map_quant',
                                'test_dataset-test_patient-test_run.quant.sf')

        eo_hdr = ''
        with open(pjoin(expected_output) , 'r') as fo:
            eo_hdr = list(fo.readlines())[0]
            print(eo_hdr)
        assert eo_hdr == 'Name\tLength\tEffectiveLength\tTPM\tNumReads\n'


        clean()

class TestSalmonAlnQuant:

    def test_salmon_aln_quant(self):
        """
        """

        cmd = ['nextflow',
               'test_salmon.nf',
               '-entry',
               'test_salmon_aln_quant',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_salmon.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'salmon_aln_quant',
                                'test_dataset-test_patient-test_run.quant.sf')

        eo_hdr = ''
        with open(pjoin(expected_output) , 'r') as fo:
            eo_hdr = list(fo.readlines())[0]
            print(eo_hdr)
        assert eo_hdr == 'Name\tLength\tEffectiveLength\tTPM\tNumReads\n'


        clean()
