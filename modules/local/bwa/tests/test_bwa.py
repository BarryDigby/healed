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

class TestBwaIndex:

    def test_bwa_index(self):
        """
        This passes, but only if you run it twice without cleaning. The issue
        seems to be that the first pass fails (since Nextflow cannot emit the
        original fasta), but the output files _are_ placed in the correct
        output directory. The second Nextflow instance finishes (as it uses the
        storeDir correctly). 
        """
        clean()

        cmd = ['nextflow',
               'test_bwa.nf',
               '-entry',
               'test_bwa_index',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_bwa.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'NC_045512.2.fa',
                                'bwa_index',
                                'NC_045512.2.fa.bwt')

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '40eb7fc45471984999df09a9366a7995'

        clean()

class TestBwaMemSamtoolsSort:

    def test_bwa_mem_samtools_sort(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_bwa.nf',
               '-entry',
               'test_bwa_mem_samtools_sort',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_bwa.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(),
                                    'work',
                                    '*',
                                    '*',
                                    '*bam'))[0]


        eo_md5 = md5()
        with open(pjoin(expected_output) , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '5f285d8a2cfc4faa6d14d80b1103c819'

        clean()
