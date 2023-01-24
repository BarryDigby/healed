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

class TestSamtoolsFaidx:

    def test_samtools_faidx(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_faidx',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = glob(pjoin(os.getcwd(),
                                    'work',
                                    '*',
                                    '*',
                                    '*fai'))[0]

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '214bc367026d784167273beee5e3e2d9'

        clean()

class TestSamtoolsCoverage:

    def test_samtools_coverage(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_coverage',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'samtools_coverage',
                                'test.bam.coverage')

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '12f2ba7ab8b7c54f34f86b78c05ca4a6'

        clean()

class TestSamtoolsIndex:

    def test_samtools_index(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_index',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                '*bai'))[0]

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '923ea4e9ea87c72ce94cb311a08ebcde'

        clean()

class TestSamtoolsRmdup:

    def test_samtools_rmdup(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_rmdup',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                '*rmdup.bam'))[0]

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '09009449780846a7100cf6005ca0e541'

        clean()

class TestSamtoolsSort:

    def test_samtools_sort(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_sort',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                '*sorted.bam'))[0]

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == 'a9f15e9707964445f4dd6b8a24ddda2e'

        clean()

class TestSamtoolsStats:

    def test_samtools_stats(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_stats',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_run',
                                'samtools_stats',
                                'test.bam.stats')

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'r') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '9e9d4229c9f7bcb64979a0607c26a147'

        clean()

class TestSamtoolsView:

    def test_samtools_view(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_samtools.nf',
               '-entry',
               'test_samtools_view',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_samtools.config',
               '--no-resume']
        subprocess.run(cmd)
        
        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                '*bam.bam'))[0]

        print(expected_output)


        eo_md5 = md5()
        with open(expected_output , 'rb') as fo:
                eo_md5.update(str(fo.read()).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '0c930b2ca8a63dc1f680e63d6b7e2853'

        clean()
