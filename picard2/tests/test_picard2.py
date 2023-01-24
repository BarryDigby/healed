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
import pysam


def clean():
    """
    """
    work_dir = pjoin(os.getcwd(), 'work')
    globbed_works = glob(pjoin(work_dir, '[A-Za-z0-9][A-Za-z0-9]'))

    for i in globbed_works:
        shutil.rmtree(i, ignore_errors=True)

    shutil.rmtree(pjoin(os.getcwd(), 'null'), ignore_errors=True)

class TestPicardMarkDuplicates:

    def test_picard_mark_duplicates(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_picard2.nf',
               '-entry',
               'test_picard_mark_duplicates',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_picard2.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                '*.mkdup.bam'))[0]

        print(expected_output)


        eo_md5 = md5()
        rows = pysam.view("-S", expected_output)
        filt = '\n'.join([x for x in rows if not(re.search("TMP_DIR", x))])
        eo_md5.update(str(filt).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '983032af86194ccfc2368011b5dbe77a'

        clean()

class TestPicardCreateSeqDict:

    def test_picard_create_seq_dict(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_picard2.nf',
               '-entry',
               'test_picard_create_seq_dict',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_picard2.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(),
                                'work',
                                '*',
                                '*',
                                'chr1.dict'))[0]

        print(expected_output)
        eo_md5 = md5()
        with open(pjoin(expected_output) , 'r') as fo:
            filt = [re.sub(r'UR:.*$', '', x) for x in fo.readlines()]
            eo_md5.update(str(''.join(filt)).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == '55fabdc886d3f04f15a7c9ded2260bf8'
        clean()
