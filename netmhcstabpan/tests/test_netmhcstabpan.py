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

class TestNetmhcstabpan:

    def test_netmhcstabpan(self):
        """
        """
        clean()

        cmd = ['nextflow',
               'test_netmhcstabpan.nf',
               '-entry',
               'test_netmhcstabpan',
               '--pwd',
               os.getcwd(),
               '-config',
               'test_netmhcstabpan.config',
               '--no-resume']
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(),
                                'null',
                                'test_dataset',
                                'test_patient',
                                'test_norm_run_test_tumor_run',
                                'netmhcstabpan',
                                'test.fa.netmhcstabpan.txt')

        eo_md5 = md5()
        with open(expected_output , 'r') as fo:
                filt = '\n'.join([x for x in fo.readlines() if not(re.search("^#", x))])
                eo_md5.update(str(filt).encode('utf-8'))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == 'a8d85d5da757115a8d1d3fcbd421cbc5'

        clean()
