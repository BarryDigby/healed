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
    """ """
    work_dir = pjoin(os.getcwd(), "work")
    globbed_works = glob(pjoin(work_dir, "[A-Za-z0-9][A-Za-z0-9]"))

    for i in globbed_works:
        shutil.rmtree(i, ignore_errors=True)

    shutil.rmtree(pjoin(os.getcwd(), "null"), ignore_errors=True)


class TestSeqtkSample:
    def test_seqtk_sample(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_seqtk.nf",
            "-entry",
            "test_seqtk_sample",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_seqtk.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_run",
            "seqtk_sample",
            "test_1.subd.subd.100K.fastq.gz",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 in [
            "5ba8e7739286216815d11bde9cd2f0bd",
            "1616037655b1028a580fec514e27cb3c",
            "4d2f51bd76e8b09383cafbd369d4dcbf",
            "05f00a4b78dc7a1618eb0de96b1b533b",
        ]

        clean()
