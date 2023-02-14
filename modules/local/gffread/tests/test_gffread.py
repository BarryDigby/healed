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


class TestGffreadMakeTxFa:
    def test_gffread_make_tx_fa(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_gffread.nf",
            "-entry",
            "test_gffread_make_tx_fa",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_gffread.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(), "work", "*", "*", "test.test.transcripts.fa"))[0]

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "4b6f1ec4d4294704ba61ebec68ad3370"

        clean()
