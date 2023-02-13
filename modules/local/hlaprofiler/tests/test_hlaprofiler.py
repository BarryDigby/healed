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


class TestHlaprofilerPredict:
    def test_hlaprofiler_predict(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_hlaprofiler.nf",
            "-entry",
            "test_hlaprofiler_predict",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_hlaprofiler.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_run",
            "hlaprofiler_predict",
            "test_dataset-test_patient-test_run.HLATypes.txt",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str("".join([i[:4] for i in fo.readlines()])).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "a6f2745a48497924094bae08e17ddd91"

        clean()
