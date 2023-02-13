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
import gzip


def clean():
    """ """
    work_dir = pjoin(os.getcwd(), "work")
    globbed_works = glob(pjoin(work_dir, "[A-Za-z0-9][A-Za-z0-9]"))

    for i in globbed_works:
        shutil.rmtree(i, ignore_errors=True)

    shutil.rmtree(pjoin(os.getcwd(), "null"), ignore_errors=True)


class TestTrimGalore:
    def test_trim_galore(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_trim_galore.nf",
            "-entry",
            "test_trim_galore",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_trim_galore.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_run",
            "trim_galore",
            "test_dataset-test_patient-test_run_1.trimmed.fq.gz",
        )

        print(expected_output)

        eo_md5 = md5()
        with gzip.open(expected_output, "r") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "037f151a2ecb4ab961943963e711523f"

        clean()
