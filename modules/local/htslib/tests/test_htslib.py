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


def clean():
    """ """
    work_dir = pjoin(os.getcwd(), "work")
    globbed_works = glob(pjoin(work_dir, "[A-Za-z0-9][A-Za-z0-9]"))

    for i in globbed_works:
        shutil.rmtree(i, ignore_errors=True)

    shutil.rmtree(pjoin(os.getcwd(), "null"), ignore_errors=True)


class TestHtslibBgzipSomatic:
    def test_htslib_bgzip_somatic(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_htslib.nf",
            "-entry",
            "test_htslib_bgzip_somatic",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_htslib.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tum_run",
            "htslib_bgzip",
            "test.vcf.gz",
        )

        eo_md5 = md5()
        with gzip.open(pjoin(expected_output), "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "2519eb82368afe4ddc1382e43f17baed"

        clean()


class TestHtslibBgzip:
    def test_htslib_bgzip(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_htslib.nf",
            "-entry",
            "test_htslib_bgzip",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_htslib.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(), "null", "test_dataset", "test_patient", "test_run", "htslib_bgzip", "test.vcf.gz"
        )

        eo_md5 = md5()
        with gzip.open(pjoin(expected_output), "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "2519eb82368afe4ddc1382e43f17baed"

        clean()


class TestHtslibTabix:
    def test_htslib_tabix(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_htslib.nf",
            "-entry",
            "test_htslib_tabix",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_htslib.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(), "work", "*", "*", "*tbi"))[0]

        eo_md5 = md5()
        with gzip.open(pjoin(expected_output), "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "cda739f03fa976bb3e171b9cf6a17ca1"

        clean()
