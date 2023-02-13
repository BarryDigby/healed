#! /usr/bin/env python

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


class TestAbra2:
    def test_abra2(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_abra2.nf",
            "-entry",
            "test_abra2",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_abra2.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_norm_output = glob(pjoin(os.getcwd(), "work", "*", "*", "*norm_abra.bam"))[0]

        print(expected_norm_output)

        eno_md5 = md5()
        with open(expected_norm_output, "rb") as fo:
            eno_md5.update(str(fo.read()).encode("utf-8"))
        eno_md5 = eno_md5.hexdigest()

        assert eno_md5 == "6e1fe0e1408102f6fc26edb616661165"

        expected_tum_output = glob(pjoin(os.getcwd(), "work", "*", "*", "*tumor_abra.bam"))[0]

        print(expected_tum_output)

        eto_md5 = md5()
        with open(expected_tum_output, "rb") as fo:
            eto_md5.update(str(fo.read()).encode("utf-8"))
        eto_md5 = eto_md5.hexdigest()

        assert eto_md5 == "c5e783307b63acbf5fed866961b197d3"
        clean()


class TestAbra2Cadabra:
    def test_abra2_cadabra(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_abra2.nf",
            "-entry",
            "test_abra2_cadabra",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_abra2.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_tumor_run_test_norm_run",
            "abra2_cadabra",
            "test_dataset-test_patient-test_tumor_run_test_norm_run.abra2.vcf",
        )

        print(expected_output)

        eo_md5 = md5()
        with open(expected_output, "r") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "002be2835174c3436a8d54efa3a70816"
        clean()
