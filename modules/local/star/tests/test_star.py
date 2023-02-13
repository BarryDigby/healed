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


class TestStarIndex:
    def test_star_index(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_star.nf",
            "-entry",
            "test_star_index",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_star.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(os.getcwd(), "null", "NC_045512.2.fa", "star_index", "index_files", "SA")

        print(expected_output)

        eo_md5 = md5()
        with open(expected_output, "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "52938577dcfc4654772a0a1c22484eb5"

        clean()


class TestStarMap:
    def test_star_map(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_star.nf",
            "-entry",
            "test_star_map",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_star.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = glob(pjoin(os.getcwd(), "work", "*", "*", "*bam"))[0]

        eo_md5 = md5()
        with open(pjoin(expected_output), "rb") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "acdadd94a3ef87054f6692811ddfcc53"

        clean()
