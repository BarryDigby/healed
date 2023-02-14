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


class TestNetmhcpan:
    def test_netmhcpan(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_netmhcpan.nf",
            "-entry",
            "test_netmhcpan",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_netmhcpan.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tumor_run",
            "netmhcpan",
            "test.fa.netmhcpan.txt",
        )

        eo_md5 = md5()
        with open(expected_output, "r") as fo:
            filt = "\n".join([x for x in fo.readlines() if not (re.search("Tmp", x))])
            eo_md5.update(str(filt).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "0333230db29c027e843d822e5e952336"

        clean()
