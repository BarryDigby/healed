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


class TestBcftoolsSort:
    def test_bcftools_sort(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_sort",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(), "null", "test_dataset", "test_patient", "test_run", "bcftools_sort", "test_1.sorted.vcf"
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r", encoding="utf-8") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "27922af2f7df690e121c1747cabde2e3"

        clean()

    def test_bcftools_sort_gzipped(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_sort_gzipped",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(), "null", "test_dataset", "test_patient", "test_run", "bcftools_sort", "test_2.sorted.vcf"
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r", encoding="utf-8") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "27922af2f7df690e121c1747cabde2e3"

        clean()


class TestBcftoolsConsensus:
    def test_bcftools_consensus_snv(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_consensus_snv",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tum_run",
            "bcftools_consensus",
            "test_dataset-test_patient-test_norm_run_test_tum_run.test.consensus.fa",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r", encoding="utf-8") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "6d38d65c236250eab32988234c158f94"

        clean()


class TestBcftoolsFilter:
    def test_bcftools_filter(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_filter",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tum_run",
            "bcftools_filter",
            "test.bfilt.vcf",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r", encoding="utf-8") as fo:
            eo_md5.update("".join([i for i in fo.readlines() if not (re.search("Date", i))]).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "9924c034d52fba438318eb1a97348c96"

        clean()


class TestBcftoolsIndex:
    @pytest.mark.skip(reason="Needs to read binary index.")
    def test_bcftools_index(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_index",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(), "null", "test_dataset", "test_patient", "test_run", "bcftools_index", "test.vcf.gz.csi"
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str(fo.read()))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "9dc652dfb036c77af60c3929792c1084"

        clean()

    @pytest.mark.skip(reason="Needs to read binary index.")
    def test_bcftools_index(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_index_somatic",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tum_run",
            "bcftools_index",
            "test.vcf.gz.csi",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str(fo.read()))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "9dc652dfb036c77af60c3929792c1084"

        clean()


class TestBcftoolsStats:
    def test_bcftools_stats(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_stats",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(), "null", "test_dataset", "test_patient", "test_run", "bcftools_stats", "test.vcf_stats"
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "aa3c701c568eaf402d6002c097398e44"

        clean()

    def test_bcftools_stats_somatic(self):
        """ """
        clean()

        cmd = [
            "nextflow",
            "test_bcftools.nf",
            "-entry",
            "test_bcftools_stats_somatic",
            "--pwd",
            os.getcwd(),
            "-config",
            "test_bcftools.config",
            "--no-resume",
        ]
        subprocess.run(cmd)

        expected_output = pjoin(
            os.getcwd(),
            "null",
            "test_dataset",
            "test_patient",
            "test_norm_run_test_tum_run",
            "bcftools_stats",
            "test.vcf_stats",
        )

        eo_md5 = md5()
        with open(pjoin(expected_output), "r") as fo:
            eo_md5.update(str(fo.read()).encode("utf-8"))
        eo_md5 = eo_md5.hexdigest()

        assert eo_md5 == "aa3c701c568eaf402d6002c097398e44"

        clean()
