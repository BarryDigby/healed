Testing is not currently set up for CI.

Testing can be performed locally, but requires:

- Nextflow
- Singularity
- Python3 and PyTest

Testing workflow:
git clone https://gitlab.com/bgv-lens/nextflow/modules/tools/hlaprofiler.git
cd hlaprofiler/tests
python3 -m pytest
