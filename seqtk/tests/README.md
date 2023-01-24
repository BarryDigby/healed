Testing is not currently set up for CI.

Note: This module is not producing consistent files. These tests may be helpful
to see expected outputs, but will fail the checksum test.

Testing can be performed locally, but requires:
- Nextflow
- Singularity
- Python3 and PyTest

Testing workflow:
git clone https://gitlab.com/bgv-lens/nextflow/modules/tools/seqtk.git
cd seqtk/tests
python3 -m pytest
