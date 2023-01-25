Testing is not currently set up for CI.

Note: Need way to access SnpEff annotation information. Annotation can be
generated following directions at
https://gitlab.com/bgv-lens/nextflow/modules/tools/lens/-/wikis/Running-LENS#populate-raft-global-references-directory.

Testing can be performed locally, but requires:
- Nextflow
- Singularity
- Python3 and PyTest

Testing workflow:
git clone https://gitlab.com/bgv-lens/nextflow/modules/tools/snpeff.git
cd snpeff/tests
python3 -m pytest
