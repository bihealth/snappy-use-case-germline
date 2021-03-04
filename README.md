# SNAPPY Use Case "Germline"

This repository contains the
[SNAPPY](https://github.com/bihealth/snappy-pipeline)
use case for germline data processing. Processing starts from FASTQ files, and the result are
files that can be loaded into
[VarFish](https://github.com/bihealth/varfish-server)

**Setup**

```shell
# Preparation
git clone --recurse-submodules git@github.com:bihealth/snappy-use-case-germline.git
cd snappy-use-case-germline
pip install -r requirements.txt

# Setup
python code/set_use_case.py


```


## Quick Info

- License: MIT
- Copyright: 2018-2021 by Core Unit Bioinformatics (CUBI), Berlin Institute of Health
