# Integration Tests

This directory contains integrations tests for ScanNeo2. 
Here, we provide a set of configuration files and expected 
output files that can be used to test the installation of 
ScanNeo2. 

Here, we provide folders for a set of different test cases.

In `custom-test` we call ScanNeo2 using a set of pre-calculated
SNVs (`data/variants/snvs.vcf`) using provided mhc-I alleles
(`data/hla/mhc-I.tsv`). In addition, the expected output is
stored in `results`.

In `indel-test` we call ScanNeo2 to detect indel-derived neoantigens
using provided mhc-I alleles (`data/hla/mhc-I.tsv`). The exptected
output is stored in `results`.

## Running the Tests

```
snakemake --cores all --configfile .tests/integration/<testcase>/config/config.yml

```
The expected output is stored in the `results` folder in the root directory of ScanNeo2.
