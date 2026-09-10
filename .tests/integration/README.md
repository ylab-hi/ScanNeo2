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

In `replicate-test` we call ScanNeo2 with two samples, one of which
carries two RNA-seq tumor replicates (`tumor_rep1`, `tumor_rep2`).
It exercises the long-format sample sheet: the per-group alignment and
variant calling fan out (once per replicate) while prioritization pools
each sample's groups into a single combined, provenance-tagged result.

## Running the Tests

```
snakemake --cores all --configfile .tests/integration/<testcase>/config/config.yml

```
The expected output is stored in the `results` folder in the root directory of ScanNeo2.
