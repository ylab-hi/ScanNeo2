# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.5.x] -2023-xx-xx (unreleased)


### Features

- NMD information (e.g., escape rule,...) is now also calculated for all variants

## [0.3.5] - 2025-10-28

### Fix 

- added 10 retries to curl calls to avoid errors when downloading files (VEP cache, plugins)
- prioritization of MHC-II variants uses now the standalone version of IEDB 
- explicitely chop HLA alleles to first two fields (e.g., HLA-A*02:01:01 becomes HLA-A*02:01)
- changed mhc-II refset to match IEDBs binding affinity prediction
- shell=True caused interactive python session instead of subprocess (when predicting binding affinities MHC-II))

### Changed

- code polishing
- install instructions for HLA-HD in README

## [0.3.4] - 2025-04-04

### Fix 

- Moved script for preparing hla input to workflow/scripts/genotyping
- Added missing pipe to log file in rule hlatyping_mhcII
- Added wrong input file to filtering of mhcII reads on SE

## [0.3.3] - 2025-04-02 

### Fix 

- Fixed parameter in `.tests/integration/config_basic/config.yaml` (align) to match general config

## [0.3.2] - 2025-03-29

### Fix 

- Removed fixed versions for bwa (0.7.19) and samtools (1.21) to allow installation 
of the latest versions.
- Updated the VEP wrapper version from "v1.31.1/bio/vep/annotate" to 
"v5.9.0/bio/vep/annotate".
- Modified functions for file type handling: adjusted indentation in get_preproc_input, 
updated get_input_filter_reads_mhcII_PE to handle BAM files directly, renamed 
get_aligned_reads to get_aligned_reads_featurecounts, added new function 
get_output_hlatyping_mhcII, and updated rule input parameters accordingly.


## [0.3.1] - 2025-03-27

### Fix 

- Changed protocol for HLA alleles reference list to https (rule get_hla_info)
- Fixed path to input files in finalize_mhcII_input.py (which caused error when using paired-end reads)
- Fixed bug in quantification/featurecounts - supported regardless of input type (PE or SE): added wrapper script
- replaced MHC-II binding affinity prediction with IEDB API
- Manual Dockerfile replaces containerized version

## [0.3.0] - 2024-08-30

### Features

- Added sequence similarity filter for MHC-I
    - self-similarity (using kernel similarity)
    - pathogen similarity (BLAST against pathogen-derived epitopes from IEDB)
    - proteome similarity (BLAST against human proteome)
- Prioritization of neoantigens is now done separately for each variant type (speeds up the process)
- Update to recent version of ScanExitron
    - this version updated to recent version of regtools (v0.5.0) - which is available on Conda
    - Singularity/Docker is not necessary anymore
    - Added option to use strand information in exitron calling
- ScanNeo2 now uses conda environments for all tools (ditched Singularity/Docker)

### Fix 
- renamed similarity fields for pathogen and protein to more descriptive names 

## [0.2.12] - 2024-08-21 

## Fix 

- Fixed missing alleles in HLA alleles reference list - [#34](https://github.com/ylab-hi/ScanNeo2/issues/34)

## [0.2.11] - 2024-08-02

## Fix 

- Updated transindel environment to recent samtools version (as --o introduced in samtools >= 1.13 required by transindel)

## [0.2.10] - 2024-07-08 

### Fix 

- Allow to combine multiple VCF files in indel detection using mutect2 (e.g., when multiple samples are provided)

## [0.2.9] - 2024-07-04

### Fix 

- Splitted rules in HLA typing to ensure better distribution of the workload
- Changed order in HLA typing rules (BAM files are now part of single-end)
    - samtools fastq is only called for BAM files
    - input of filtering directly from preprocessed/raw reads

## [0.2.8] - 2024-06-26

### Fix 

- Added threads option to samtools sort calls to speed up the process
- Fixed wrong call to optitype within the wrapper script

## [0.2.7] - 2024-06-23 

### Fix 

- Separated samtools, bcftools and realign environments to avoid conflicts
- Changed order of genotyping rules to catch errors when no alleles can be found
    - Alleles are merged according to nartype (e.g., DNA, RNA) and then combined
- Force concat of VCF files in genotyping to avoid errors when no variants are found
- Added optitype wrapper to avoid errors when empty BAM files are provided / no HLA reads

## [0.2.6] - 2024-06-20

### Fix 

- Added routines to catch errors when rnaseq data is not provided but exitron/alternative splicing calling is activated
- Added reference genome index as input to germline indel calling (necessary when only indel calling is activated)
- removed -C from BWA mem call (on DNAseq data) to avoid error on Illumina identifiers

## [0.2.5] - 2024-06-19

### Fix 

- Wrong indentation in HLAtyping caused error when providing no normal sample (NoneType was being iterated)
- Fixed missing input in get_reads_hlatyping_PE rule (tmp folder) that caused error when using paired-end reads
- Added else case in get_input_hlatyping function (rule get_reads_hlatyping_PE) for input reads when preproc is deactivated

## [0.2.4] - 2024-05-19

### Fix

- Added concurrency to splAdder call
- Added routines that lets ScanNeo2 finish (even when splAdder results are empty or faulty)

## [0.2.3] - 2024-03-02

### Features (somewhat)

- Added paramter nonchr in reference attribute to exclude non-chromosomal contigs from the reference genome

### Fix

- Fixed wrong path in quality control for single-end reads

## [0.2.2] - 2024-03-01

### Fix 

- Conda instal wheel caused error on the spladder environment
    - pysam requires exactly python=3.6

## [0.2.1] - 2024-02-29

### Fix

- Removed hlahd path from config and hlatyping - needs to be installed in $PATH


## [0.2.0] - 2024-02-25

### Features

- ScanNeo2 supports Snakmake>=8 
    - --use-conda replaced by --software-deployment-method conda
    - --use-singularity replaced by --software-deployment-method apptainer
- Gather/scatter of the indel calling speeds up ScanNeo2 on multiple cores
    - added script to split bamfiles by chromosome (scripts/split_bam_by_chr.py)
    - haplotypecaller first/final round is done per chromosome and later merged
    - mutect2 is done per chromosome and later merged
- Genotyping MHC-II works now on both single-end and paired-end
- User-defined HLA alleles are matched against the hla refset
- Added multiple routine to catch errors when only custom variants are provided
- Added additional parameters in config file

### Fix 

- When using BAMfiles the HLA typing wrongly expected single-end reads and performed preprocessing
- Each environment is no thoroughly versioned to ensure interoperability
- Missing immunogenicity calculation on certain values of MHC-I fixed
- Fixed prediction of binding affinity in MHC-II (as the columns are different from MHC-I)


## [0.1.6] - 2024-02-13

### Fix 

- linked rules for prediction of binding affinities and immunogenicity to input of prioritization


## [0.1.5] - 2024-02-13

### Fix 

- fixed wrong reference genome in exitron2vcf call (which forced ScanNeo2 to use alternative rules)
- removed redundant rules for alternative genomes (ScanNeo2 now uses ensembl globally)

## [0.1.4] - 2024-02-12

### Added

- added input directive in rule `prepare_cds` (exitron rules). Makes sure that annotations are present if exitron calling is executed first

## [0.1.3] - 2024-02-10

### Added

- added routines/fixed issues when no normal sample is provided

## [0.1.2] - 2024-01-17

### Added

- Added alternative link for VEP cache to improve download speed
- Added missing scripts to modify the ensembl header
- Modularized rule for long indel detection

## [0.1.1] - 2024-01-16

### Added

- Fixed errors when providing custom input for MHC alleles
- Refactoring of genotyping scripts 
- Added more detailed instructions in README

## [0.1.0] - 2023-08-17

### Added

- Comprehensive workflow with different modules to detect variants from sequencing data
- Different modules for each step
- Support data in single-end, paired-end .fastq or BAM files
- preprocessing, alignment
- genotyping
- alternative splicing
- gene fusion, exitron, SNVs and indels
