# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2024-03-07

### Features

- Added sequence similarity filter for MHC-I
    - self-similarity (using kernel similarity)
    - pathogen similarity (BLAST against pathogen-derived epitopes from IEDB)
    - proteome similarity (BLAST against human proteome)
- Prioritization of neoantigens is now done separately for each variant type (speeds up the process)
- NMD information (e.g., escape rule,...) is now also calculated for all variants

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
