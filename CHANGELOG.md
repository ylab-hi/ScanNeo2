# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
