# Output

ScanNeo2 returns its output files in the `results/<name/of/sample>` folder (as specified in the config file). The final output is located within `prioritization`.

```
results/<name/of/sample>
  - dnaseq
    - align
    - qualitycontrol
    - reads
    - indel 
  - rnaseq 
    - align
    - altsplicing
    - indel
    - exitron
    - qualitycontrol
    - reads
  - variants
  - annotation
  - hla
  - prioritization
```

Note: the final results are located in the prioritization folder. In addition, these individual folders contain separate results for each <group> (as specified in the configuration file) but are merged in later stages of ScanNeo2. 

## PRE-PROCESSING



### Quality Control

### Pre-processed reads

## HLA

This folder contains the results of the HLA genotyping. Here, the files `mhc-I.txt` and `mhc-II.txt` contain the alleles for MHC class I and class II, respectively.

e.g., 
```
custom    HLA-A*68:01
custom    HLA-B*15:07
custom    HLA-C*07:04
custom    HLA-A*02:01
custom    HLA-C*03:03
custom    HLA-B*44:02
```

Note: `custom` corresponds to the 


This is the same format as in 



## ALIGNMENT

## VARIANT CALLING

ScanNeo2 calls different variants (according to the configuration) and then collects the results in VCF in the folder `results/<name/of/sample>/variants/` (except for gene fusion events). This can be helpful to gain insights into the variants of each type.

### ALTERNATIVE SPLICING

### EXITRON

Exitron events are called using [ScanExitron](https://github.com/ylab-hi/ScanExitron) and the results are stored in `results/<name/of/sample>/rnaseq/exitron/`. Most importantly, ScanExitron generates the <group>.exitron file which contains all the predicted exitron events. Please consult [ScanExitron](https://github.com/ylab-hi/ScanExitron) for a detailed description of the data fields. In addition, the intermediate results (*.janno) are also kept. ScanNeo2 takes the output of ScanExitron and first converts it into VCF (`<group>_exitron.vcf`). In the next step, this file is augmented with information about the `<group>` and source (`exitron`), which is stored in the keys GRP and SRC of the INFO field, respectively. The file `<group>_exitrons_augmented.vcf` is generated for this. Finally, the files are sorted (`<group>_exitrons.vcf.gz`), and merged into `results/<name/of/sample>/variants/exitrons.vcf.gz`.

### INDEL/SNVs

## PRIORITIZATION

In the prioritization, the output files are generated in `results/<name/of/sample>/prioritization/`. For each variant type, this includes the `<variant_type>_variant_effects.tsv` in which the effects of each variant are listed, and for each MHC class the file `<variant_type>_<mhc_class>_neoepitopes.txt` which contains the detected neoepitopes. In addition, `<mhc_class>_neoepitopes_all.txt` contains all detected neoepitopes in one file. 

The folder structure looks this this (if all modules were activated)
```
﻿- altsplicing_mhc-I_neoepitopes.txt
- altsplicing_variant_effects.tsv
- exitrons_mhc-I_neoepitopes.txt
- exitrons_variant_effects.tsv
- fusions_mhc-I_neoepitopes.txt
- fusions_variant_effects.tsv
- long.indels_mhc-I_neoepitopes.txt
- long.indels_variant_effects.tsv
- somatic.short.indels_mhc-I_neoepitopes.txt
- somatic.short.indels_variant_effects.tsv
- somatic.snvs_mhc-I_neoepitopes.txt
- somatic.snvs_variant_effects.tsv
- custom_protein_mhc-I_neoepitopes.txt
- custom_protein_variant_effects.tsv
- mhc-I_neoepitopes_all.txt

﻿- altsplicing_mhc-II_neoepitopes.txt
- altsplicing_variant_effects.tsv
- exitrons_mhc-II_neoepitopes.txt
- exitrons_variant_effects.tsv
- fusions_mhc-II_neoepitopes.txt
- fusions_variant_effects.tsv
- long.indels_mhc-II_neoepitopes.txt
- long.indels_variant_effects.tsv
- somatic.short.indels_mhc-II_neoepitopes.txt
- somatic.short.indels_variant_effects.tsv
- somatic.snvs_mhc-II_neoepitopes.txt
- somatic.snvs_variant_effects.tsv
- mhc-II_neoepitopes_all.txt                     
```
These include the files `<vartype>_variant_effects.tsv` and include `variant_effects.txt` and individual files for predicted MHC classes (e.g., `mhc-I_neoepitopes.txt` and `mhc-II_neoepitopes.txt`). The former is an intermediate file that contains the variants and their effects on the protein sequence. It can be used as a reference and provides more information about the variants. The following table describes its content.

| Field         | Value     | Description |
|--------------|-----------|------------|
| chrom | String | Chromosome in which the variant occurs. In the case of fusion events, this describes the chromosome of each segment, separated by `\|` (e.g., `chr1\|chr2`) |
| start | Integer | Reference position (0-based) of the variant. In the case of fusion events, this describes the start position of each segment, separated by `\|` (e.g., `1341234\|418728`) |
| end | Integer | End position of the variant |
| gene_id | String | Gene ID the variant occurs in |
| gene_name | String | Corresponding gene name the variant occurs in |
| transcript_id | String | Correspond transcript id in the variant occurs in |
| source | String | The source of the variant (e.g., SNV, Indel,..) | 
| group | String | The group of the variant |
| var_type | String | The effect of the variant (e.g., inframe deletion,...) |
| wt_subseq | String | Wildtype protein sequence (flanking left and right) of the variant |
| mt_subseq | String | Mutant protein sequence (flanking left and right) of the variant |
| var_start | Number | 0-based start position of the variant in the annotation |
| aa_var_start | Integer | Start position (0-based) of the variant within the subsequence |
| aa_var_end | Integer | Start position (0-based) of the variant within the subsequence |
| vaf | Float | Corresponding variant allele frequency (if available) | 
| ao | Float | Observed alleles/reads that support the variant |
| dp | Float | Sequencing depth at the position of the variant |
| TPM | Float | Transcript per Million (TPM) for the corresponding transcript |
| NMD | String | Indicates if the variant is involved in the nonsense-mediated decay (NMD) pathway |
| PTC_dist_ejc | Integer | Distance of the premature stop codon (PTC) to the next exon junction |
| PTC_exon_number | Integer | Exon number the PTC occurs in |
| NMD_escape_rule | Integer | Rule used to escape the NMD pathway (if applicable) |

In addition, the `mhc-I_neoepitopes.txt` is partly redundant to `variant_effects.txt` and contains the following fields:

| Field         | Value     | Description |
|--------------|-----------|------------|
| chrom | String | Chromosome in which the variant occurs. In the case of fusion events, this describes the chromosome of each segment, separated by `\|` (e.g., `chr1\|chr2`) |
| start | Integer | Reference position (0-based) of the variant. In the case of fusion events, this describes the start position of each segment, separated by `\|` (e.g., `1341234\|418728`) |
| end | Integer | End position of the variant |
| allele | String | HLA allele |
| gene_id | String | Gene ID the variant occurs in |
| gene_name | String | Corresponding gene name the variant occurs in |
| transcript_id | String | Correspond transcript id in the variant occurs in |
| source | String | The source of the variant (e.g., SNV, Indel,..) | 
| group | String | The group of the variant |
| var_type | String | The effect of the variant (e.g., inframe deletion,...) |
| var_start | Number | 0-based start position of the variant in the annotation |
| wt_epitope_seq | String | wildtype sequence of the epitope |
| wt_epitope_seq_ic50 | Float | binding affinity of the wildtype sequence |
| wt_epitope_rank | Float | rank of the wildtype epitope | 
| mt_epitope_seq | String | mutant sequence of the epitope |
| mt_epitope_seq_ic50 | Float | binding affinity of the mutant sequence |
| mt_epitope_rank | Float | rank of the mutant epitope | 
| vaf | Float | variant allele frequency |
| supporting | Integer | reads supporting the variant |
| TPM | Float | Transcripts per Million |
| agretopicity | Float | agretopicity score defined mt_affinity/wt_affinity |
| NMD | String | Indicates if the variant is involved in the nonsense-mediated decay (NMD) pathway. Populated for all frameshift variants (SNV / short-indel / long-indel / exitron / alt-splicing / fusion). Values: `NMD_variant`, `NMD_escaping_variant`, or empty (`.`) when no PTC can be determined. |
| PTC_dist_ejc | Integer | Distance of the premature stop codon (PTC) to the next exon junction |
| PTC_exon_number | Integer | Exon number the PTC occurs in |
| NMD_escape_rule | Integer | Rule used to escape the NMD pathway (if applicable) |
| wt_immunogenicity | Float | Immunogenicity score of the wildtype epitope. A higher score indicates a greater probability of eliciting an immune response |
| mt_immunogenicity | Float | Immunogenicity score of the mutant epitope. A higher score indicates a greater probability of eliciting an immune response |
| self-similarity | Float | Similarity measure between the wildtype and mutant epitope. Float values between 0 and 1. `0` Indicates no similarity or a complete difference between the WT and MT sequences. `1` Indicates perfect similarity, meaning the WT and MT sequences are identical in terms of their k-mer similarities.
| pathogen_similarity | Float | Similarity measure between the mutant epitope and known pathogens - more details below |
| pathogen_evalue | Float | BLAST e-value for the pathogen similarity |
| pathogen_bitscore | Float | BLAST bitscore for the pathogen similarity |
| pathogen | String | Name of the detected (similar) pathogen |
| proteome_similarity | Float | Similarity measure between the mutant epitope and the (human) proteome |
| proteome_evalue | Float | BLAST e-value for the proteome similarity |
| proteome_bitscore | Float | BLAST bitscore for the proteome similarity |
| protein | String | Name/ID of the detected (similar) protein |

### pathogen/proteome similarity 

The sequence similarity `ssim` is defined as:
```math
ssim = \frac{\text{identity}}{100}*aligncov
```
where `aligncov` is defined as:
```math
aligncov = \frac{\text{length of alignment}}{\text{length of mutant epitope}}
```





# Logs












