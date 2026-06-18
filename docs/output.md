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

Two pre-processing stages run per sample, both gated on `preproc.activate: true`:

### Quality Control

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports land in `results/<sample>/{dnaseq,rnaseq}/qualitycontrol/` per read group, with separate forward / reverse reports for paired-end inputs. Use these to spot adapter contamination, GC bias, and per-base quality drop-off before trimming.

### Pre-processed reads

[fastp](https://github.com/OpenGene/fastp) writes the trimmed reads to `results/<sample>/{dnaseq,rnaseq}/reads/` as `<group>_preproc.fq.gz` (single-end) or `<group>_preproc_r1.fq.gz` / `<group>_preproc_r2.fq.gz` (paired-end). Sliding-window trimming and minimum-length filtering follow the `preproc.slidingwindow` and `preproc.minlen` settings in the config.

## HLA

This folder contains the results of the HLA genotyping. The files `mhc-I.tsv` and `mhc-II.tsv` carry the typed alleles for MHC class I and class II respectively. Each row is `<source>\t<allele>`:

```
DNA       HLA-A*02:01
RNA       HLA-A*02:01
custom    HLA-A*68:01
custom    HLA-B*15:07
```

The first column records where the allele came from:

- `DNA` — predicted from DNA-seq reads (OptiType for class I, HLA-HD for class II).
- `RNA` — predicted from RNA-seq reads.
- `custom` — user-supplied via `data.custom.hlatyping.MHC-{I,II}` in the config; useful when alleles are already known or when running on a sample where read-based typing isn't appropriate.

Multiple sources for the same allele are kept (e.g. an allele typed independently from both DNA and RNA appears twice). Downstream binding-affinity prediction operates on the deduplicated allele set.

## ALIGNMENT

Aligned BAMs land in `results/<sample>/dnaseq/align/` and `results/<sample>/rnaseq/align/`. The two paths use different aligners by design:

- **DNA-seq** is aligned directly with **BWA-MEM** — sufficient for variant calling against a reference genome.
- **RNA-seq** is first aligned with **STAR** in chimeric-aware mode (the `align.chim*` config parameters control chimeric-segment thresholds; see the STAR manual). The STAR BAM is then re-aligned with BWA via the `realign` rule because the downstream RNA-variant callers (transIndel, ScanExitron) need a BWA-style CIGAR. Both intermediates and the final BAM are kept.

`postproc_bam_index` writes `.bai` index files alongside each BAM.

## VARIANT CALLING

ScanNeo2 calls different variants (according to the configuration) and then collects the results in VCF in the folder `results/<name/of/sample>/variants/` (except for gene fusion events). This can be helpful to gain insights into the variants of each type.

### ALTERNATIVE SPLICING

[SplAdder](https://spladder.readthedocs.io/) detects alternative-splicing events from the RNA-seq alignment; intermediate splice-graph files land in `results/<sample>/rnaseq/altsplicing/`. The `splicing_to_vcf` rule converts SplAdder's output to per-group VCFs that are then sorted, augmented with `GRP` / `SRC` INFO keys (same convention as the exitron path below), and merged into `results/<sample>/variants/altsplicing.vcf.gz`.

### EXITRON

Exitron events are called using [ScanExitron](https://github.com/ylab-hi/ScanExitron) and the results are stored in `results/<name/of/sample>/rnaseq/exitron/`. Most importantly, ScanExitron generates the <group>.exitron file which contains all the predicted exitron events. Please consult [ScanExitron](https://github.com/ylab-hi/ScanExitron) for a detailed description of the data fields. In addition, the intermediate results (*.janno) are also kept. ScanNeo2 takes the output of ScanExitron and first converts it into VCF (`<group>_exitron.vcf`). In the next step, this file is augmented with information about the `<group>` and source (`exitron`), which is stored in the keys GRP and SRC of the INFO field, respectively. The file `<group>_exitrons_augmented.vcf` is generated for this. Finally, the files are sorted (`<group>_exitrons.vcf.gz`), and merged into `results/<name/of/sample>/variants/exitrons.vcf.gz`.

### INDEL/SNVs

Two callers feed this path:

- **transIndel** for long indels. The `detect_long_indel_ti_build_DNA` and `detect_long_indel_ti_build_RNA` rules each build a remapped BAM (with redefined CIGAR) before `detect_long_indel_ti_call` extracts the indels; per-group VCFs are augmented with `GRP` / `SRC=long_indel` INFO and merged into `results/<sample>/variants/long.indels.vcf.gz`. Whether the DNA build, RNA build, or both run is set by `indel.mode`.
- **GATK Mutect2** for short indels and SNVs. `detect_short_indels_m2` runs per split BAM (parallel across read-groups), `filter_short_indels_m2` applies Mutect2's own learned filters, and the augment / merge / select rules separate the per-VCF SNV / short-indel streams. Final results land in `results/<sample>/variants/somatic.short.indels.vcf.gz` and `results/<sample>/variants/somatic.snvs.vcf.gz`.

The `indel.type` config key selects which callers run (`short`, `long`, or `all`); `indel.mode` selects the input modality (`DNA`, `RNA`, or `BOTH`) where applicable.

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












