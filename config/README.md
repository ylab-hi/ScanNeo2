# Configuring ScanNeo2

ScanNeo2 is configured through a single YAML file, `config/config.yaml`. The
default file ships with sensible values for every section; you typically only
need to edit the `data:` section to point at your inputs.

Run the workflow with the default configuration:

```bash
snakemake --cores all --software-deployment-method conda
```

…or with a custom configuration file at any path:

```bash
snakemake --cores all --software-deployment-method conda \
          --configfile /path/to/my-config.yaml
```

Paths in the config file are resolved **relative to the directory from which
you invoke `snakemake`**.

The sections below mirror the structure of `config/config.yaml`. The schema is
defined in [`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml);
its parameter table is rendered automatically by the Snakemake Workflow Catalog.

---

## `reference` — reference genome & annotation

| Key | Type | Default | Description |
|---|---|---|---|
| `reference.release` | int | `111` | Ensembl release used to download the reference genome and annotation. |
| `reference.nonchr` | bool | `false` | Include non-chromosomal / scaffold contigs in the reference. |

## `threads`, `mapq`, `basequal` — global limits

| Key | Type | Default | Description |
|---|---|---|---|
| `threads` | int | `30` | Upper bound on threads any single rule may use. Effective threads per rule are `min(rule.threads, --cores)`. |
| `mapq` | int | `30` | Minimum read mapping quality (MAPQ, Phred-scaled) used during filtering. |
| `basequal` | int | `20` | Minimum base-call quality (Phred-scaled) used during filtering. |

## `data` — input samples

The `data` section is the only block that **must** be edited.

| Key | Type | Description |
|---|---|---|
| `data.name` | str | Run name. Results are written to `results/<name>/`. |
| `data.dnaseq` | mapping | DNA-seq inputs, one entry per group. Key = group name; value = one path (single-end) or two space-separated paths (paired-end). Accepted extensions: `.fastq` / `.fq` / `.bam`. |
| `data.rnaseq` | mapping | RNA-seq inputs, same format as `data.dnaseq`. |
| `data.normal` | str | Group name of the matched normal/control sample (used for somatic calling). Leave empty if no matched normal exists. |
| `data.custom.variants` | path | Optional path to a user-supplied VCF to prioritize directly (bypassing variant calling). |
| `data.custom.hlatyping.MHC-I` | path | Optional path to a file listing MHC class I alleles (used when `hlatyping.MHC-I_mode: custom`). |
| `data.custom.hlatyping.MHC-II` | path | Optional path to a file listing MHC class II alleles (used when `hlatyping.MHC-II_mode: custom`). |

Example `data:` block:

```yaml
data:
  name: my_run
  dnaseq:
    tumor:  /path/to/dna_tumor_R1.fq.gz /path/to/dna_tumor_R2.fq.gz
    normal: /path/to/dna_normal.bam
  rnaseq:
    tumor:  /path/to/rna_tumor_R1.fq.gz /path/to/rna_tumor_R2.fq.gz
  normal: normal
```

## `preproc` — fastq pre-processing

Applied only to FASTQ inputs (BAM inputs are not re-trimmed).

| Key | Type | Default | Description |
|---|---|---|---|
| `preproc.activate` | bool | `true` | Whether to run pre-processing. |
| `preproc.minlen` | int | `10` | Discard reads shorter than this (bp) after trimming. |
| `preproc.slidingwindow.activate` | bool | `true` | Enable sliding-window quality trimming. |
| `preproc.slidingwindow.wsize` | int | `3` | Sliding-window size (bp). |
| `preproc.slidingwindow.wqual` | int | `20` | Mean base quality required within the window (Phred-scaled). |

## `align` — STAR chimeric alignment

Parameters passed to STAR for RNA-seq gene-fusion detection.

| Key | Type | Default | Description |
|---|---|---|---|
| `align.chimSegmentMin` | int | `20` | Minimum length of each chimeric segment (0 disables chimeric detection). |
| `align.chimScoreMin` | int | `10` | Minimum total chimeric-alignment score. |
| `align.chimJunctionOverhangMin` | int | `10` | Minimum overhang length for a chimeric junction. |
| `align.chimScoreDropMax` | int | `30` | Maximum allowed drop of chimeric score below read length. |
| `align.chimScoreSeparation` | int | `10` | Minimum score separation between best and next-best chimeric alignment. |

## `altsplicing` — alternative splicing events

| Key | Type | Default | Description |
|---|---|---|---|
| `altsplicing.activate` | bool | `true` | Include alternative-splicing events. |
| `altsplicing.confidence` | int | `3` | Confidence level (1–3) for filtering input alignments. |
| `altsplicing.iterations` | int | `5` | Number of intron-edge addition iterations (sensitivity vs runtime). |
| `altsplicing.edgelimit` | int | `250` | Maximum edges in the splice graph. |

## `exitronsplicing` — exitron splicing events

| Key | Type | Default | Description |
|---|---|---|---|
| `exitronsplicing.activate` | bool | `true` | Include exitron-splicing events. |
| `exitronsplicing.ao` | int | `3` | Allele observation count. |
| `exitronsplicing.pso` | float | `0.05` | Percent spliced-out. |
| `exitronsplicing.strand` | int | `1` | Library strand specificity (0=unstranded, 1=forward, 2=reverse). |

## `genefusion` — Arriba gene-fusion calling

| Key | Type | Default | Description |
|---|---|---|---|
| `genefusion.activate` | bool | `true` | Include gene-fusion events. |
| `genefusion.maxevalue` | float | `0.3` | Maximum E-value. |
| `genefusion.suppreads` | int | `2` | Minimum supporting reads (fusions below this are discarded). |
| `genefusion.maxsuppreads` | int | `1000` | Maximum supporting reads. |
| `genefusion.maxidentity` | float | `0.3` | Genes with identity above this fraction are treated as homologs and discarded. |
| `genefusion.hpolymerlen` | int | `6` | Remove breakpoints adjacent to homopolymers of this length. |
| `genefusion.readthroughdist` | int | `10000` | Distance (bp) below which adjacent breakpoints are considered read-through. |
| `genefusion.minanchorlen` | int | `20` | Discard fusions whose segments are shorter than this. |
| `genefusion.splicedevents` | int | `4` | Fusions between genes require at least this many spliced breakpoints. |
| `genefusion.maxkmer` | float | `0.6` | Remove reads where a repetitive 3-mer makes up more than this fraction. |
| `genefusion.fraglen` | int | `200` | Mean fragment length. |
| `genefusion.maxmismatch` | float | `0.01` | Maximum mismatch fraction. |

## `indel` — small variant calling

| Key | Type | Default | Description |
|---|---|---|---|
| `indel.activate` | bool | `true` | Include indels & SNVs. |
| `indel.type` | str | `all` | `long`, `short`, or `all`. |
| `indel.mode` | str | `BOTH` | Call from `DNA`, `RNA`, or `BOTH`. |
| `indel.strategy` | str | `OPTIMAL_F_SCORE` | Posterior-probability threshold strategy: `OPTIMAL_F_SCORE`, `FALSE_DISCOVERY_RATE`, or `CONSTANT`. |
| `indel.fscorebeta` | float | `1.0` | Relative weight of recall to precision (used with `OPTIMAL_F_SCORE`). |
| `indel.fdr` | float | `0.05` | False-discovery rate target (used with `FALSE_DISCOVERY_RATE`). |
| `indel.sliplen` | int | `8` | Minimum reference bases to suspect a slippage event. |
| `indel.sliprate` | float | `0.1` | Suspected slippage frequency. |

## `quantification` — expression quantification

| Key | Type | Default | Description |
|---|---|---|---|
| `quantification.mode` | str | `BOTH` | Quantify from `DNA`, `RNA`, or `BOTH`. |

## `hlatyping` — HLA typing

| Key | Type | Default | Description |
|---|---|---|---|
| `hlatyping.class` | str | `BOTH` | HLA class to type: `I`, `II`, or `BOTH`. |
| `hlatyping.MHC-I_mode` | str | `DNA, RNA` | Source(s) for MHC-I typing: `DNA`, `RNA`, or `custom` (use `data.custom.hlatyping.MHC-I`). |
| `hlatyping.MHC-II_mode` | str | `DNA, RNA` | Source(s) for MHC-II typing. |
| `hlatyping.freqdata` | path | `./hlahd_files/freq_data/` | HLA-HD frequency-data directory (only required when class II is enabled). |
| `hlatyping.split` | path | `./hlahd_files/HLA_gene.split.txt` | HLA-HD gene-split file. |
| `hlatyping.dict` | path | `./hlahd_files/dictionary/` | HLA-HD dictionary directory. |

HLA-HD must be installed and on `PATH` if class II is enabled — see the
top-level [README](../README.md) for installation notes.

## `prioritization` — epitope prediction

| Key | Type | Default | Description |
|---|---|---|---|
| `prioritization.class` | str | `I` | Epitope MHC class to predict: `I`, `II`, or `BOTH`. |
| `prioritization.lengths.MHC-I` | str | `8,9,10,11` | Comma-separated MHC-I epitope lengths (aa). |
| `prioritization.lengths.MHC-II` | str | `13,14,15` | Comma-separated MHC-II epitope lengths (aa). |

---

For tutorials and worked examples, see the [ScanNeo2 wiki](https://github.com/ylab-hi/ScanNeo2/wiki).
