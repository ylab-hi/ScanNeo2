# Config File

In the following the parameters in the `config.yml` are discussed. In principle, the `config.yml` consists of different blocks. ScanNeo2 always utilizes the `config.yml` that is located in `config/config.yml`. However, the option `--configfile` allows specifying custom config files. It is to be noted that merely overwrites the `config/config.yml` and should therefore include all parameters to prevent using settings from multiple config files.

## GENERAL
On the top level, the parameters are applied system-wide when applicable. This includes the reference genome within the `reference` attribute.  The `release` option corresponds to the ENSEMBL release version and `nonchr` indicates to whether (=true) or not (=false) DNA sequence that is not assigned to chromosomes should be included in the analysis. Other options include the number of cores (`threads`), the mapping quality (`mapq`), and the average Phred scores (`basequal`). 

```
threads: 30
mapq: 30  
basequal: 20
```

## SAMPLES

Per-sample inputs live in a separate **sample sheet** (TSV), referenced by `config.yaml`:

```yaml
samples: config/samples.tsv
```

One `snakemake` invocation processes **many samples in parallel** — Snakemake fans out the per-sample DAGs automatically (issue [#93](https://github.com/ylab-hi/ScanNeo2/issues/93)).

The sheet is **long-format**: one row per sequencing input, identified by `(sample, seqtype, group)`. Columns (TAB-separated):

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | sample name; results are written to `results/<sample>/` and used as the `{sample}` wildcard. Repeated across the sample's rows. |
| `seqtype` | seq rows | `dnaseq` or `rnaseq`. Leave empty on a custom-only row. |
| `group` | seq rows | group name, **unique within `(sample, seqtype)`**. Becomes the `{group}` wildcard and the provenance label in the output. |
| `type` | seq rows | `tumor` or `normal`. `normal` marks the matched normal/control (excluded from HLA typing). |
| `reads` | seq rows | one path (single-end) or two space-separated paths (paired-end); `.fq` / `.fastq` / `.bam`. |
| `custom_variants` | no | path to a user-supplied VCF (**sample-level** — see below) |
| `custom_proteins` | no | path to a TSV of `(wildtype, mutant)` protein pairs (sample-level) |
| `custom_hla_I` | no | path to a file listing MHC-I alleles, used when `hlatyping.MHC-I_mode` contains `custom` (sample-level) |
| `custom_hla_II` | no | same for MHC-II, when `hlatyping.MHC-II_mode` contains `custom` (sample-level) |

```tsv
sample	seqtype	group	type	reads	custom_variants	custom_proteins	custom_hla_I	custom_hla_II
P1	dnaseq	tumor	tumor	t_R1.fq.gz t_R2.fq.gz
P1	dnaseq	normal	normal	n_R1.fq.gz n_R2.fq.gz
P1	rnaseq	tumor_rep1	tumor	rep1.bam
P1	rnaseq	tumor_rep2	tumor	rep2.bam
P2					vcfB.vcf.gz			hlaB.tsv
```

`P1` above has a matched normal and **two RNA-seq replicates**; `P2` is a **custom-only** sample (predefined VCF + HLA list, no sequencing).

### Replicates and groups

A **replicate is just another `group`**. Give each replicate a distinct group name within its `(sample, seqtype)` (e.g. `tumor_rep1`, `tumor_rep2`). Every group is aligned and variant-called as an **independent, parallel** branch of the DAG, then all of a sample's candidates are **pooled** at prioritization into one combined list (a `bcftools concat` union — nothing is deduplicated or consensus-filtered). Provenance is preserved: each row of the final neoepitope table carries a `group` column (which replicate/condition it came from) and a `source` column (which caller), so you can pool or split by replicate yourself. On a cluster, the per-group calling distributes across nodes; the final per-sample prioritization is a single (multi-threaded) job.

All groups within one `(sample, seqtype)` must share the same read type (SE/PE) and file type (`.fq`/`.bam`); a mismatch is rejected at load. Put a differing input in its own sample.

### Custom inputs and custom-only samples

The `custom_*` columns are **sample-level**: place them on any of a sample's rows (they must be blank or identical across that sample's rows). A sample with *only* custom inputs (no sequencing) is a single row with `seqtype`/`group`/`type`/`reads` left empty, as `P2` above.

### Migrating from the old `data:` block

The single-sample `data:` block is removed in v0.5.0. To migrate:

1. Create `config/samples.tsv`; use what was `data.name` as the `sample` value.
2. Turn each `data.dnaseq.<group>` / `data.rnaseq.<group>` entry into a row: `seqtype` = `dnaseq`/`rnaseq`, `group` = the old key, `reads` = its path(s), `type` = `normal` for the group named in `data.normal` else `tumor`.
3. Map `data.custom.{variants,proteins}` → `custom_variants` / `custom_proteins` and `data.custom.hlatyping.MHC-{I,II}` → `custom_hla_I` / `custom_hla_II` on any of the sample's rows.
4. Replace the entire `data:` block in `config.yaml` with `samples: config/samples.tsv`.

In `custom_variants`, predefined variants in VCF format can be provided. When available ScanNeo2 utilizes specific INFO keys, which are used in the [results](https://github.com/ylab-hi/ScanNeo2/wiki/Output#prioritization). These include `AO`, `DP`, `AF` which correspond to the observed alleles (supporting reads), the depth of the variant, and the variant allele frequency, respectively.

In `custom_proteins` a TSV of `(wildtype, mutant)` protein pairs can be provided, bypassing variant calling and VEP entirely. This is useful for neoantigen candidates from sources ScanNeo2 does not natively call — other variant callers, RNA editing, proteogenomics, or hand-curated candidates — and for benchmarking with known peptides. The TSV needs a header line. **Required columns**: `id`, `wildtype_protein`, `mutant_protein`. **Optional columns** (each defaults sensibly if absent): `vaf`, `ao`, `dp`, `gene_id`, `gene_name`, `transcript_id`, `chrom`, `group`, `var_type`. Column order is not fixed. Both `wildtype_protein` and `mutant_protein` must be non-empty per row — mutant-only rows are rejected with a clear error because there is no variant region to detect and no wildtype contrast for binding-affinity comparison or self-similarity scoring. Both sequences are truncated at the first `*` or `X` (project convention for a stop codon) before downstream processing. The genomic / transcript / expression output columns (`chrom`, `gene_id`, `TPM`, `NMD`, `PTC_*`, `NMD_escape_rule`) are intentionally left empty for protein input, since they are not recoverable from a raw protein pair. An example TSV ships at [`.tests/integration/data/proteins/proteins.tsv`](https://github.com/ylab-hi/ScanNeo2/blob/main/.tests/integration/data/proteins/proteins.tsv).

In `custom_hla_I` / `custom_hla_II`, user-defined class I and class II alleles can be provided in tab-delimited format. See the [hla section](https://github.com/ylab-hi/ScanNeo2/wiki/Output#hla) in the output wiki page for more information.

These columns are *additive* — ScanNeo2 augments the standard analysis with the user-supplied data unless the corresponding pipeline component (hlatyping, indel, ...) is deactivated.


## PRE-PROCESSING

```
preproc: 
  activate: true  
  minlen: 10
  slidingwindow:
    activate: true
    wsize: 3
```

ScanNeo2 provides an optional pre-processing procedure that is only applied to raw sequencing data. Here, `activate: true` enables the pre-processing, that can be combined with a window trimming from the 3'-end with a defined window size (`wsize`).  Other parameters include the minimum length of the sequencing reads (`minlen`). Note: the globally defined base quality is also applied here.

## ALIGNMENT

```
align:
  chimSegmentMin: 20
  chimScoreMin: 10
  chimJunctionOverhangMin: 10
  chimScoreDropMax: 30
  chimScoreSeparation: 10
```

In principle, the alignment procedure is done differently for DNA- and RNA-seq data. In the case of DNA-seq, the reads are directly aligned using BWA. For RNA-seq data, the sequencing reads are first aligned using STAR followed by realignment with BWA. The reason for that is the variant calling on the transcriptome (e.g., gene fusion, alternative splicing, exitron) which requires splice-aware alignments. Consequently, the parameters in this section control the chimeric alignments and are identical to STAR v2.7.10b. Please refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/STAR_2.7.10b_alpha_230301/doc/STARmanual.pdf) for details. 

| Option | Description |
| ------- | ---------- |
| `chimSegmentMin` | minimum length of chimeric segment length, if ==0, no chimeric output |
| `chimScoreMin` | minimum total (summed) score of the chimeric segments |
| `chimJunctionOverhangMin` | minimum overhang for a chimeric junction |
| `chimScoreDropMax` | max drop (difference) of the chimeric score (the sum of scores of all chimeric segments) from the read length |
| `chimScoreSeparation` | minimum difference (separation) between the best chimeric score and the next one |

## VARIANT CALLING

Each module in the variant calling can be switched on/off using the `activate` (`true` or `false`) property. 

### ALTERNATIVE SPLICING

```
altsplicing:
  activate: true 
  confidence: 3  
  iterations: 5 
  edgelimit: 250  
```

In the detection of alternative splicing events, the parameter `confidence` determines how strongly input alignments are filtered before new nodes and edges are added to the splicing graphs. There are four confidence levels, with confidence increasing from 0 to 3. `iterations` add new intron edges into the splicing a certain number of times. Increasing the value increases the sensitivity, but also the runtime. Shouldn't be set lower than 5. `edgelimit` sets an upper boundary for the maximum number of edges (to reduce its complexity) and limit the runtime.

| Parameter | Value | Description |
| --------- | ----- | ----------- |
| `confidence` | 0-3 | Confidence Interval for the SplAdder with 0 the lowest and 3 the highest confidence |
| `iterations` | 5- | Number of iterations to add new intron edges into the splicing graph |
| `edgelimit` | 250- | Limit the number of edges in the splicing graph |

Please refer to the [SplAdder documentation](https://spladder.readthedocs.io/en/latest/spladder_modes.html) for more details.

### EXITRON SPLICING
```
exitronsplicing:
  activate: true 
  ao: 3  
  pso: 0.05  
  strand: 1 # 0=unstranded, 1=forward, 2=reverse
```

In the exitron splicing, the reported exitrons are controlled by `ao` (allele observed) and `pso` (percent spliced in). The former describes the minimum number of reads that support the exitron, and the latter the minimum cutoff for the exon-exclusion rate. In other words, the ratio of the relative abundance of all isoforms missing a certain exon over the relative abundance of all isoforms of the gene missing the exon.

| Parameter | Value | Description |
| --------- | ----- | ----------- |
| ao | 0- | Minimum number of supporting reads for an exitron event |
| pso | 0.0-1.0 | Minimum cutoff for the percent spliced out index | 


### GENE FUSION
```
genefusion:
  activate: true 
  maxevalue: 0.3
  suppreads: 2  
  maxsuppreads: 1000
  maxidentity: 0.3  
  hpolymerlen: 6  
  readthroughdist: 10000  
  minanchorlen: 20  
  splicedevents: 4  
  maxkmer: 0.6  
  fraglen: 200 
  maxmismatch: 0.01
```
In the detection of gene fusion events, ScanNeo2 utilizes [`Arriba`](https://github.com/suhrig/arriba). Please refer to the tools website for more information. ScanNeo2 allows to specify the following parameters. `maxevalue` is a cutoff for the the number of supporting reads that are expected to have occurred by chance (e-value). Fusion events that exceed this cutoff are discarded. It is to be noted that a high e-value both reports more false positives and increases the runtime dramatically. `suppreads` is a cutoff for the minimal number of supporting reads. This means that fusion events that fall short of this value are discarded. Similarly, `maxsuppreads` defines the maximal number of supporting reads for fusion events which are discarded when exceeding this value. `Arriba` issues a warning when the threshold has been hit (see [logs](https://github.com/ylab-hi/ScanNeo2/wiki/Output#logs)). 


### INDELs/SNVs
```
indel:
  activate: true
  type: short  # long, short, all
  mode: BOTH
  strategy: OPTIMAL_F_SCORE # OPTIMAL_F_SCORE, FALSE_DISCOVERY_RATE, CONSTANT
  fscorebeta: 1.0
  fdr: 0.05
  sliplen: 8
  sliprate: 0.1

```

## HLA GENOTYPING
```
hlatyping:
  class: I # I, II or BOTH
  # the mode (origin) for the genotyping of each class (comma-separated list)
  MHC-I_mode: RNA # DNA, RNA, custom (if empty alleles have to be specified in custom)
  MHC-II_mode: RNA # DNA, RNA, custom (if empty alleles have to be specified in custom)
  # specific path for class II hlatyping (only required when class: II, or BOTH)
  freqdata: ./hlahd_files/freq_data/
  split: ./hlahd_files/HLA_gene.split.txt
  dict: ./hlahd_files/dictionary/

```

`class` selects which MHC classes to type; this also gates the prioritization block below. The two `*_mode` entries take a comma-separated list of read sources used to predict alleles:

| Mode | Behaviour |
|---|---|
| `DNA` | Type from DNA-seq reads (OptiType for class I, HLA-HD for class II) |
| `RNA` | Type from RNA-seq reads (same tools) |
| `custom` | Skip read-based typing; use the file at `data.custom.hlatyping.MHC-{I,II}` |

Combinations are allowed (e.g. `DNA, RNA` runs typing from both sources and merges the alleles). Combining a read-based mode with `custom` (e.g. `DNA, custom`) keeps both sets.

`freqdata`, `split`, and `dict` point at the [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) reference data directories and are only required when class II is being typed via `DNA` or `RNA`. HLA-HD itself must be installed separately and accessible on `PATH` as `hlahd.sh` — the `hlahd.yml` conda env only carries its `bowtie2` dependency. ScanNeo2 checks for all of these at workflow-load time and emits a `[config error]` if anything is missing.

## PRIORITIZATION

```
prioritization:
  class: I # I, II or BOTH
  lengths:
    MHC-I: 8,9,10,11
    MHC-II: 13,14,15
```

`class` selects which MHC classes go through binding-affinity prediction and downstream scoring. It must be a subset of `hlatyping.class` — you can't prioritize a class you haven't typed; ScanNeo2 catches that mismatch at workflow-load time.

`lengths.MHC-I` / `lengths.MHC-II` set the epitope k-mer lengths submitted to the binding predictor (netMHCpan / netMHCIIpan), one prediction job per `(allele, length, wt|mt)` cell. The default ranges (`8,9,10,11` for class I, `13,14,15` for class II) cover the dominant binding-affinity windows for each class. Both discrete lists (`8,9,10,11`) and ranges (`8-11`) are accepted.

## CLUSTER EXECUTION (SLURM)

By default `snakemake --cores all --sdm conda` runs every job locally, in one machine or one interactive allocation. To distribute jobs across a SLURM cluster, ScanNeo2 ships a generic profile at `workflow/profiles/slurm/` (it needs the `snakemake-executor-plugin-slurm`, already in `environment.yml`):

```bash
snakemake --workflow-profile workflow/profiles/slurm --configfile config/config.yaml
```

Snakemake then submits each job with `sbatch`, translating each rule's `threads` into `--cpus-per-task` and the per-rule `runtime` / `mem_mb` from the profile into walltime / memory. With a sample sheet that has multiple samples or replicate groups, the independent per-group branches (alignment, variant calling) fan out across nodes; the per-sample prioritization stays a single job.

The profile is deliberately **cluster-agnostic**: it sets no account and no partition, so jobs land on your cluster's default partition under your default account. Override for your site without editing the file:

```bash
snakemake --workflow-profile workflow/profiles/slurm \
  --default-resources slurm_account=<account> slurm_partition=<partition>
```

or copy the profile and edit `default-resources` / `set-resources`. The memory and runtime tiers (e.g. STAR at 64 GB) are sized for human-scale data as starting points — tune them to your inputs.







