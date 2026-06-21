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

```
samples: config/samples.tsv
```

This pattern lets one `snakemake` invocation process **many samples in parallel** — Snakemake fans out the per-sample DAGs automatically (issue [#93](https://github.com/ylab-hi/ScanNeo2/issues/93)). For a single sample, use a one-row sheet.

The sheet is wide-format, one row per sample. Columns (TAB-separated):

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | sample name; results are written to `results/<sample>/` and used as the `{sample}` wildcard |
| `dnaseq_tumor` | no | DNA-seq tumor reads — one path (single-end) or two space-separated paths (paired-end); `.fq` / `.fastq` / `.bam` |
| `dnaseq_normal` | no | DNA-seq matched normal/control reads — same format |
| `rnaseq` | no | RNA-seq reads — same format |
| `custom_variants` | no | path to a user-supplied VCF |
| `custom_proteins` | no | path to a TSV of `(wildtype, mutant)` protein pairs |
| `custom_hla_I` | no | path to a file listing MHC-I alleles (used when `hlatyping.MHC-I_mode` is `custom`) |
| `custom_hla_II` | no | same, for MHC-II |

Leave a cell empty when not applicable. Empty cells become `None` in the per-sample data dict.

```
sample	dnaseq_tumor	dnaseq_normal	rnaseq	custom_variants	custom_proteins	custom_hla_I	custom_hla_II
sampleA	tumor_R1.fq.gz tumor_R2.fq.gz	normal_R1.fq.gz normal_R2.fq.gz	rna.bam				
sampleB	tumorB.bam			vcfB.vcf.gz			
sampleC				proteins.tsv		hlaC.tsv	
```

The previous single-sample `data:` block is removed in v0.5.0. To migrate an existing config:

1. Create `config/samples.tsv` (or another path) with one row keyed on what was previously `data.name`.
2. Map old keys to columns: `data.dnaseq.<group>` → `dnaseq_tumor` / `dnaseq_normal`; `data.rnaseq.<group>` → `rnaseq`; `data.custom.{variants,proteins}` → `custom_variants` / `custom_proteins`; `data.custom.hlatyping.MHC-{I,II}` → `custom_hla_I` / `custom_hla_II`.
3. Replace the entire `data:` block in `config.yaml` with `samples: config/samples.tsv`.

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







