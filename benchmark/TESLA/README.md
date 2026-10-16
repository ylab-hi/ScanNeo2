# TESLA benchmark data

Reference data for benchmarking ScanNeo2's neoepitope recall against experimentally
validated peptides. Kept here rather than under `resources/`, which ScanNeo2 generates
at runtime and is not tracked.

Source: Wells et al., *Key Parameters of Tumor Epitope Immunogenicity Revealed Through a
Consortium Approach Improve Neoantigen Prediction*, Cell 183:818-834.e13 (2020),
<https://doi.org/10.1016/j.cell.2020.09.015>. Open Archive. Both files are the published
supplementary tables, unmodified apart from the filenames.

## The two cohorts are distinct

TESLA reports two separate sets of peptides. They must not be pooled without saying so,
and their column schemas differ.

| | `tesla_main_cohort_tableS4.tsv` | `tesla_validation_cohort_tableS7.tsv` |
|---|---|---|
| supplementary table | S4 | S7 |
| cohort | main TESLA consortium | independent validation |
| subjects | 6 — melanoma 1, 2, 3; NSCLC 10, 12, 16 | 3 melanoma — 4, 8, 9 |
| pMHC tested | 608 | 310 |
| immunogenic | 37 | 4 |
| assay | pMHC multimer, PBMC (melanoma) / TIL (NSCLC) | tetramer, patient-matched PBMC |
| `VALIDATED` encoding | `TRUE` / `FALSE` (column 16) | `1` / `0` (column 13) |
| binding affinity column | `MEASURED_BINDING_AFFINITY` | `PREDICTED_BINDING_AFFINITY` |

S7 also lacks `MHC`, `FRAC_HYDROPHOBIC`, `NUMBER_PREDICTING`, `TCR_FLOW_I` and
`TCR_NANOPARTICLE`; the allele is recoverable from the `PMHC` column prefix.

## Reachable denominator

The benchmark runs over five of the six main-cohort subjects plus all three
validation-cohort subjects.

The subject-to-sample mapping is the `SAMPLES` dict at the top of
`benchmark/tesla_recall.py`, with instructions to edit it. It has to be set by hand: the
sample directory names come from whatever sample sheet configured the run, and nothing
under `results/` records that a given directory is a particular TESLA subject. Delete any
subject whose data you do not have — the denominators are derived from whatever is
listed, so a partial cohort scores correctly against its own totals rather than looking
like a low score. A subject absent from the tables, or one whose
`results/<sample>/prioritization/mhc-I_neoepitopes_all.txt` is missing, is reported and
the run stops before printing any figures.

**Subject 10 is not reachable**: it holds 3 immunogenic pMHC, but its WES/RNA is not
available via Synapse (<https://www.synapse.org/#!Synapse:syn21048999>), so it is excluded,
as it is in published comparisons.

```
main cohort   37 immunogenic - 3 (subject 10) = 34 reachable
validation                                      4 reachable
                                        total  38
```

## Scoring, and why the two classes point in opposite directions

A pMHC counts as found when the mutant peptide **and** the HLA allele both match a row of
`results/<sample>/prioritization/mhc-I_neoepitopes_all.txt`. All 37 + 4 validated peptides
are MHC class I, 8-11 aa, so none falls outside what ScanNeo2 predicts.

- The **immunogenic** pMHC are the positives. Recall over them is the benchmark: higher is
  better.
- The **non-immunogenic** pMHC were tested and found not immunogenic, so they are confirmed
  **true negatives**. Proposing one is a false positive: **lower is better.** Retaining few
  of them means they were already discarded before any prioritisation filter ran. A low
  number here is not a sensitivity gap, and should not be chased.

Reproduce with

```
python benchmark/tesla_recall.py \
    --export benchmark/TESLA/scanneo2_tesla_overlap.tsv
```

Run it from the repository root; the table and `results/` paths are relative to it.

The script prints the denominators it derived from the mapping before scoring anything —
for the eight subjects above, 845 assessed pMHC of which 38 are immunogenic — so a
mangled table or a short mapping is visible in the output rather than silently lowering
the score.

## `scanneo2_tesla_overlap.tsv`

Every ScanNeo2 prediction that overlaps a TESLA-assessed pMHC, committed alongside the
reference tables so the numbers below can be inspected without rerunning the pipeline.
3322 rows: `tesla_subject` and `tesla_validated` prepended to the full neoepitope record,
so each hit carries its gene, source, IC50, VAF, TPM, agretopicity, self-similarity,
NMD and pathogen/proteome-similarity fields.

210 rows are immunogenic and 3112 non-immunogenic; a pMHC can appear on several rows when
the peptide arises from more than one transcript, so de-duplicate on
`(tesla_subject, mt_epitope_seq)` to recover the 36 and 586 counted below.

## Version

Measured with **ScanNeo2 v0.7.1 plus the then-unreleased changes from #236**, which is
`ecfd314` on `main`. Recorded because the score is only meaningful against a known
pipeline version. These numbers should be stable: variant calling, annotation and
prediction are fixed, so only a change that alters which peptides are proposed — a new
neoepitope source, a filter change, a different caller or HLA typer — should move them.
If a routine change does move them, that is worth investigating rather than accepting.

## Results

| | ScanNeo2 | nextNEOpi (Rieder et al. 2022) |
|---|---|---|
| immunogenic recall *(higher better)* | **36/38 (94.7%)** | 36/38 (94.7%) |
| true negatives retained *(lower better)* | **586/807 (72.6%)** | 687/807 (85.1%) |
| true negatives excluded | **221** | 120 |

Identical sensitivity, with 101 more confirmed true negatives discarded before filtering.

Per-subject immunogenic recall is identical to nextNEOpi in every subject: 9/9, 4/4, 12/13,
4/4, 3/4 (main cohort) and 1/1, 1/1, 2/2 (validation cohort). nextNEOpi's figures are its
*unfiltered* candidate set, which is the comparable stage; after its relaxed prioritisation
filter it reports 32/38.

Per-subject true negatives retained: 64/88, 91/104, 74/84, 61/75, 51/107, 64/124, 76/85,
105/140 (subjects 1, 2, 3, 4, 8, 9, 12, 16). Subject 8 is the widest margin — 51 of 107
retained against nextNEOpi's 99 — while still recovering its one immunogenic peptide.

Matching on peptide alone instead of peptide + allele finds 636/845 rather than 622/845, so
14 pMHC match the peptide but not the allele.

The two main-cohort misses are reachability limits rather than pipeline defects:

- **`KIVEMSTSK`** (subject 3) — absent from our output entirely; no supporting variant.
- **`YLNEAVFNFV`** (subject 16) — we call `chr11:100959900 GG>AT`, a phased dinucleotide
  substitution giving `ATG` = Met, so we emit `YLNEAMFNFV`. TESLA's peptide needs only the
  second base changed (`GTG` = Val). A pipeline that calls adjacent substitutions
  independently reproduces TESLA's peptide and scores a match; one that phases them does
  not. nextNEOpi phases too, and misses the same peptide.

Note that nextNEOpi's 36/38 is its *unfiltered* candidate set; after its relaxed
prioritisation filter it reports 32/38. ScanNeo2's `mhc-I_neoepitopes_all.txt` is the
comparable unfiltered set.
