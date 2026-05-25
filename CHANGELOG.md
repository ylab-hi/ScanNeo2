# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Changed

- **Audit hardcoded `threads:` values; cap and couple to `{threads}`** — 10 rules across `align.smk`, `altsplicing.smk`, `hlatyping_mhcI.smk`, `hlatyping_mhcII.smk`, `quantification.smk`. The most impactful change is **`hlatyping_mhcI_SE/PE` dropping from `threads: 64` to `threads: 1`** — OptiType is single-threaded (ILP solver; the wrapper does not parallelise across cores), so reserving 64 cores per invocation was serializing the downstream workflow for nothing. Other changes: `split_bamfile_RG` capped at `min(10, config["threads"])` (samtools split is I/O-bound; closes #127); `spladder` → `config["threads"]`; `filter_reads_mhcII_SE/PE` → `max(2, config["threads"])` to enforce bowtie2's two-thread minimum; `countfeatures` → `threads: 4`; `dnaseq_postproc` → `threads: 4` (samtools threading plateaus past ~4); and `rnaseq_postproc_fixmate` / `rnaseq_postproc_markdup` / `dnaseq_postproc` have their literal `-@N` shell arguments coupled to `-@ {threads}` so the directive value is what samtools actually receives. ([#127](https://github.com/ylab-hi/ScanNeo2/issues/127), [#138](https://github.com/ylab-hi/ScanNeo2/pull/138))

## [0.3.14] - 2026-05-24

### Refactored

- **Reformat the remaining 13 Snakemake files with `snakefmt 2.0.0`**: mechanical pass (indent 2→4 space, single→double quotes, directive sort, trailing commas, long-line wrap) across `workflow/Snakefile` and every rule file not already cleaned by #124. With this PR the whole 17-file workflow tree is `snakefmt --check` clean. No behavior change; verified end-to-end by `snakefmt --check workflow/`, `snakemake --lint`, and `snakemake --dry-run`. ([#123](https://github.com/ylab-hi/ScanNeo2/issues/123), [#126](https://github.com/ylab-hi/ScanNeo2/pull/126))

### Changed

- **Bump `workflow/envs/samtools.yml` from `samtools=1.14` to `samtools=1.20`** — aligns the env with the version pinned by `v4.0.0` snakemake-wrappers (adopted for bcftools in #129), so samtools usage across the workflow is internally consistent regardless of whether a rule reaches samtools via the env file or a wrapper. Affects 14 rules across `align.smk`, `germline.smk`, `hlatyping_mhcI.smk`, `hlatyping_mhcII.smk`, `hlatyping_prep.smk`, `indel.smk`, and `annotation.smk`. Part of #105 (samtools half of item 3); the deferred wrapper-conversion of the 8 single-command samtools rules is tracked in #133. ([#105](https://github.com/ylab-hi/ScanNeo2/issues/105), [#132](https://github.com/ylab-hi/ScanNeo2/pull/132))
- **Convert 8 single-command bcftools shell rules to `v4.0.0/bio/bcftools/{index,concat}` wrappers** (germline.smk × 5, indel.smk × 3): 6 `bcftools index -t` and 2 `bcftools concat -a` rules drop their `conda: ../envs/bcftools.yml` + `shell:` blocks and gain a `wrapper:` line. The wrapper's env pins **bcftools=1.20** — same as our `bcftools.yml` — so behavior is unchanged. Verified on a full SLURM run: the mutect2 merged VCF (the concat wrapper) is byte-identical to the pre-PR run (728,873 records, body MD5 `98a80f3b861d46adb54cfc741fc75d10`); the merged `.tbi` files cover the same 25 chromosomes and return matching tabix region-query row counts. `bcftools.yml` itself stays for the 15 remaining pipeline rules (`sort | view`, `concat | sort`). Part of #105 (item 3, bcftools subset). ([#105](https://github.com/ylab-hi/ScanNeo2/issues/105), [#129](https://github.com/ylab-hi/ScanNeo2/pull/129))
- **Consolidate conda envs (items 1-2 of #105)**: delete `workflow/envs/picard.yml` — orphan; the one picard call (`rule create_sequence_dictionary` in `ref.smk`) uses the `bio/picard/createsequencedictionary` wrapper which ships its own env — and merge `workflow/envs/realign.yml` into `workflow/envs/basic.yml`. `realign.yml` (channels: `conda-forge`, `bioconda`; deps: `bwa=0.7.17`, `samtools=1.16.1`) was a strict subset of `basic.yml`; the two rules that referenced it (`realign` and `bwa_align_dnaseq` in `align.smk`) only invoke bwa + samtools, so they work identically under `basic.yml`. Items 3 (bcftools/samtools shell → wrapper conversion) and 4 (samtools version sprawl) remain tracked under #105. ([#105](https://github.com/ylab-hi/ScanNeo2/issues/105), [#128](https://github.com/ylab-hi/ScanNeo2/pull/128))
- **Make ScanNeo2 catalog standardized-usage compliant**: `.snakemake-workflow-catalog.yml` was malformed (missing the top-level `usage:` key, deprecated `singularity` naming, `report:` at the wrong nesting level) — almost certainly the reason the Snakemake Workflow Catalog has been flagging us as non-compliant. Rewritten to the current schema with conservative `conda`-only deployment claims, and `desc:`/`flags:` populated. `config/README.md` rewritten from a wiki-pointer into a self-contained parameter reference (per-section tables for every block in `config.yaml`, example `data:` block, invocation snippets); the catalog scrapes this file as the workflow's *Configuration* section. Closes the description / standardized-usage part of #115. ([#115](https://github.com/ylab-hi/ScanNeo2/issues/115), [#125](https://github.com/ylab-hi/ScanNeo2/pull/125))

### Fixed

- **Friendlier config errors at workflow load time**: `handle_seqfiles` and `data_structure` now fail fast with a clear `[config error] ... file(s) not found: <path> (paths are case-sensitive on Linux — check for typos).` message when a configured input path (a `data.dnaseq` / `data.rnaseq` replicate, `data.custom.variants`, or `data.custom.hlatyping.MHC-{I,II}`) doesn't exist on disk. Previously a typo'd or case-mismatched path cascaded into a confusing `MissingInputException` on a downstream rule with the actual root cause buried at the bottom; now snakemake exits before DAG construction with the missing-file message. The exit is immediate (hard fail on the first missing file) — this prevents the mixed-valid-invalid replicate case from silently dropping the broken replicate and proceeding with only the valid one. Same `--lint` guard as PR #122 — linting still completes on the placeholder default config. Closes #130. ([#130](https://github.com/ylab-hi/ScanNeo2/issues/130), [#134](https://github.com/ylab-hi/ScanNeo2/pull/134))
- **Run-summary banner now matches the configured MHC class**: the `HLA typing` and `Prioritization` lines in `print_run_summary` previously printed both `MHC-I` and `MHC-II` columns unconditionally even when only one class was configured (`hlatyping.class: I` / `prioritization.class: I`). Now gated on the active class with the same `("I", "BOTH")` / `("II", "BOTH")` membership pattern already used elsewhere in `common.smk`. Closes #131. ([#131](https://github.com/ylab-hi/ScanNeo2/issues/131), [#134](https://github.com/ylab-hi/ScanNeo2/pull/134))
- **`preproc_paired_end` had a typo in its `params:` block — `dapters=""` instead of `adapters=""`** — so the `adapters` argument was silently dropped before reaching the fastp wrapper. Practical impact zero (the value is `""`, i.e. fastp defaults), but the typo was a real bug. Surfaced during the snakefmt-reformat review of `preproc.smk`. ([#126](https://github.com/ylab-hi/ScanNeo2/pull/126))
- **Make `snakemake --lint` pass on the default config**: the Snakemake Workflow Catalog reported a critical lint error — `IndexError` in `handle_seqfiles` on the default placeholder config — and the workflow could not be statically analysed. Four small changes together let `snakemake --lint --configfile config/config.yaml` exit 0 cleanly: a final-return guard in `handle_seqfiles` against an empty `filetype` list; `data_structure` skips its `sys.exit(1)` when `--lint` is in `sys.argv` (the friendly *"no valid sequence files"* message still prints; real-run UX unchanged); the misleading `<path/to/...>` placeholders in `config/config.yaml` are replaced with commented-out examples; and `wildcard_constraints no=r"\d+"` becomes `r"\d{1,}"` to avoid Snakemake's lint scanner falsely flagging the literal `+` as path composition. Part of #115. ([#115](https://github.com/ylab-hi/ScanNeo2/issues/115), [#122](https://github.com/ylab-hi/ScanNeo2/pull/122))
- **Fix `snakefmt` parse errors and replace super-linter with a pinned `snakefmt` CI step**: three parse-error sites prevented `snakefmt --check` (and the Snakemake Workflow Catalog) from statically analysing the workflow. The multi-line `params: resources={...}` dict in `germline.smk` was extracted to a module-level `VQSR_RESOURCES` constant; the backslash-continued multi-line `params: extra=f"""..."""` in `genefusion.smk` was moved to an `_arriba_extra()` helper in `common.smk` (keeps `genefusion.smk` rule-only so `snakemake --lint` no longer flags *"mixed rules and functions"*); a `((` in `ref.smk` was disambiguated to `( (` and one ```backtick``` command-sub modernised to `$(...)` so shfmt stops parsing the line as arithmetic. The four touched files (`common.smk`, `genefusion.smk`, `germline.smk`, `ref.smk`) were reformatted with `snakefmt 2.0.0`; the remaining 13 are deferred to #123. The CI `formatting` job is now a dedicated step pinning `snakefmt==2.0.0` via `pip` and running `--check` only on Snakemake files changed vs the base branch — useful going forward without forcing a tree-wide reformat now; replaces super-linter, which had been configured to run only snakefmt. Also adds `permissions: contents: read` at workflow level and `persist-credentials: false` on the `actions/checkout` steps as defensive hardening. Part of #115. ([#115](https://github.com/ylab-hi/ScanNeo2/issues/115), [#124](https://github.com/ylab-hi/ScanNeo2/pull/124))
- **`valid_paired_end()` was matching the second filename against `file1` instead of `file2`**, breaking the paired-end pair-validation logic. Pre-existing bug, surfaced during the `snakefmt` review of `common.smk`. The function is currently unused (its sole caller is commented out), so no production impact, but the intent is restored. ([#124](https://github.com/ylab-hi/ScanNeo2/pull/124))

## [0.3.13] - 2026-05-22

### Changed

- **Parallelize binding-affinity prediction in one thread pool**: `collect_binding_affinities` previously ran `len(alleles) × len(epilens)` tasks — each processing its FASTA batches sequentially — and was called once for wt then once for mt, capping concurrency at `min(threads, alleles × epilens)` (e.g. 12 for a 3-allele, 4-length run) regardless of the allocated thread count. It and `calc_binding_affinities` are replaced by a single orchestrator that enumerates every `(group, allele, epitope-length, FASTA-batch)` unit and runs them all in one `ThreadPoolExecutor`, so concurrency scales with `min(threads, total_batches)`. The per-batch computation, global-seqnum keying and cross-allele merge are unchanged — output is identical, only the scheduling differs. ([#113](https://github.com/ylab-hi/ScanNeo2/issues/113), [#118](https://github.com/ylab-hi/ScanNeo2/pull/118))

### Fixed

- **Re-enable the variant-overlap filter in prioritization output**: `prediction.py`'s output-assembly loop emitted every epitope of a variant subsequence — including ones lying entirely up-/downstream of the mutation, which are pure wildtype (`mt_epitope_seq == wt_epitope_seq`) and not neoantigens. The filter meant to drop them was commented out due to a coordinate-frame mismatch. It is re-enabled in the `mt_subseq` frame: every occurrence of the epitope in `mt_subseq` is enumerated, and the epitope is kept only if some occurrence spans `[aa_var_start, aa_var_end]` — that occurrence also feeds the `wt_epitope_seq` slice. Enumerating (rather than first-occurrence `find()`) correctly handles a k-mer that repeats within `mt_subseq`. ([#108](https://github.com/ylab-hi/ScanNeo2/issues/108), [#114](https://github.com/ylab-hi/ScanNeo2/pull/114))

## [0.3.12] - 2026-05-22

### Fixed

- **Fix epitope-to-variant mis-association from batched prediction FASTA**: `predict_binding.py` numbers the sequences of each FASTA file it receives starting from 1, ignoring the `>N` header. `split_fasta_into_batches` split the prediction FASTA into 500-sequence batches, so every batch was numbered independently and `calc_binding_affinities` collapsed every global record `K, K+500, …` onto key `K`. The output-assembly loop looks variants up by their global sequence number, so variants past the first 500 had their binding affinities — and epitope sequences — attached to the wrong variant. `split_fasta_into_batches` now returns each batch's global offset and `calc_binding_affinities` translates the per-batch sequence number to a global one (`offset + per-batch number`) before keying. ([#109](https://github.com/ylab-hi/ScanNeo2/issues/109), [#111](https://github.com/ylab-hi/ScanNeo2/pull/111))

## [0.3.11] - 2026-05-21

### Changed

- **Pin loose conda dependencies**: Pinned the remaining unpinned packages across 6 `workflow/envs/*.yml` files (`curl`, `bowtie2`, `pip`, `python`, `perl`, `perl-env`, `bwa`, `samtools`, `pysam`) to versions verified by dry-run solving each env. `realign.yml` also gained the `conda-forge` channel — without it the solver starved modern samtools of its dependencies and fell back to samtools 1.3.1 — and `manipulate_vcf.yml` had its channels reordered to `conda-forge, bioconda` (dropping `defaults`, which had resolved `pysam` to 0.9.1 and conflicted on `libdeflate`). ([#64](https://github.com/ylab-hi/ScanNeo2/issues/64), [#101](https://github.com/ylab-hi/ScanNeo2/pull/101))
- **Improve documentation**: Fixed the self-contradictory Snakemake-version note and a `config.yml`→`config.yaml` typo in the README, and clarified that HLA-HD's `hlahd.sh` must be on `PATH`. Added inline comments to the previously-undocumented `config/config.yaml` options (`reference`, `data`, `preproc`, the STAR chimeric-alignment parameters, `quantification`, the HLA-HD reference-data paths, `prioritization` lengths). Added module-level docstrings to all 8 `workflow/scripts/prioritization/` modules. ([#71](https://github.com/ylab-hi/ScanNeo2/issues/71), [#102](https://github.com/ylab-hi/ScanNeo2/pull/102))
- **Standardize STAR, samtools, and fastqc wrapper versions**: Converged the Snakemake `wrapper:` versions for the inconsistent families — all three STAR wrappers to `v2.2.1` (`star/align` had been pinned at two different versions) and all samtools wrappers to `v2.3.0` — and bumped the 5 fastqc wrappers to `v9.8.0`. Verified tool-version-neutral by comparing each wrapper's `environment.yaml`: `star` stays 2.7.10b, `samtools` 1.17, `fastqc` 0.12.1, so only the wrapper glue script changes. ([#69](https://github.com/ylab-hi/ScanNeo2/issues/69), [#103](https://github.com/ylab-hi/ScanNeo2/pull/103))
- **Structured startup run summary and clearer config errors**: Replaced the raw `print(config['data'])` dict dump with `print_run_summary()` — a boxed, aligned block printed once at the start of a run showing the run name, reference release, threads, output directory, per-mode input data (sample IDs, file type, read type, resolved paths), variant-calling modes, HLA typing class/modes, and prioritization class/epitope lengths. `handle_seqfiles()` now reports an invalid input as `data.<mode>.<replicate>: ...` and notes when a path is likely an unfilled config placeholder; the no-data abort message is likewise made specific and all config errors go to stderr. ([#77](https://github.com/ylab-hi/ScanNeo2/issues/77), [#106](https://github.com/ylab-hi/ScanNeo2/pull/106))
- **Add JSON Schema validation for the config file**: New `workflow/schemas/config.schema.yaml` (draft-07) validates all 14 config sections — field types, required keys, `enum`s for the fixed-value fields (`indel.type`/`mode`/`strategy`, `quantification.mode`, the `class` fields), nullable types for the optional `data.*` paths, and `additionalProperties: false` to catch mistyped keys. The Snakefile calls `validate(config, ...)` immediately after loading the config, so a missing/mistyped key or wrong value type fails fast with a clear message instead of a cryptic downstream `KeyError`. Completes the schema half of #65 (the wildcard-constraints half landed in #81). ([#65](https://github.com/ylab-hi/ScanNeo2/issues/65), [#107](https://github.com/ylab-hi/ScanNeo2/pull/107))

### Fixed

- **Strip the wildtype-padding sentinel before it reaches epitope output**: `prioritization/effects.py` pads the wildtype protein with `$` so it length-matches the mutant for position-wise variant-boundary detection. `prediction.py`'s FASTA-writing loop stripped `$`, but its output-assembly loop sliced the raw padded `wt_subseq` straight into the `wt_epitope_seq` column — leaking `$` into the neoepitopes table. IEDB's `predict_immunogenicity.py` rejects any non-standard residue with `sys.exit(1)` (reported on stdout), so this surfaced as a `prioritization` rule failure once #92 added `check=True`. `$` is now stripped once when the variant-effects tsv line is parsed, so both loops see a clean `wt_subseq`; `SequenceSimilarity.self_similarity` skips on a wt/mt length mismatch instead of detecting `$`; and `calc_immunogenicity_mhcI` filters its input to non-empty standard-amino-acid peptides — a fully novel epitope now yields an empty wt sequence, which would still trip the IEDB tool's empty-list `UnboundLocalError`. ([#107](https://github.com/ylab-hi/ScanNeo2/pull/107))

## [0.3.10] - 2026-05-20

### Fixed

- **Drop `tmp/` as rule input to stop spurious reruns**: Three rules (`rnaseq_postproc_markdup`, `dnaseq_postproc`, `get_reads_hlatyping_BAM`) declared `tmp="tmp/"` as input, but `samtools sort -T tmp/` and intermediate BAMs in those same rules wrote into the directory — bumping its mtime so Snakemake's mtime trigger re-fired the rules on every invocation. Each rule now `mkdir -p tmp/` itself, removing the input dependency entirely; the now-orphan `create_tmp_folder` rule and `workflow/rules/prelim.smk` are deleted, and the `include:` line in `workflow/Snakefile` is removed. Each `samtools sort -T` now uses a wildcard-tagged prefix (`tmp/sort_{sample}_{group}_`, `tmp/sort_{sample}_{group}_{nartype}_`) so concurrent jobs' shards are visibly attributable. ([#31](https://github.com/ylab-hi/ScanNeo2/issues/31), [#91](https://github.com/ylab-hi/ScanNeo2/pull/91))
- **Wrap shell pipelines so stderr from upstream commands lands in the rule log**: Bash binds `> {log} 2>&1` only to the last command in a pipeline, so upstream stderr (`bwa mem`, `samtools collate`/`fastq`, `bcftools sort`/`concat`, `yara_mapper`) was leaking to Snakemake's main output and obscuring its own progress messages. Wrapped 21 multi-command pipelines in `( ... ) > {log} 2>&1` across 8 rule files: 4 in `align.smk` (`rnaseq_postproc_fixmate`, `realign`, `bwa_align_dnaseq`, `dnaseq_postproc`), 7 in `indel.smk`, 3 in `germline.smk`, 2 each in `altsplicing.smk`, `exitron.smk`, `hlatyping_mhcI.smk`, and 1 in `custom.smk`. The two `(curl ... | tar xz ...)` pipelines in `prioritization.smk` were already correctly wrapped. Verified by re-running `bwa_align_dnaseq` + `dnaseq_postproc` on the integration test: main output clean, every command's stderr now in `{log}`. ([#88](https://github.com/ylab-hi/ScanNeo2/issues/88), [#90](https://github.com/ylab-hi/ScanNeo2/pull/90))
- **Fix two NaN crashes in `prioritization/filtering.py`**: (1) Removed a broken-and-unused `$`-filter in `Immunogenicity.__init__` — `wt_epitope_seq.str.contains("\$")` returns `NaN` for rows where the wildtype peptide didn't meet length requirements, and `~NaN` raises `TypeError: bad operand type for unary ~: 'float'`. The filter's result was also never used (next line called `calc_immunogenicity_mhcI(wt_epitope_seq)` with the unfiltered original) and the `$` marker is already stripped upstream in `prediction.py`. (2) Added a `pd.isna(wt_seq) or pd.isna(mt_seq)` guard to `SequenceSimilarity.self_similarity` — `'$' in NaN` raised `TypeError: argument of type 'float' is not iterable`. NaN entries now get -1 (matching the existing skip value for `$`-found entries). Both surfaced by the integration QC test for #87. ([#87](https://github.com/ylab-hi/ScanNeo2/pull/87))

### Changed

- **Mark large intermediates as `temp()`**: Annotation-only sweep across rule files (no logic changes). `temp()` on STAR aligned/fixmate BAMs, transindel build BAMs, BWA per-chr split directory + indexes, and per-chromosome htcaller/mutect2 VCFs and their bgzipped/indexed siblings — Snakemake auto-deletes once no downstream rule needs them. ([#67](https://github.com/ylab-hi/ScanNeo2/issues/67), [#87](https://github.com/ylab-hi/ScanNeo2/pull/87)) (`protected()` on the final outputs was also added here but reverted before release — see [#92](https://github.com/ylab-hi/ScanNeo2/pull/92).)
- **Harden Python file handling and subprocess calls**: Converted bare `open()` to `with` / `contextlib.ExitStack` across the prioritization, genotyping, and quantification scripts, and made `prioritization/effects.py:VariantEffects` a context manager so file descriptors are released even on a mid-script exception. Added `check=True` + `stderr=subprocess.PIPE` + error logging to the `filtering.py` subprocess calls (blastp, makeblastdb, IEDB immunogenicity) that previously swallowed non-zero exits, with an empty-input guard in `calc_immunogenicity_mhcI` so a variant type with zero neoepitopes no longer crashes the IEDB immunogenicity tool (which raises `UnboundLocalError` on an empty peptide list). Deleted the dead top-level `predict_affinities.py` / `predict_immunogenicity.py`, superseded by the `prioritization/` package. Also fixed a missing-newline bug in `finalize_mhcII_input.py` that produced malformed reverse-read FASTQ (the sequence line concatenated with the following `+` separator). ([#68](https://github.com/ylab-hi/ScanNeo2/issues/68), [#92](https://github.com/ylab-hi/ScanNeo2/pull/92))

## [0.3.9] - 2026-05-17

### Added

- **Pytest scaffolding and unit tests for `combine_all_alleles.py`**: First testing-infrastructure PR. Adds `pytest.ini`, `.github/workflows/test.yml` (runs `pytest .tests/unit/` on Python 3.12), `.gitignore`, and 6 unit tests covering the script's happy path, allele-field truncation, off-refset rejection, empty-input diagnosis, malformed-line rejection, and mixed-source breakdown. Locks in the diagnostic behavior added in PR #83. Hand-rolled rather than auto-generated because naive `snakemake --generate-unit-tests` produces 6+ GB fixtures (vendored IEDB tools declared as rule inputs). ([#85](https://github.com/ylab-hi/ScanNeo2/pull/85))

### Changed

- **Harden `_run_prediction` subprocess call**: Added `timeout=3600s` (per-batch wall-clock cap), `stderr=subprocess.PIPE`, and `check=True` to the netMHCpan/netMHCIIpan invocation in `prediction.py`. Wraps the call in `try/except subprocess.TimeoutExpired / subprocess.CalledProcessError`, logs a warning, and returns an empty dict so a single hung or failing batch doesn't stall the whole prioritization run. ([#59](https://github.com/ylab-hi/ScanNeo2/issues/59), [#75](https://github.com/ylab-hi/ScanNeo2/pull/75))
- **Switch `optitype_wrapper.py` to multiple positional BAM args**: Replaced the single space-joined `sys.argv[1]` BAM string with variadic trailing positional args (`<nartype> <prefix> <outpath> <bam1> [bam2] [...]`) and updated both `hlatyping_mhcI_SE` / `hlatyping_mhcI_PE` rule call sites in `hlatyping_mhcI.smk` accordingly. Added a usage-string guard for misinvocations. Wrapper-side half of issue #76; the shell-side `:q` quoting needed to complete the spaces-in-paths goal is tracked in [#80](https://github.com/ylab-hi/ScanNeo2/issues/80). ([#76](https://github.com/ylab-hi/ScanNeo2/issues/76), [#78](https://github.com/ylab-hi/ScanNeo2/pull/78))
- **Add global `wildcard_constraints` to Snakefile**: Constrained the 11 wildcards used across the workflow — binary enums for `seqtype`/`nartype`/`readtype`/`readpair`, `\d+` for the GATK split-index `no`, and `[^/]+` (path-segment-only) for `sample`/`group`/`vartype`/`chr`/`rg`/`replicate`. Prevents ambiguous rule resolution (e.g., greedy `{group}_{nartype}` matches) and path-traversal-shaped wildcard values reaching outputs. Wildcards half of [#65](https://github.com/ylab-hi/ScanNeo2/issues/65); schema validation still open. ([#65](https://github.com/ylab-hi/ScanNeo2/issues/65), [#81](https://github.com/ylab-hi/ScanNeo2/pull/81))
- **Clean up dead code, unused imports, and a temp-file leak**: Removed commented-out blocks from `add_infos_to_vcf.py`, `predict_affinities.py`, and `prioritization/prediction.py`; dropped unused `import pdb` from `prioritization/variants.py`. Replaced `NamedTemporaryFile(delete=False) + samtools view -o` in `featurecounts_wrapper.py` with `samtools view -c` — removes the temp-file leak source entirely and avoids writing a potentially large filtered BAM just to check whether it's empty. Replaced `> {log} 2>&1 | exit 0` with `> {log} 2>&1 || true` in `hlatyping_mhcII.smk:95` — the original form doesn't actually swallow exit codes under Snakemake's `set -euo pipefail` strict-mode shell. ([#70](https://github.com/ylab-hi/ScanNeo2/issues/70), [#84](https://github.com/ylab-hi/ScanNeo2/pull/84))

### Fixed

- **Improve `combine_all_alleles.py` diagnostics on zero-alleles failure**: When the script finds zero valid alleles across all input sources, the error output now includes a per-source breakdown (file path, byte size, line count, refset match count) and concrete upstream-path pointers (HLA-filtered BAMs, per-batch OptiType results, user-allele config). Distinguishes "0-byte input" (upstream produced nothing) from "rows present but none matched the refset". Hard-fail behavior preserved — running prioritization with zero alleles produces no meaningful output, so the workflow stop is correct; the fix is purely about making the rule log self-explanatory. ([#79](https://github.com/ylab-hi/ScanNeo2/issues/79), [#83](https://github.com/ylab-hi/ScanNeo2/pull/83))
- **Apply `:q` to all BAM placeholders in `hlatyping_mhcI.smk`**: Added Snakemake's shell-quoting formatter to every `{input.*}` / `{output.*}` BAM placeholder in shell blocks across the file (12 placeholder edits in 8 sites across 6 rules: `filter_reads_mhcI_{SE,PE}`, `sort_reads_mhcI_SE` / `sort_and_index_reads_mhcI_PE`, `split_reads_mhcI_{SE,PE}`, `hlatyping_mhcI_{SE,PE}`). Closes the second half of the spaces-in-paths fix that #78 started — without `:q`, Snakemake substitutes paths into the shell string literally and bash tokenizes on whitespace, so a BAM path with a space would still break even after #78's wrapper-side argv change. ([#80](https://github.com/ylab-hi/ScanNeo2/issues/80), [#82](https://github.com/ylab-hi/ScanNeo2/pull/82))
- **Remove spurious `.bai` input from `split_reads_mhcI_SE`**: The checkpoint declared a `.bai` input for `*_flt_SE_sorted.bam`, but no rule in the workflow produces one — `MissingInputException` blocked DAG resolution. The PE counterpart (`split_reads_mhcI_PE`) declares no `.bai` and runs the same `gatk SplitSamByNumberOfReads` call fine; the upstream sort uses `samtools sort -n` (queryname), so a coordinate index cannot exist for this BAM anyway, and `SplitSamByNumberOfReads` reads sequentially without needing one. Surfaced by the integration QC test for #78. ([#78](https://github.com/ylab-hi/ScanNeo2/pull/78))
- **Replace shell=True subprocess calls with list-based arguments**: Refactored 5 wrapper scripts (`optitype_wrapper`, `finalize_mhcII_input`, `get_readgroups`, `compile`, `featurecounts_wrapper`) to use list-based `subprocess.run(..., check=True)` instead of shell=True string commands. Eliminates path-with-spaces breakage, surfaces tool stderr, and lets failures fail loudly. `compile.combine_neoepitopes` rewritten to concatenate files in Python; `get_readgroups.scan_bamfile` filters `@RG` headers in Python instead of piping samtools to grep. ([#62](https://github.com/ylab-hi/ScanNeo2/issues/62), [#74](https://github.com/ylab-hi/ScanNeo2/pull/74))
- **Fix dead VQSR download URLs and add --fail to all curl rules**: The Broad's hg38 GATK resources moved from `genomics-public-data/resources/broad/hg38/v0/` (now access-denied) to `gcp-public-data--broad-references/hg38/v0/`; the old URLs caused `curl` to silently save a 298-byte XML 403 page as the VCF, and GATK VariantRecalibrator then failed with "no suitable codecs found". Updated the 4 affected URLs and added `--fail` to every curl call across 5 rule files (23 invocations) so download failures surface at the download step. Contributes to [#63](https://github.com/ylab-hi/ScanNeo2/issues/63) (`--fail` half complete; checksum verification still open). ([#74](https://github.com/ylab-hi/ScanNeo2/pull/74))

## [0.3.8] - 2026-03-26

### Fixed

- **Batch binding affinity predictions**: Split large FASTA files into batches of 500 sequences before sending to NetMHCpan/NetMHCIIpan, preventing indefinite hangs on large inputs (e.g., 3,700+ alternative splicing events). Also fixed `as_completed` loop running outside the executor context and added detailed progress logging. ([#58](https://github.com/ylab-hi/ScanNeo2/pull/58))
- **Fix malformed params and string comparisons**: Fixed `align.smk` params that passed a literal string instead of the config MAPQ value, and replaced `is not ""` with `!=` across `compile.py` and `merge_counttables.py` (future Python 3.12 SyntaxError). ([#61](https://github.com/ylab-hi/ScanNeo2/issues/61), [#72](https://github.com/ylab-hi/ScanNeo2/pull/72))

### Refactored

- **Standardize log paths and fix missing log redirects**: Unified log directory structure across all 138 rules (`logs/ref/`, `logs/download/`, `logs/{sample}/{stage}/`), added missing `> {log} 2>&1` redirects to ~40 shell rules, and fixed several log path bugs (typos, duplicates, SE/PE mismatches). ([#60](https://github.com/ylab-hi/ScanNeo2/pull/60))

## [0.3.7] - 2025-11-01

### Fix
- Improved error logging and output capture for sample processing steps.
- Optimized the workflow by sorting input files by QNAME before splitting 
in mhc-I genotyping, and removing the indexing step

## [0.3.6] - 2025-10-30

### Fixes
- Fixed bug in prioritization of MHC-II variants. Added wrong perl path to IEDB tools
- Enhanced validation to skip malformed amino acid change entries

### Chores
- Updated workflow environment configuration with additional dependencies


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
