import os
import shutil
import sys
import glob
from pathlib import Path


def _arriba_extra():
    """Build the `extra=` argument for the arriba wrapper.

    Held here (rather than in genefusion.smk) so genefusion.smk stays
    rule-only; mixing functions and rules in the same .smk trips
    `snakemake --lint`.
    """
    return " ".join(
        [
            "-G 'gene_name=gene_name|gene_id gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS'",
            f"-E {config['genefusion']['maxevalue']}",
            f"-S {config['genefusion']['suppreads']}",
            f"-L {config['genefusion']['maxidentity']}",
            f"-H {config['genefusion']['hpolymerlen']}",
            f"-R {config['genefusion']['readthroughdist']}",
            f"-A {config['genefusion']['minanchorlen']}",
            f"-M {config['genefusion']['splicedevents']}",
            f"-K {config['genefusion']['maxkmer']}",
            f"-F {config['genefusion']['fraglen']}",
            f"-V {config['genefusion']['maxmismatch']}",
        ]
    )


########### CONFIG ##########
def per_sample_data(row):
    """Build the per-sample data dict for one row of the sample sheet.

    Returns a dict matching the legacy `config['data']` shape (dnaseq,
    dnaseq_filetype, dnaseq_readtype, rnaseq, rnaseq_filetype, rnaseq_readtype,
    normal, custom.{variants,proteins,hlatyping.{MHC-I,MHC-II}}) so every
    existing per-rule consumer can swap `config['data'][X]` for
    `SAMPLES[wildcards.sample][X]` without further changes.

    Wide-sheet group convention: the `dnaseq_tumor` cell becomes a `dna_tumor`
    entry in the dnaseq dict; `dnaseq_normal` becomes `dna_normal`; `rnaseq`
    becomes `rna_tumor`. Cells may carry one path (single-end) or two
    space-separated paths (paired-end), matching `handle_seqfiles` parsing.
    """
    sample = row["sample"]

    dnaseq_raw = {}
    if row.get("dnaseq_tumor", ""):
        dnaseq_raw["dna_tumor"] = row["dnaseq_tumor"]
    if row.get("dnaseq_normal", ""):
        dnaseq_raw["dna_normal"] = row["dnaseq_normal"]

    rnaseq_raw = {}
    if row.get("rnaseq", ""):
        rnaseq_raw["rna_tumor"] = row["rnaseq"]

    # handle_seqfiles is pure; it validates + normalises one seqdict at a time
    dnaseq, dnaseq_filetype, dnaseq_readtype = handle_seqfiles(
        dnaseq_raw or None, f"sample {sample!r} dnaseq"
    )
    rnaseq, rnaseq_filetype, rnaseq_readtype = handle_seqfiles(
        rnaseq_raw or None, f"sample {sample!r} rnaseq"
    )

    normal = "dna_normal" if "dna_normal" in dnaseq else None

    def cell(name):
        v = row.get(name, "")
        return v if v else None

    custom_variants = cell("custom_variants")
    custom_proteins = cell("custom_proteins")
    custom_hla_I = cell("custom_hla_I")
    custom_hla_II = cell("custom_hla_II")

    # check that any user-supplied non-fastq paths actually exist on disk
    custom_paths = [
        (f"sample {sample!r} custom_variants", custom_variants),
        (f"sample {sample!r} custom_proteins", custom_proteins),
        (f"sample {sample!r} custom_hla_I", custom_hla_I),
        (f"sample {sample!r} custom_hla_II", custom_hla_II),
    ]
    for label, path in custom_paths:
        if path is not None and not Path(path).is_file():
            print(
                f"[config error] {label}: file not found: '{path}' "
                f"(paths are case-sensitive on Linux — check for typos).",
                file=sys.stderr,
            )
            if "--lint" not in sys.argv:
                sys.exit(1)

    # abort if no data could be found for this sample
    if len(dnaseq) == 0 and len(rnaseq) == 0:
        if custom_variants is None and custom_proteins is None:
            print(
                f"[config error] sample {sample!r}: no valid sequence files found and "
                "no custom variants or proteins provided -- nothing to run for this "
                "sample. Check the corresponding row in the sample sheet.",
                file=sys.stderr,
            )
            # skip the abort under `snakemake --lint` so static rule analysis can
            # still complete on a placeholder/empty default config (e.g. for the
            # Snakemake Workflow Catalog re-scan)
            if "--lint" not in sys.argv:
                sys.exit(1)

    return {
        "dnaseq": dnaseq,
        "dnaseq_filetype": dnaseq_filetype,
        "dnaseq_readtype": dnaseq_readtype,
        "rnaseq": rnaseq,
        "rnaseq_filetype": rnaseq_filetype,
        "rnaseq_readtype": rnaseq_readtype,
        "normal": normal,
        "custom": {
            "variants": custom_variants,
            "proteins": custom_proteins,
            "hlatyping": {
                "MHC-I": custom_hla_I,
                "MHC-II": custom_hla_II,
            },
        },
    }


def handle_seqfiles(seqdata, mode):
    readtype = []
    filetype = []

    # create new dictionary for modified information
    mod_seqdata = {}

    if seqdata is not None:
        # iterate over replicates

        for rpl in list(seqdata.keys()):
            # make sure to ignore keys with empty values
            if seqdata[rpl] is not None:
                files = [Path(file) for file in seqdata[rpl].split(" ")]
                missing = [str(f) for f in files if not f.is_file()]
                if missing:
                    print(
                        f"[config error] {mode}.{rpl}: file(s) not found: "
                        f"{', '.join(missing)} (paths are case-sensitive on Linux "
                        f"— check for typos).",
                        file=sys.stderr,
                    )
                    if "--lint" not in sys.argv:
                        sys.exit(1)
                elif len(files) == 1:  # SE
                    f1_ext = get_file_extension(files[0])
                    if f1_ext in [".fq", ".fastq", ".bam"]:
                        mod_seqdata[rpl] = files[0]
                        filetype.append(f1_ext)
                        readtype.append("SE")
                    else:
                        print(
                            f"[config error] {mode}.{rpl}: '{files[0]}' is not a valid "
                            f"input file (expected .fq/.fastq/.bam -- this is likely an "
                            f"unfilled placeholder from the default config/config.yaml).",
                            file=sys.stderr,
                        )
                elif len(files) == 2:  # PE
                    f1_ext = get_file_extension(files[0])
                    f2_ext = get_file_extension(files[1])
                    # check if file extensions are the same
                    if f1_ext == f2_ext:
                        # if(valid_paired_end(files[0], files[1])):
                        mod_seqdata[rpl] = files
                        filetype.append(f1_ext)
                        readtype.append("PE")

                        # else:
                        # print('files not in valid PE format')
                    else:
                        ext1 = f1_ext if f1_ext else "?"
                        ext2 = f2_ext if f2_ext else "?"
                        print(
                            f"[config error] {mode}.{rpl}: paired-end files have "
                            f"different extensions ({ext1} vs {ext2}).",
                            file=sys.stderr,
                        )
                        return mod_seqdata, None, None

                # check if filetype and readtype are the same
        #        if all_identical(filetype) and all_identical(readtype):
        if not filetype:
            return mod_seqdata, None, None
        return mod_seqdata, filetype[0], readtype[0]

    else:
        return mod_seqdata, None, None


# determines the file extension for a given file - excludes .gz
def get_file_extension(path):
    filename = path.name
    extpat = r"\.(fastq|fq|bam)(\.gz)?$"
    res = re.search(extpat, filename)
    file_ext = ""
    if res is not None:
        if res.group(0).endswith(".gz"):
            file_ext = filename[res.start() : -3]
        else:
            file_ext = filename[res.start() :]
    return file_ext


# returns the reads (raw/preprocessed) for a given sample
def get_reads(wildcards):
    if config["preproc"]["activate"]:
        if SAMPLES[wildcards.sample][f"{wildcards.readtype}_readtype"] == "SE":
            return SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.replicate]
        elif SAMPLES[wildcards.sample][f"{wildcards.readtype}_readtype"] == "PE":
            return {
                "r1": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
                "r2": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz",
            }


# check if files are a valid paired-end pair
def valid_paired_end(path1, path2):
    valid = False

    # only consider filename
    file1 = path1.name
    file2 = path2.name

    # check if first file contains _R1, _1 or _fwd
    pattern = r"\_(R1|R2|1|2|fwd|rev)\.(fastq|fq){1}(\.gz)?$"
    f1_se = re.search(pattern, file1)
    f2_se = re.search(pattern, file2)

    # patterns needs to be found in both files
    if f1_se is not None and f2_se is not None:
        if file1[: f1_se.start()] == file2[: f2_se.start()]:
            valid = True
        else:
            print("{} and {} have different filestem ".format(file1, file2))
    else:
        print("{} and {} are not valid PE files".format(file1, file2))

    return valid


# check if files in list are identical
def all_identical(l):
    if l.count(l[0]) == len(l):
        return True
    else:
        return False


def print_run_summary(config, samples):
    """Print a one-time, human-readable summary of the resolved run configuration.

    Header carries the shared workflow-level settings (reference, threads,
    variant calling toggles, hlatyping/prioritization classes); a per-sample
    block then shows each sample's input data + custom inputs. Output goes to
    stderr so it sits above the Snakemake job table.
    """
    bar = "=" * 70

    def seq_rows(seqdict, filetype, readtype):
        if not seqdict:
            return ["(none)"]
        ft = filetype if filetype else "?"
        rt = readtype if readtype else "?"
        rows = []
        for rpl, path in seqdict.items():
            if isinstance(path, list):
                loc = " , ".join(str(p) for p in path)
            else:
                loc = str(path)
            rows.append(f"{rpl}  [{ft}, {rt}]  {loc}")
        return rows

    def mode_state(active, detail):
        state = "on" if active else "off"
        if detail:
            return f"{state}  ({detail})"
        return state

    # resolve workflow-level display values up front; the f-strings below only
    # ever substitute a plain variable (Snakemake's .smk parser mishandles
    # f-string replacement fields containing operators or nested quotes)
    release = str(config["reference"]["release"])
    threads = str(config["threads"])
    indel_type = str(config["indel"]["type"])
    indel_mode = str(config["indel"]["mode"])
    indel_detail = f"type {indel_type}, mode {indel_mode}"
    hla_class = str(config["hlatyping"]["class"])
    mhc1_mode = str(config["hlatyping"]["MHC-I_mode"])
    mhc2_mode = str(config["hlatyping"]["MHC-II_mode"])
    pri_class = str(config["prioritization"]["class"])
    len1 = str(config["prioritization"]["lengths"]["MHC-I"])
    len2 = str(config["prioritization"]["lengths"]["MHC-II"])

    n_samples = len(samples)
    lines = [
        bar,
        f"  ScanNeo2 -- {n_samples} sample(s)",
        bar,
        f"  Reference        : Ensembl release {release}",
        f"  Threads per rule : {threads}",
        "",
        "  Variant calling",
    ]
    variant_modes = [
        ("Indels", config["indel"]["activate"], indel_detail),
        ("Alt. splicing", config["altsplicing"]["activate"], ""),
        ("Exitrons", config["exitronsplicing"]["activate"], ""),
        ("Gene fusion", config["genefusion"]["activate"], ""),
    ]
    for label, active, detail in variant_modes:
        pad = label.ljust(14)
        state = mode_state(active, detail)
        lines.append(f"    {pad}: {state}")
    lines.append("")
    hla_parts = []
    if hla_class in ("I", "BOTH"):
        hla_parts.append(f"MHC-I: {mhc1_mode}")
    if hla_class in ("II", "BOTH"):
        hla_parts.append(f"MHC-II: {mhc2_mode}")
    lines.append(f"  HLA typing       : class {hla_class}  ({' | '.join(hla_parts)})")
    pri_parts = []
    if pri_class in ("I", "BOTH"):
        pri_parts.append(f"MHC-I: {len1}")
    if pri_class in ("II", "BOTH"):
        pri_parts.append(f"MHC-II: {len2}")
    lines.append(
        f"  Prioritization   : class {pri_class}  (epitope lengths -- {' | '.join(pri_parts)})"
    )
    lines.append(bar)
    lines.append("  Samples")

    for sample_name, d in samples.items():
        lines.append("")
        lines.append(f"    {sample_name}  (results/{sample_name}/)")
        normal_disp = str(d["normal"]) if d["normal"] else "(none)"
        custom_disp = (
            str(d["custom"]["variants"])
            if d["custom"]["variants"]
            else "(no custom variants)"
        )
        proteins_disp = (
            str(d["custom"]["proteins"])
            if d["custom"]["proteins"]
            else "(no custom proteins)"
        )
        for label, key in (("DNAseq", "dnaseq"), ("RNAseq", "rnaseq")):
            ft_key = f"{key}_filetype"
            rt_key = f"{key}_readtype"
            rows = seq_rows(d[key], d.get(ft_key), d.get(rt_key))
            pad = label.ljust(9)
            first = rows[0]
            lines.append(f"      {pad}: {first}")
            for extra in rows[1:]:
                lines.append(f"               {extra}")
        npad = "Normal".ljust(9)
        cpad = "Custom".ljust(9)
        ppad = "Proteins".ljust(9)
        lines.append(f"      {npad}: {normal_disp}")
        lines.append(f"      {cpad}: {custom_disp}")
        lines.append(f"      {ppad}: {proteins_disp}")
        for cls in ("MHC-I", "MHC-II"):
            custom_alleles = d["custom"]["hlatyping"].get(cls)
            if custom_alleles:
                ca = str(custom_alleles)
                lines.append(f"               custom {cls} alleles: {ca}")
    lines.append(bar)

    print("\n".join(lines), file=sys.stderr)


def check_vendored_scripts(config):
    """Surface missing vendored entry-point scripts at workflow-load time so
    the user gets a clear actionable message instead of a cryptic mid-run
    failure deep inside indel / exitron / prioritization rules.

    Two failure modes are covered:
      (a) git submodule directory present but empty (clone without
          `--recurse-submodules`) — applies to `workflow/scripts/transindel/`
          and `workflow/scripts/scanexitron/`. Gated on the feature using
          the submodule being active so a user who has indel/exitron
          disabled isn't blocked.
      (b) IEDB tool directory present but its entry-point script missing
          (partial download). Snakemake's missing-output detection looks at
          the directory only, so a half-populated dir doesn't re-trigger the
          download rule. We only flag when the parent dir exists; if it's
          absent the download rule will run cleanly.
    """
    missing = []

    if config["indel"]["activate"] and config["indel"]["type"] in ("long", "all"):
        p = Path("workflow/scripts/transindel/transIndel_build_DNA.py")
        if not p.is_file():
            missing.append(
                (
                    str(p),
                    "git submodule not initialised — "
                    "run `git submodule update --init --recursive`",
                )
            )

    if config["exitronsplicing"]["activate"]:
        p = Path("workflow/scripts/scanexitron/ScanExitron.py")
        if not p.is_file():
            missing.append(
                (
                    str(p),
                    "git submodule not initialised — "
                    "run `git submodule update --init --recursive`",
                )
            )

    iedb_targets = [
        Path("workflow/scripts/mhc_i/src/predict_binding.py"),
        Path("workflow/scripts/mhc_ii/mhc_II_binding.py"),
        Path("workflow/scripts/immunogenicity/predict_immunogenicity.py"),
    ]
    for p in iedb_targets:
        if p.parent.is_dir() and not p.is_file():
            missing.append(
                (
                    str(p),
                    f"IEDB tool directory exists but is incomplete — "
                    f"remove `{p.parent}` to let Snakemake re-trigger the download rule",
                )
            )

    if not missing:
        return

    print("[config error] required vendored script(s) missing:", file=sys.stderr)
    for path, hint in missing:
        print(f"  {path}  ({hint})", file=sys.stderr)
    # skip the abort under `snakemake --lint` so static rule analysis still
    # completes on a fresh clone where the submodules / IEDB dirs are absent
    # (e.g. for the Snakemake Workflow Catalog re-scan).
    if "--lint" not in sys.argv:
        sys.exit(1)


def check_hlahd_setup(config):
    """Surface missing HLA-HD prerequisites at workflow-load time so the user
    gets a clear actionable message instead of a mid-run failure deep inside
    the class-II HLA typing rules.

    HLA-HD itself is installed by the user (not by conda — `hlahd.yml` only
    carries bowtie2, HLA-HD's runtime dependency). The reference data files
    are likewise downloaded out-of-band from the HLA-HD distribution.

    All checks are gated on class-II typing being active via DNA or RNA
    reads — if `MHC-II_mode` is purely `custom` (alleles supplied directly
    by the user), HLA-HD isn't invoked and we don't block.
    """
    if config["hlatyping"]["class"] not in ("II", "BOTH"):
        return
    mhc2_mode = config["hlatyping"]["MHC-II_mode"]
    if "DNA" not in mhc2_mode and "RNA" not in mhc2_mode:
        return

    errors = []

    if shutil.which("hlahd.sh") is None:
        errors.append(
            "hlahd.sh not found on PATH — HLA-HD is not pulled in by "
            "`hlahd.yml` (which only carries bowtie2). Install HLA-HD from "
            "https://www.genome.med.kyoto-u.ac.jp/HLA-HD/ and ensure "
            "`hlahd.sh` is on PATH."
        )

    freqdata = Path(config["hlatyping"]["freqdata"])
    if not freqdata.is_dir():
        errors.append(
            f"hlatyping.freqdata: directory not found: '{freqdata}' — "
            "expected the HLA-HD `freq_data/` directory shipped with the tool."
        )

    split = Path(config["hlatyping"]["split"])
    if not split.is_file():
        errors.append(
            f"hlatyping.split: file not found: '{split}' — "
            "expected HLA-HD's `HLA_gene.split.txt`."
        )

    dict_dir = Path(config["hlatyping"]["dict"])
    if not dict_dir.is_dir():
        errors.append(
            f"hlatyping.dict: directory not found: '{dict_dir}' — "
            "expected the HLA-HD `dictionary/` directory shipped with the tool."
        )

    if not errors:
        return

    print(
        "[config error] HLA-HD prerequisites missing (class II hlatyping "
        "is active via DNA/RNA mode):",
        file=sys.stderr,
    )
    for err in errors:
        print(f"  {err}", file=sys.stderr)
    # skip the abort under `snakemake --lint` so static rule analysis still
    # completes on a fresh clone where HLA-HD isn't installed (e.g. for the
    # Snakemake Workflow Catalog re-scan).
    if "--lint" not in sys.argv:
        sys.exit(1)


def check_cross_field_consistency(config, samples):
    """Surface cross-field config inconsistencies that would otherwise produce
    empty / nonsensical results (no crash) or a mid-run failure inside an HLA
    typing rule. Both cases waste a full pipeline run before the user notices.

    Two checks:
      1. `prioritization.class` must be a subset of `hlatyping.class` — every
         MHC class the user wants to prioritize must also be typed. Otherwise
         the prioritization step finishes with no candidate alleles for the
         missing class and the user sees nothing to look at. Workflow-level
         (sample-independent).
      2. If `hlatyping.MHC-{I,II}_mode` contains `custom`, each sample must
         provide the matching `custom_hla_{I,II}` cell. The per-sample
         `custom_paths` check in `per_sample_data` already validates "if a
         path is given, it exists"; this one catches the complementary case
         of "custom mode promised but no path given". Gated on the class
         actually covering the relevant MHC.
    """
    errors = []

    p_class = config["prioritization"]["class"]
    h_class = config["hlatyping"]["class"]

    requires_I = p_class in ("I", "BOTH")
    requires_II = p_class in ("II", "BOTH")
    types_I = h_class in ("I", "BOTH")
    types_II = h_class in ("II", "BOTH")

    if requires_I and not types_I:
        errors.append(
            f"prioritization.class='{p_class}' asks for MHC-I neoepitopes "
            f"but hlatyping.class='{h_class}' does not type MHC-I — set "
            "hlatyping.class to 'I' or 'BOTH', or drop MHC-I from prioritization."
        )
    if requires_II and not types_II:
        errors.append(
            f"prioritization.class='{p_class}' asks for MHC-II neoepitopes "
            f"but hlatyping.class='{h_class}' does not type MHC-II — set "
            "hlatyping.class to 'II' or 'BOTH', or drop MHC-II from prioritization."
        )

    custom_I = types_I and "custom" in config["hlatyping"]["MHC-I_mode"]
    custom_II = types_II and "custom" in config["hlatyping"]["MHC-II_mode"]
    for sample_name, d in samples.items():
        if custom_I and d["custom"]["hlatyping"].get("MHC-I") is None:
            errors.append(
                f"sample {sample_name!r}: hlatyping.MHC-I_mode contains 'custom' "
                "but the custom_hla_I cell is empty — provide a path to a TSV of "
                "MHC-I alleles for this sample, or drop 'custom' from MHC-I_mode."
            )
        if custom_II and d["custom"]["hlatyping"].get("MHC-II") is None:
            errors.append(
                f"sample {sample_name!r}: hlatyping.MHC-II_mode contains 'custom' "
                "but the custom_hla_II cell is empty — provide a path to a TSV of "
                "MHC-II alleles for this sample, or drop 'custom' from MHC-II_mode."
            )

    if not errors:
        return

    print("[config error] cross-field inconsistency:", file=sys.stderr)
    for err in errors:
        print(f"  {err}", file=sys.stderr)
    # skip the abort under `snakemake --lint` so static rule analysis still
    # completes on a fresh / placeholder config (e.g. for the Snakemake
    # Workflow Catalog re-scan).
    if "--lint" not in sys.argv:
        sys.exit(1)


# build the SAMPLES map from the sheet loaded in Snakefile and run the
# workflow-load validators (issue #93)
if len(SHEET) == 0:
    print(
        f"[config error] sample sheet {config['samples']!r} has no rows; add at "
        "least one sample row before running.",
        file=sys.stderr,
    )
    if "--lint" not in sys.argv:
        sys.exit(1)

dup_samples = SHEET["sample"][SHEET["sample"].duplicated()].tolist()
if dup_samples:
    print(
        f"[config error] sample sheet has duplicate sample name(s): "
        f"{sorted(set(dup_samples))}. Each row must have a unique 'sample' value.",
        file=sys.stderr,
    )
    if "--lint" not in sys.argv:
        sys.exit(1)

SAMPLES = {row["sample"]: per_sample_data(row) for _, row in SHEET.iterrows()}
check_vendored_scripts(config)
check_hlahd_setup(config)
check_cross_field_consistency(config, SAMPLES)
print_run_summary(config, SAMPLES)


########### PREPROCESSING ##########
def get_raw_reads(wildcards):
    if SAMPLES[wildcards.sample]["rnaseq_readtype"] == "SE":
        return dict(zip(["sample"], SAMPLES[wildcards.sample]["rnaseq"][wildcards.group]))


# returns raw reads for a given sample
def get_qc_input(wildcards):
    return SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group]


def get_qc_input_fwd(wildcards):
    return SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group][0]


def get_qc_input_rev(wildcards):
    return SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group][1]


# returns the reads (raw/preprocessed) for a given sample
def get_preproc_input(wildcards):
    if config["preproc"]["activate"]:
        if SAMPLES[wildcards.sample][f"{wildcards.seqtype}_readtype"] == "SE":
            return {SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group]}

        elif SAMPLES[wildcards.sample][f"{wildcards.seqtype}_readtype"] == "PE":
            return {
                "sample": [
                    SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group][0],
                    SAMPLES[wildcards.sample][wildcards.seqtype][wildcards.group][1],
                ]
            }


########### HLA GENOTYPING ##########
def get_input_reads_hlatyping_BAM(wildcards):
    seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
    return SAMPLES[wildcards.sample][seqtype][wildcards.group]


def get_input_filtering_hlatyping_SE(wildcards):
    seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
    if SAMPLES[wildcards.sample][f"{seqtype}_filetype"] == ".bam":
        return expand(
            "results/{sample}/{seqtype}/reads/{group}_flt_BAM.fq",
            sample=wildcards.sample,
            seqtype=seqtype,
            group=wildcards.group,
        )
    else:
        if config["preproc"]["activate"]:
            return expand(
                "results/{sample}/{seqtype}/reads/{group}_preproc.fq.gz",
                sample=wildcards.sample,
                seqtype=seqtype,
                group=wildcards.group,
            )
        else:
            return SAMPLES[wildcards.sample][seqtype][wildcards.group]


def get_input_filtering_hlatyping_PE(wildcards):
    if config["preproc"]["activate"]:
        return expand(
            "results/{sample}/{seqtype}/reads/{group}_{readpair}_preproc.fq.gz",
            sample=wildcards.sample,
            seqtype="dnaseq" if wildcards.nartype == "DNA" else "rnaseq",
            group=wildcards.group,
            nartype=wildcards.nartype,
            readpair=wildcards.readpair,
        )
    else:
        seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
        return SAMPLES[wildcards.sample][f"{wildcards.seqtype}"][wildcards.group]


def aggregate_mhcI_SE(wildcards):
    checkpoint_output = checkpoints.split_reads_mhcI_SE.get(**wildcards).output[0]
    return expand(
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE/{no}_result.tsv",
        sample=wildcards.sample,
        group=wildcards.group,
        nartype=wildcards.nartype,
        no=glob_wildcards(os.path.join(checkpoint_output, "R_{no}.bam")).no,
    )


def aggregate_mhcI_PE(wildcards):
    checkpoint_output = checkpoints.split_reads_mhcI_PE.get(**wildcards).output[0]
    return expand(
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE/{no}_result.tsv",
        sample=wildcards.sample,
        group=wildcards.group,
        nartype=wildcards.nartype,
        no=glob_wildcards(os.path.join(checkpoint_output, "R1_{no}.bam")).no,
    )


def get_all_mhcI_alleles(wildcards):
    values = []

    if "DNA" in config["hlatyping"]["MHC-I_mode"]:
        if len(SAMPLES[wildcards.sample]["dnaseq"]) != 0:
            if (
                SAMPLES[wildcards.sample]["dnaseq_readtype"] == "SE"
                or SAMPLES[wildcards.sample]["dnaseq_filetype"] == ".bam"
            ):
                for key in SAMPLES[wildcards.sample]["dnaseq"].keys():
                    values += expand(
                        "results/{sample}/hla/mhc-I/genotyping/{group}_DNA_flt_SE.tsv",
                        sample=wildcards.sample,
                        group=key,
                    )
            elif SAMPLES[wildcards.sample]["dnaseq_readtype"] == "PE":
                for key in SAMPLES[wildcards.sample]["dnaseq"].keys():
                    values += expand(
                        "results/{sample}/hla/mhc-I/genotyping/{group}_DNA_flt_PE.tsv",
                        sample=wildcards.sample,
                        group=key,
                    )

    if "RNA" in config["hlatyping"]["MHC-I_mode"]:
        if len(SAMPLES[wildcards.sample]["rnaseq"]) != 0:
            if (
                SAMPLES[wildcards.sample]["rnaseq_readtype"] == "SE"
                or SAMPLES[wildcards.sample]["rnaseq_filetype"] == ".bam"
            ):
                for key in SAMPLES[wildcards.sample]["rnaseq"].keys():
                    values += expand(
                        "results/{sample}/hla/mhc-I/genotyping/{group}_RNA_flt_SE.tsv",
                        sample=wildcards.sample,
                        group=key,
                    )
            elif SAMPLES[wildcards.sample]["rnaseq_readtype"] == "PE":
                for key in SAMPLES[wildcards.sample]["rnaseq"].keys():
                    values += expand(
                        "results/{sample}/hla/mhc-I/genotyping/{group}_RNA_flt_PE.tsv",
                        sample=wildcards.sample,
                        group=key,
                    )

    if "custom" in config["hlatyping"]["MHC-I_mode"]:
        values += [SAMPLES[wildcards.sample]["custom"]["hlatyping"]["MHC-I"]]

    if len(values) == 0:
        print(
            "No HLA data found. Please check the config file for correct specification of data and HLA genotyping mode"
        )
        sys.exit(1)

    return values


##### MHC CLASS II #####
def get_input_filter_reads_mhcII_SE(wildcards):
    if config["preproc"]["activate"]:
        return expand(
            "results/{sample}/{seqtype}/reads/{group}_preproc.fq.gz",
            sample=wildcards.sample,
            seqtype="dnaseq" if wildcards.nartype == "DNA" else "rnaseq",
            group=wildcards.group,
        )
    else:
        seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
        return SAMPLES[wildcards.sample][f"{seqtype}"][wildcards.group]


def get_input_filter_reads_mhcII_PE(wildcards):
    seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
    """ if the reads are in BAM format, we consider them as processed and forward
  them to the filtering rule. If the reads are in FASTQ format, we may preprocess"""
    if SAMPLES[wildcards.sample][f"{seqtype}_filetype"] == ".bam":
        return expand(
            "results/{sample}/hla/reads/{group}_{nartype}_BAM.fq",
            sample=wildcards.sample,
            group=wildcards.group,
            nartype=wildcards.nartype,
        )

    else:
        if config["preproc"]["activate"]:
            return expand(
                "results/{sample}/{seqtype}/reads/{group}_{readpair}_preproc.fq.gz",
                sample=wildcards.sample,
                seqtype="dnaseq" if wildcards.nartype == "DNA" else "rnaseq",
                group=wildcards.group,
                readpair=["R1", "R2"],
            )
        else:
            return SAMPLES[wildcards.sample][f"{seqtype}"][wildcards.group]


def get_output_hlatyping_mhcII(wildcards):
    seqtype = "dnaseq" if wildcards.nartype == "DNA" else "rnaseq"
    if (
        SAMPLES[wildcards.sample][f"{seqtype}_filetype"] == ".bam"
        or SAMPLES[wildcards.sample][f"{seqtype}_readtype"] == "PE"
    ):
        # we consider the reads in BAM files as PE (so we can use the rule filter_reads_mhcII_PE)
        return expand(
            "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_PE.bam",
            sample=wildcards.sample,
            group=wildcards.group,
            nartype=wildcards.nartype,
        )
    else:
        if SAMPLES[wildcards.sample][f"{seqtype}_readtype"] == "SE":
            return expand(
                "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_SE.fq",
                sample=wildcards.sample,
                group=wildcards.group,
                nartype=wildcards.nartype,
            )


def get_predicted_mhcII_alleles(wildcards):
    values = []

    # routines to genotype from DNA
    if "DNA" in config["hlatyping"]["MHC-II_mode"]:
        if SAMPLES[wildcards.sample]["dnaseq"] is not None:
            for key in SAMPLES[wildcards.sample]["dnaseq"].keys():

                if SAMPLES[wildcards.sample]["normal"] is not None:
                    if key in SAMPLES[wildcards.sample]["normal"]:
                        continue

                values += expand(
                    "results/{sample}/hla/mhc-II/genotyping/{group}_{nartype}/result/{group}_{nartype}_final.result.txt",
                    sample=wildcards.sample,
                    group=key,
                    nartype="DNA",
                )

        else:  # if no dnaseq data is specified, but mode is DNA or BOTH, then ignore
            print(
                "dnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is DNA or BOTH -- will be ignored"
            )

    # routines to genotype from RNA
    if "RNA" in config["hlatyping"]["MHC-II_mode"]:
        if SAMPLES[wildcards.sample]["rnaseq"] is not None:
            for key in SAMPLES[wildcards.sample]["rnaseq"].keys():

                # exclude normal samples (if specified)
                if SAMPLES[wildcards.sample]["normal"] is not None:
                    normal = SAMPLES[wildcards.sample]["normal"].split(" ")
                    if key in normal:
                        continue

                values += expand(
                    "results/{sample}/hla/mhc-II/genotyping/{group}_{nartype}/result/{group}_{nartype}_final.result.txt",
                    sample=wildcards.sample,
                    group=key,
                    nartype="RNA",
                )

        else:  # if no rnaseq data is specified, but mode is RNA or BOTH, then ignore
            print(
                "rnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is RNA or BOTH -- will be ignored"
            )

    if len(values) == 0:
        print(
            "No data found. Check config file for correct specification of data and hla genotyping mode"
        )
        sys.exit(1)

    return values


def get_all_mhcII_alleles(wildcards):
    values = []

    if (
        "DNA" in config["hlatyping"]["MHC-II_mode"]
        or "RNA" in config["hlatyping"]["MHC-II_mode"]
    ):
        values += expand(
            "results/{sample}/hla/mhc-II/genotyping/mhc-II.tsv", sample=wildcards.sample
        )

    if "custom" in config["hlatyping"]["MHC-II_mode"]:
        if SAMPLES[wildcards.sample]["custom"]["hlatyping"]["MHC-II"] is not None:
            values += [SAMPLES[wildcards.sample]["custom"]["hlatyping"]["MHC-II"]]
        else:
            print("No custom alleles specified in config file for MHC-II hlatyping")
            sys.exit(1)

    if len(values) == 0:
        print(
            "No hla data found. Check config file for correct specification of data and hla genotyping mode"
        )
        sys.exit(1)

        print(values)

    return values


########### ALIGNMENT ##########
def get_rnaseq_star_bam(wildcards):
    """Return the STAR-aligned BAM path for a sample's rnaseq filetype.

    Lets star_align_fastq and merge_alignment_results coexist as
    always-defined rules with distinct outputs (issue #93). The downstream
    rnaseq_postproc_fixmate rule reads through this helper instead of a
    hard-coded path; the helper consults SAMPLES[sample]['rnaseq_filetype']
    and returns the corresponding rule's output.
    """
    ft = SAMPLES[wildcards.sample]["rnaseq_filetype"]
    if ft == ".bam":
        return (
            f"results/{wildcards.sample}/rnaseq/align/"
            f"{wildcards.group}_aligned_STAR.from_bam.bam"
        )
    return (
        f"results/{wildcards.sample}/rnaseq/align/"
        f"{wildcards.group}_aligned_STAR.bam"
    )


def get_star_input(wildcards):
    if SAMPLES[wildcards.sample]["rnaseq_filetype"] == ".bam":
        return dict(zip(["bam"], [SAMPLES[wildcards.sample]["rnaseq"][wildcards.group]]))

    elif (
        SAMPLES[wildcards.sample]["rnaseq_filetype"] == ".fq"
        or SAMPLES[wildcards.sample]["rnaseq_filetype"] == ".fastq"
    ):
        if config["preproc"]["activate"]:
            if SAMPLES[wildcards.sample]["rnaseq_readtype"] == "SE":
                return dict(
                    zip(
                        ["fq1"],
                        expand(
                            "results/{sample}/rnaseq/reads/{group}_preproc.fq.gz",
                            sample=wildcards.sample,
                            group=wildcards.group,
                        ),
                    )
                )
            elif SAMPLES[wildcards.sample]["rnaseq_readtype"] == "PE":  # PE
                return dict(
                    zip(
                        ["fq1", "fq2"],
                        expand(
                            "results/{sample}/rnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                            readtype=["R1", "R2"],
                            group=wildcards.group,
                            sample=wildcards.sample,
                        ),
                    )
                )
            else:  # no pre-processing
                return SAMPLES[wildcards.sample]["rnaseq"][wildcards.group]


# collect the individual alignments from splitted bamfiles
def aggregate_aligned_rg(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.split_bamfile_RG.get(**wildcards).output[0]
    return expand(
        "results/{sample}/rnaseq/align/{group}/{rg}.bam",
        sample=wildcards.sample,
        group=wildcards.group,
        rg=glob_wildcards(os.path.join(checkpoint_output, "{rg}.bam")).rg,
    )


def get_readgroups_input(wildcards):
    # return only bam from STAR align
    if SAMPLES[wildcards.sample][f"{wildcards.seqtype}_filetype"] in [".fq", ".fastq"]:
        return expand(
            "results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
            sample=wildcards.sample,
            seqtype=wildcards.seqtype,
            group=wildcards.group,
        )

    #    return ["results/{sample}/{seqtype}/align/{group}_final_STAR.bam".format(**wildcards)]

    elif SAMPLES[wildcards.sample][f"{wildcards.seqtype}_filetype"] in [".bam"]:
        val = []

        # For DNAseq
        if wildcards.seqtype == "dnaseq":
            val.append(str(SAMPLES[wildcards.sample]["dnaseq"][wildcards.group]))
        elif wildcards.seqtype == "rnaseq":
            # needs both the raw data and star aligned bam
            val.append(str(SAMPLES[wildcards.sample]["rnaseq"][wildcards.group]))
            val += expand(
                "results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
                sample=wildcards.sample,
                seqtype="rnaseq",
                group=wildcards.group,
            )

        return val


# input for alignment w/ DNAseq data
def get_realign_input(wildcards):
    # For DNAseq use the (raw) reads defined in config
    if wildcards.seqtype == "dnaseq":
        return SAMPLES[wildcards.sample]["dnaseq"][wildcards.group]
    elif wildcards.seqtype == "rnaseq":
        return expand(
            "results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
            sample=wildcards.sample,
            seqtype="rnaseq",
            group=wildcards.group,
        )


def get_align_input_dnaseq(wildcards):
    if config["preproc"]["activate"]:
        if SAMPLES[wildcards.sample]["dnaseq_readtype"] == "SE":
            return expand("results/{sample}/dnaseq/reads/inputreads.fq.gz", **wildcards)
        else:  # PE
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/{sample}/dnaseq/reads/inputreads_{readtype}.fq.gz",
                        readtype=["r1", "r2"],
                        **wildcards,
                    ),
                )
            )
    else:  # no pre-processing
        if SAMPLES[wildcards.sample]["dnaseq_readtype"] == "SE":
            return dnaseq_input[wildcards.sample]
        else:  # PE
            return dict(zip(["fq1", "fq2"], dnaseq_input[wildcards.sample]))


def get_dna_align_input(wildcards):
    if SAMPLES[wildcards.sample]["dnaseq_filetype"] in [".fq", ".fastq"]:
        if config["preproc"]["activate"]:
            if SAMPLES[wildcards.sample]["dnaseq_readtype"] == "SE":
                return expand(
                    "results/{sample}/dnaseq/reads/{group}_preproc.fq.gz", **wildcards
                )
            elif SAMPLES[wildcards.sample]["dnaseq_readtype"] == "PE":  # PE
                return expand(
                    "results/{sample}/dnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                    readtype=["R1", "R2"],
                    group=wildcards.group,
                    sample=wildcards.sample,
                )
        else:  # no pre-processing
            return SAMPLES[wildcards.sample]["dnaseq"][wildcards.group]

    # is bamfile
    else:  # no pre-processing has been performed
        return SAMPLES[wildcards.sample]["dnaseq"][wildcards.group]


########### GENE EXPRESSION ##########
def get_aligned_reads_featurecounts(wildcards):
    if wildcards.seqtype == "dnaseq":
        if SAMPLES[wildcards.sample]["dnaseq_filetype"] == ".bam":
            return SAMPLES[wildcards.sample]["dnaseq"][wildcards.group]
        else:
            return expand(
                "results/{sample}/{seqtype}/align/{group}_aligned_BWA.bam",
                sample=wildcards.sample,
                seqtype="dnaseq",
                group=wildcards.group,
            )

    elif wildcards.seqtype == "rnaseq":
        if SAMPLES[wildcards.sample]["rnaseq_filetype"] == ".bam":
            return SAMPLES[wildcards.sample]["rnaseq"][wildcards.group]
        else:
            return expand(
                "results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
                sample=wildcards.sample,
                seqtype="rnaseq",
                group=wildcards.group,
            )


def get_counts(wildcards):
    counts = []

    if SAMPLES[wildcards.sample]["dnaseq"] is not None:
        if config["quantification"]["mode"] in ["DNA", "BOTH"]:
            counts += expand(
                "results/{sample}/{seqtype}/quantification/{group}_counts.txt",
                sample=wildcards.sample,
                seqtype="dnaseq",
                group=list(SAMPLES[wildcards.sample]["dnaseq"].keys()),
            )

    else:
        message = (
            f"dnaseq data has not been specified in the config file"
            f"but specified mode for counts in config file is RNA or BOTH"
            f" -- will be ignored"
        )
        print(message)

    if SAMPLES[wildcards.sample]["rnaseq"] is not None:
        if config["quantification"]["mode"] in ["RNA", "BOTH"]:
            counts += expand(
                "results/{sample}/{seqtype}/quantification/{group}_counts.txt",
                sample=wildcards.sample,
                seqtype="rnaseq",
                group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
            )
    else:
        message = (
            f"rnaseq data has not been specified in the config file"
            f"but specified mode for counts in config file is DNA or BOTH"
            f" -- will be ignored"
        )
        print(message)

    if len(counts) == 0:
        print("No counts will be generated!")

    return counts


########### CUSTOM VARIANTS ##########
def get_custom_variants(wildcards):
    return SAMPLES[wildcards.sample]["custom"]["variants"]


########### INDEL CALLING ##########
def get_longindels(wildcards):
    indels = []
    if SAMPLES[wildcards.sample]["dnaseq"] is not None:
        if config["indel"]["mode"] in ["DNA", "BOTH"]:
            indels += expand(
                "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz",
                sample=wildcards.sample,
                seqtype="dnaseq",
                group=list(SAMPLES[wildcards.sample]["dnaseq"].keys()),
            )

    if SAMPLES[wildcards.sample]["rnaseq"] is not None:
        if config["indel"]["mode"] in ["RNA", "BOTH"]:
            indels += expand(
                "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz",
                sample=wildcards.sample,
                seqtype="rnaseq",
                group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
            )

    return indels


# htc first round collect/combine results
def aggregate_vcf_htc_first_round(wildcards):
    checkpoint_output = checkpoints.split_bam_htc_first_round.get(**wildcards).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd/{chr}.vcf.gz",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


def aggregate_idx_htc_first_round(wildcards):
    checkpoint_output = checkpoints.split_bam_htc_first_round.get(**wildcards).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd/{chr}.vcf.gz.tbi",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


# htc final round collect/combine results
def aggregate_vcf_htc_final_round(wildcards):
    checkpoint_output = checkpoints.split_bam_htc_final_round.get(**wildcards).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final/{chr}.vcf.gz",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


def aggregate_idx_htc_final_round(wildcards):
    checkpoint_output = checkpoints.split_bam_htc_final_round.get(**wildcards).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final/{chr}.vcf.gz.tbi",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


# mutect2 collect/combine results
def aggregate_vcf_mutect2(wildcards):
    checkpoint_output = checkpoints.split_bam_detect_short_indels_m2.get(
        **wildcards
    ).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


def aggregate_idx_mutect2(wildcards):
    checkpoint_output = checkpoints.split_bam_detect_short_indels_m2.get(
        **wildcards
    ).output[0]
    return expand(
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz.tbi",
        sample=wildcards.sample,
        seqtype=wildcards.seqtype,
        group=wildcards.group,
        chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.bam")).chr,
    )


def get_shortindels(wildcards):
    indels = []

    if config["indel"]["mode"] in ["DNA", "BOTH"]:
        if SAMPLES[wildcards.sample]["dnaseq"] is not None:
            indels += expand(
                "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
                sample=wildcards.sample,
                seqtype="dnaseq",
                group=list(SAMPLES[wildcards.sample]["dnaseq"].keys()),
            )
        else:
            print(
                "dnaseq data has not been specified in the config file, but specified mode for indel calling in config file is DNA or BOTH - skipping..."
            )

    if config["indel"]["mode"] in ["RNA", "BOTH"]:
        if SAMPLES[wildcards.sample]["rnaseq"] is not None:
            indels += expand(
                "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
                sample=wildcards.sample,
                seqtype="rnaseq",
                group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
            )
        else:
            print(
                "rnaseq data has not been specified in the config file, but specified mode for indel calling in config file is RNA or BOTH"
            )

    if len(indels) == 0:
        print(
            "No data found. Check config file for correct specification of data and indel calling mode"
        )
        sys.exit(1)

    return indels


def get_snvs(wildcards):
    snvs = []
    if config["indel"]["mode"] in ["RNA", "BOTH"]:
        snvs += expand(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz",
            sample=wildcards.sample,
            seqtype="rnaseq",
            group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
        )

    if config["indel"]["mode"] in ["DNA", "BOTH"]:
        snvs += expand(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz",
            sample=wildcards.sample,
            seqtype="dnaseq",
            group=list(SAMPLES[wildcards.sample]["dnaseq"].keys()),
        )

    return snvs


########### EXITRON CALLING ##########
def get_exitrons(wildcards):
    return expand(
        "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf.gz",
        sample=wildcards.sample,
        group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
    )


########### GENE FUSIONS ##########
def get_fusions(wildcards):
    fusions = []
    if SAMPLES[wildcards.sample]["rnaseq"] is not None:
        if config["genefusion"]["activate"]:
            fusions += expand(
                "results/{sample}/rnaseq/genefusion/{group}_fusions.tsv",
                sample=wildcards.sample,
                group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
            )

    return fusions


########### ALT SPLICING ##########
def get_altsplicing(wildcards):
    altsplicing = []
    if SAMPLES[wildcards.sample]["rnaseq"] is not None:
        if config["altsplicing"]["activate"]:
            altsplicing += expand(
                "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf.gz",
                sample=wildcards.sample,
                group=list(SAMPLES[wildcards.sample]["rnaseq"].keys()),
            )
    return altsplicing


########### NEOANTIGEN PRIORIZATION ##########
def get_prioritization_snvs(wildcards):
    snv = []
    if config["indel"]["activate"]:
        if config["indel"]["type"] in ["short", "all"]:
            snv += expand(
                "results/{sample}/annotation/somatic.snvs.vcf",
                sample=wildcards.sample,
            )

    return snv


def get_prioritization_indels(wildcards):
    indels = []
    if config["indel"]["activate"]:
        if config["indel"]["type"] in ["short", "all"]:
            indels += expand(
                "results/{sample}/annotation/somatic.short.indels.vcf",
                sample=wildcards.sample,
            )

    return indels


def get_prioritization_long_indels(wildcards):
    long_indels = []
    if config["indel"]["activate"]:
        if config["indel"]["type"] in ["long", "all"]:
            long_indels += expand(
                "results/{sample}/annotation/long.indels.vcf",
                sample=wildcards.sample,
            )

    return long_indels


def get_prioritization_exitrons(wildcards):
    exitrons = []
    if config["exitronsplicing"]["activate"]:
        if len(SAMPLES[wildcards.sample]["rnaseq"]) != 0:
            exitrons += expand(
                "results/{sample}/annotation/exitrons.vcf",
                sample=wildcards.sample,
            )
        else:
            print(
                "rnaseq data has not been specified in the config file, but exitron calling is activated - skipping..."
            )
    return exitrons


def get_prioritization_altsplicing(wildcards):
    altsplicing = []
    if config["altsplicing"]["activate"]:
        if len(SAMPLES[wildcards.sample]["rnaseq"]) != 0:
            altsplicing += expand(
                "results/{sample}/annotation/altsplicing.vcf",
                sample=wildcards.sample,
            )
        else:
            print(
                "rnaseq data has not been specified in the config file, but alternative splicing calling is activated - skipping..."
            )
    return altsplicing


def get_prioritization_custom(wildcards):
    custom = []
    if SAMPLES[wildcards.sample]["custom"]["variants"] is not None:
        custom += expand(
            "results/{sample}/annotation/custom.vcf", sample=wildcards.sample
        )

    return custom


def get_prioritization_proteins(wildcards):
    """User-supplied TSV of (wildtype, mutant) protein pairs, consumed
    verbatim — no upstream rule, no VEP annotation."""
    if SAMPLES[wildcards.sample]["custom"]["proteins"] is None:
        return []
    return [SAMPLES[wildcards.sample]["custom"]["proteins"]]


def get_prioritization_mhcI(wildcards):
    alleles = []
    if config["prioritization"]["class"] in ["I", "BOTH"]:
        alleles += expand(
            "results/{sample}/hla/mhc-I.tsv", sample=wildcards.sample
        )

    return alleles


def get_prioritization_mhcII(wildcards):
    alleles = []
    if config["prioritization"]["class"] in ["II", "BOTH"]:
        alleles += expand(
            "results/{sample}/hla/mhc-II.tsv", sample=wildcards.sample
        )
    return alleles


def get_prioritization_counts(wildcards):
    counts = []
    # counts can only be generated if either RNAseq or DNAseq data is provided
    if len(SAMPLES[wildcards.sample]["rnaseq"]) != 0 or len(SAMPLES[wildcards.sample]["dnaseq"]) != 0:

        # make sure indels are called from rnaseq/dnaseq data

        counts += expand(
            "results/{sample}/quantification/allcounts.txt",
            sample=wildcards.sample,
        )
    return counts
