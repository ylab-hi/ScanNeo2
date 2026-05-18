"""Tests for workflow/scripts/genotyping/combine_all_alleles.py.

The script is invoked by Snakemake as:
    python combine_all_alleles.py '<space-joined-tsvs>' <class> <output>

It reads each input TSV (lines of `source\tHLA-allele`), filters against
the class refset (mhc-I_refset.txt / mhc-II_refset.txt), merges identical
alleles from multiple sources, and writes a `source[,source]...\tallele`
output. On zero matched alleles it exits 1 with a per-source diagnostic
breakdown (see PR #83).
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow/scripts/genotyping/combine_all_alleles.py"


def run_script(input_arg, mhc_class, output_path):
    """Invoke the script with the same argv shape Snakemake uses.

    cwd must be the repo root because the script opens
    `workflow/scripts/genotyping/<class>_refset.txt` relative to cwd.
    """
    return subprocess.run(
        [sys.executable, str(SCRIPT), input_arg, mhc_class, str(output_path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def test_happy_path_merges_sources(tmp_path):
    in_tsv = tmp_path / "in.tsv"
    in_tsv.write_text("dna_tumor\tHLA-A*01:01\nrna_tumor\tHLA-A*01:01\n")
    out_tsv = tmp_path / "out.tsv"

    result = run_script(str(in_tsv), "mhc-I", out_tsv)

    assert result.returncode == 0, result.stderr + result.stdout
    assert out_tsv.read_text() == "dna_tumor,rna_tumor\tHLA-A*01:01\n"


def test_allele_chopped_to_two_fields(tmp_path):
    """HLA-A*02:01:01:01 → HLA-A*02:01 (only first two fields kept)."""
    in_tsv = tmp_path / "in.tsv"
    in_tsv.write_text("dna_tumor\tHLA-A*01:01:01:01\n")
    out_tsv = tmp_path / "out.tsv"

    result = run_script(str(in_tsv), "mhc-I", out_tsv)

    assert result.returncode == 0, result.stderr + result.stdout
    assert out_tsv.read_text() == "dna_tumor\tHLA-A*01:01\n"


def test_off_refset_allele_rejected_with_diagnostic(tmp_path):
    """Input has rows but no MHC values match the refset → exit 1."""
    in_tsv = tmp_path / "in.tsv"
    in_tsv.write_text("dna_tumor\tHLA-FAKE*99:99\n")
    out_tsv = tmp_path / "out.tsv"

    result = run_script(str(in_tsv), "mhc-I", out_tsv)

    assert result.returncode == 1
    assert "ERROR: No valid alleles" in result.stdout
    assert "0 matched the mhc-I reference set" in result.stdout


def test_empty_input_file_diagnosed(tmp_path):
    """0-byte input → exit 1, breakdown identifies the source as empty."""
    in_tsv = tmp_path / "in.tsv"
    in_tsv.write_text("")
    out_tsv = tmp_path / "out.tsv"

    result = run_script(str(in_tsv), "mhc-I", out_tsv)

    assert result.returncode == 1
    assert "EMPTY (0 bytes)" in result.stdout
    assert "upstream rule produced no output" in result.stdout


def test_malformed_line_rejected(tmp_path):
    """Row without exactly 2 tab-separated columns → exit 1."""
    in_tsv = tmp_path / "in.tsv"
    in_tsv.write_text("only-one-column-no-tabs\n")
    out_tsv = tmp_path / "out.tsv"

    result = run_script(str(in_tsv), "mhc-I", out_tsv)

    assert result.returncode == 1
    assert "Invalid input file" in result.stdout


def test_mixed_empty_and_populated_sources(tmp_path):
    """Two input files: one empty, one with off-refset content. Diagnostic
    should distinguish the two failure modes per source."""
    empty_tsv = tmp_path / "empty.tsv"
    empty_tsv.write_text("")
    bad_tsv = tmp_path / "bad.tsv"
    bad_tsv.write_text("rna_tumor\tHLA-FAKE*99:99\n")
    out_tsv = tmp_path / "out.tsv"

    # Snakemake passes multiple inputs as a single space-joined string.
    input_arg = f"{empty_tsv} {bad_tsv}"
    result = run_script(input_arg, "mhc-I", out_tsv)

    assert result.returncode == 1
    assert "EMPTY (0 bytes)" in result.stdout
    assert "0 matched the mhc-I reference set" in result.stdout
