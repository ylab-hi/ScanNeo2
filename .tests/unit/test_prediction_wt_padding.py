"""Tests for the wildtype padding in BindingAffinities.start (#218).

Where the mutant epitope runs past the end of the wildtype protein, the
wildtype presents nothing at those positions. effects.py marks them with '$'
(adjust_wildtype) and the reported wt epitope keeps the marker, while the
padding must never reach the prediction tool.

These drive the producing code: a variant-effects row with a '$'-padded
wildtype goes through start(), and the assertions read the neoepitope table
and the sequence database it wrote.
"""

import sqlite3
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))
import prediction  # noqa: E402

EPILEN = 9
WT = "ACDEFGHIKLMNPQRSTVWYC"  # 21 residues; the mutant insertion extends past it


def insertion(insert_len):
    """(wt padded as effects.py writes it, mt, aa_var_start, aa_var_end)."""
    mt = WT[:10] + "W" * insert_len + WT[10:]
    return WT.ljust(len(mt), "$"), mt, 10, 10 + insert_len - 1


def effects_row(transcript, wt, mt, var_start, var_end):
    return "\t".join(["chr1", "100", "101", "G1", "GENE1", transcript, "DNA", "g1",
                      "inframe_INS", wt, mt, "100", str(var_start), str(var_end),
                      "0.5", "10", "20", "1.0", "NA", "NA", "NA", "NA"])


def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = int(line[1:])
        else:
            seqs[name] = line
    return seqs


@pytest.fixture
def run_start(tmp_path, monkeypatch):
    """Run start() on one variant-effects row; return (written windows, rows)."""
    def every_kmer_binds(call, fa_file, mhc_class):
        allele, L = call[3], int(call[4])
        return {num: {seq[i:i + L]: (allele, i, i + L - 1, 50.0, 0.1)
                      for i in range(len(seq) - L + 1)}
                for num, seq in read_fasta(fa_file).items()}

    monkeypatch.setattr(prediction.BindingAffinities, "_run_prediction",
                        staticmethod(every_kmer_binds))
    alleles = tmp_path / "mhc-I.tsv"
    alleles.write_text("HLA-A*02:01\tA*02:01\n")

    def run(insert_len):
        wt, mt, a, b = insertion(insert_len)
        (tmp_path / "exitrons_variant_effects.tsv").write_text(
            "\t".join(f"col{i}" for i in range(22)) + "\n"
            + effects_row("T1", wt, mt, a, b) + "\n")
        prediction.BindingAffinities(1).start(
            str(alleles), str(EPILEN), str(tmp_path), "mhc-I", "exitrons")

        with sqlite3.connect(tmp_path / "exitrons_mhc-I_predictions.sqlite") as db:
            windows = [s for (s,) in db.execute("SELECT sequence FROM sequences")]
        out = (tmp_path / "exitrons_mhc-I_neoepitopes.txt").read_text().splitlines()
        cols = out[0].split("\t")
        rows = [dict(zip(cols, line.split("\t"))) for line in out[1:]]
        return windows, rows

    return run


def test_padding_marks_positions_past_the_end_of_the_wildtype(run_start):
    _, rows = run_start(4)

    # the epitope starting at the last position of the insertion reaches one
    # residue past the wildtype's end, so its last position has no wildtype
    padded = [r for r in rows if "$" in r["wt_epitope_seq"]]
    assert padded, [r["wt_epitope_seq"] for r in rows]
    for r in padded:
        assert len(r["wt_epitope_seq"]) == EPILEN          # padded, not truncated
        assert r["wt_epitope_seq"].rstrip("$") == r["wt_epitope_seq"].replace("$", "")
        assert r["wt_epitope_ic50"] == "."                 # nothing to predict
        assert r["agretopicity"] == "."
    assert ("WMNPQRSTV", "QRSTVWYC$") in [
        (r["mt_epitope_seq"], r["wt_epitope_seq"]) for r in padded]


def test_padding_never_reaches_the_prediction_input(run_start):
    windows, _ = run_start(4)

    assert windows
    assert all("$" not in w for w in windows)


def test_rows_with_a_full_length_wildtype_keep_their_affinity(run_start):
    _, rows = run_start(4)

    full = [r for r in rows if "$" not in r["wt_epitope_seq"]]
    assert full
    for r in full:
        assert len(r["wt_epitope_seq"]) == EPILEN
        assert r["wt_epitope_ic50"] not in (".", "")


def test_all_padding_when_the_wildtype_has_no_residue_there(run_start):
    # a long insertion puts a whole epitope past the end of the wildtype
    _, rows = run_start(12)

    assert "$" * EPILEN in [r["wt_epitope_seq"] for r in rows]
