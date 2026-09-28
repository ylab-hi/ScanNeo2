"""Tests for the non-standard-residue guard in BindingAffinities.start (#214).

The IEDB tools reject a whole batch file if any sequence in it has a residue
outside the 20 standard amino acids, so one bad window used to take its
batch-mates down with it. start() now keeps such a window out of the fasta:
that row loses only the affected window's predictions.

netMHCpan is replaced by a stub that, like the real tool, fails the whole
batch when any sequence carries a non-standard residue, and otherwise reports
every k-mer as a binder.
"""

import sqlite3
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))
import prediction  # noqa: E402

EPILEN = 9
VAR_POS = 10

WT = "ACDEFGHIKLMNPQRSTVWYC"
MT = WT[:VAR_POS] + "W" + WT[VAR_POS + 1:]
WT_X = WT[:VAR_POS - 2] + "X" + WT[VAR_POS - 1:]  # X inside the wt window
MT_OTHER = WT[:VAR_POS] + "Y" + WT[VAR_POS + 1:]

IEDB_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")
HEADER = "\t".join(f"col{i}" for i in range(22))


def effects_row(transcript, wt, mt):
    cols = ["chr1", "100", "101", "G1", "GENE1", transcript, "DNA", "g1",
            "SNV", wt, mt, "100", str(VAR_POS), str(VAR_POS), "0.5", "10",
            "20", "1.0", "NA", "NA", "NA", "NA"]
    return "\t".join(cols)


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
    def iedb_like(call, fa_file, mhc_class):
        allele, L = call[3], int(call[4])
        seqs = read_fasta(fa_file)
        if any(not IEDB_RESIDUES.issuperset(s) for s in seqs.values()):
            return None  # the tool rejects the whole file
        return {num: {seq[i:i + L]: (allele, i, i + L - 1, 50.0, 0.1)
                      for i in range(len(seq) - L + 1)}
                for num, seq in seqs.items()}

    monkeypatch.setattr(prediction.BindingAffinities, "_run_prediction",
                        staticmethod(iedb_like))

    alleles = tmp_path / "mhc-I.tsv"
    alleles.write_text("HLA-A*02:01\tA*02:01\n")

    def run(rows):
        (tmp_path / "exitrons_variant_effects.tsv").write_text(
            "\n".join([HEADER, *rows]) + "\n")
        prediction.BindingAffinities(1).start(
            str(alleles), str(EPILEN), str(tmp_path), "mhc-I", "exitrons")

        with sqlite3.connect(tmp_path / "exitrons_mhc-I_predictions.sqlite") as db:
            windows = [seq for (seq,) in db.execute("SELECT sequence FROM sequences")]
        rows_out = {}
        out = tmp_path / "exitrons_mhc-I_neoepitopes.txt"
        header, *lines = out.read_text().splitlines()
        cols = header.split("\t")
        for line in lines:
            r = dict(zip(cols, line.split("\t")))
            rows_out.setdefault(r["transcript_id"], []).append(r)
        return windows, rows_out

    return run


def test_nonstandard_window_does_not_fail_its_batch_mates(run_start):
    windows, rows = run_start([effects_row("T1", WT_X, MT),
                               effects_row("T2", WT, MT_OTHER)])

    assert all(IEDB_RESIDUES.issuperset(w) for w in windows)
    # T2 shares both batches with T1 and is fully predicted
    assert rows["T2"]
    assert all(r["wt_epitope_ic50"] != "." for r in rows["T2"])
    # T1 loses only its withheld wt window: neoepitopes, but no wt affinity
    assert rows["T1"]
    assert {r["wt_epitope_ic50"] for r in rows["T1"]} == {"."}
