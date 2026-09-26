"""Tests for per-allele binding affinities in BindingAffinities (#212).

Prediction runs one unit per (allele, epitope length, wt/mt, batch). Several
alleles can bind the same epitope with different affinities, so the neoepitope
table carries one row per binding allele, and each row's wt affinity must come
from that same allele: agretopicity and the ranking score compare the two.

_run_prediction is replaced by a stub that answers from a fixed table, so the
real merge in collect_binding_affinities runs. With a single worker, units
complete in submission order, i.e. in allele-file order, which lets the tests
check that the output does not depend on it.
"""

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
MT_EPITOPE = MT[2:2 + EPILEN]  # spans the variant
WT_EPITOPE = WT[2:2 + EPILEN]  # same coordinates in the wildtype

A, B = "HLA-A*02:01", "HLA-B*07:02"

HEADER = "\t".join(f"col{i}" for i in range(22))
ROW = "\t".join(["chr1", "100", "101", "G1", "GENE1", "T1", "DNA", "g1",
                 "SNV", WT, MT, "100", str(VAR_POS), str(VAR_POS), "0.5",
                 "10", "20", "1.0", "NA", "NA", "NA", "NA"])


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
    """Run start() with the given {(group, allele): {epitope: ic50}} answers.

    Returns the output rows as dicts keyed by header name.
    """
    def run(answers, allele_order=(A, B)):
        def stub(call, fa_file, group, mhc_class):
            allele = call[3]
            table = answers.get((group, allele), {})
            result = {}
            for num, seq in read_fasta(fa_file).items():
                for epitope, ic50 in table.items():
                    start = seq.find(epitope)
                    if start != -1:
                        result.setdefault(num, {})[epitope] = (
                            allele, start, start + len(epitope) - 1, ic50, 1.0)
            return result

        monkeypatch.setattr(prediction.BindingAffinities, "_run_prediction",
                            staticmethod(stub))

        alleles = tmp_path / "mhc-I.tsv"
        alleles.write_text("".join(f"{a}\t{a}\n" for a in allele_order))
        (tmp_path / "exitrons_variant_effects.tsv").write_text(
            f"{HEADER}\n{ROW}\n")

        prediction.BindingAffinities(1).start(
            str(alleles), str(EPILEN), str(tmp_path), "mhc-I", "exitrons")

        lines = (tmp_path / "exitrons_mhc-I_neoepitopes.txt").read_text().splitlines()
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    return run


BOTH_BIND = {
    ("mt", A): {MT_EPITOPE: 100.0},
    ("mt", B): {MT_EPITOPE: 200.0},
    ("wt", A): {WT_EPITOPE: 1000.0},
    ("wt", B): {WT_EPITOPE: 4000.0},
}


def test_every_binding_allele_gets_its_own_row(run_start):
    rows = run_start(BOTH_BIND)

    assert [r["allele"] for r in rows] == [A, B]
    assert all(r["mt_epitope_seq"] == MT_EPITOPE for r in rows)
    assert all(r["wt_epitope_seq"] == WT_EPITOPE for r in rows)


def test_wt_affinity_comes_from_the_same_allele(run_start):
    rows = {r["allele"]: r for r in run_start(BOTH_BIND)}

    assert float(rows[A]["mt_epitope_ic50"]) == 100.0
    assert float(rows[A]["wt_epitope_ic50"]) == 1000.0
    assert float(rows[A]["agretopicity"]) == pytest.approx(0.1)

    assert float(rows[B]["mt_epitope_ic50"]) == 200.0
    assert float(rows[B]["wt_epitope_ic50"]) == 4000.0
    assert float(rows[B]["agretopicity"]) == pytest.approx(0.05)


def test_wt_without_a_prediction_for_that_allele_is_empty(run_start):
    """A wt result for allele A must not stand in for allele B."""
    answers = {k: v for k, v in BOTH_BIND.items() if k != ("wt", B)}
    rows = {r["allele"]: r for r in run_start(answers)}

    assert float(rows[A]["wt_epitope_ic50"]) == 1000.0
    assert float(rows[B]["mt_epitope_ic50"]) == 200.0
    # None is written as '.'
    assert rows[B]["wt_epitope_ic50"] == "."
    assert rows[B]["wt_epitope_rank"] == "."
    assert rows[B]["agretopicity"] == "."


def test_output_does_not_depend_on_completion_order(run_start):
    assert run_start(BOTH_BIND, allele_order=(A, B)) == \
        run_start(BOTH_BIND, allele_order=(B, A))
