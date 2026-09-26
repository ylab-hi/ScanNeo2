"""Tests for epitope-window deduplication in BindingAffinities.start (PR #211).

A variant is annotated against every overlapping transcript, and isoforms
usually share the residues around it, so several rows of the variant-effects
table cut the same window. start() writes each distinct window to the
per-length fasta once and lets the other rows share its sequence number.

netMHCpan is replaced by a stub that reads the fasta start() wrote and reports
every k-mer of every sequence as a binder. Output rows are therefore fully
determined by which sequence each row's number points at, so a row wired to
the wrong window loses its epitopes or gains foreign ones.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))
import prediction  # noqa: E402

EPILEN = 9
VAR_POS = 10  # 0-based position of the variant residue in the subsequence

WT = "ACDEFGHIKLMNPQRSTVWYC"
MT_A = WT[:VAR_POS] + "W" + WT[VAR_POS + 1:]  # M11W
MT_B = WT[:VAR_POS] + "Y" + WT[VAR_POS + 1:]  # M11Y

# variant-effects columns start() reads, by index: 5 transcript_id,
# 9 wt_subseq, 10 mt_subseq, 12/13 aa_var_start/end, 14 vaf; 0-21 all present
HEADER = "\t".join(f"col{i}" for i in range(22))


def effects_row(transcript, wt, mt):
    cols = ["chr1", "100", "101", "G1", "GENE1", transcript, "DNA", "g1",
            "SNV", wt, mt, "100", str(VAR_POS), str(VAR_POS), "0.5", "10",
            "20", "1.0", "NA", "NA", "NA", "NA"]
    return "\t".join(cols)


def expected_epitopes(mt):
    """Every EPILEN-mer of mt that spans the variant position."""
    return {mt[p:p + EPILEN]
            for p in range(len(mt) - EPILEN + 1)
            if p <= VAR_POS <= p + EPILEN - 1}


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
    """Run start() on the given rows; return (fastas, epitopes per transcript)."""
    fastas = {}

    def every_kmer_binds(alleles, fnames, epilens, mhc_class, threads):
        res = {}
        for grp in ("wt", "mt"):
            res[grp] = {}
            for L in epilens:
                seqs = read_fasta(fnames[grp][L])
                fastas[(grp, L)] = seqs
                res[grp][L] = {
                    num: {seq[i:i + L]: ("HLA-A*02:01", i, i + L - 1, 50.0, 0.1)
                          for i in range(len(seq) - L + 1)}
                    for num, seq in seqs.items()
                }
        return res["wt"], res["mt"]

    monkeypatch.setattr(prediction.BindingAffinities,
                        "collect_binding_affinities",
                        staticmethod(every_kmer_binds))

    alleles = tmp_path / "mhc-I.tsv"
    alleles.write_text("HLA-A*02:01\tA*02:01\n")

    def run(rows):
        (tmp_path / "exitrons_variant_effects.tsv").write_text(
            "\n".join([HEADER, *rows]) + "\n")
        prediction.BindingAffinities(1).start(
            str(alleles), str(EPILEN), str(tmp_path), "mhc-I", "exitrons")

        epitopes = {}
        out = tmp_path / "exitrons_mhc-I_neoepitopes.txt"
        for line in out.read_text().splitlines()[1:]:
            cols = line.split("\t")
            epitopes.setdefault(cols[6], set()).add(cols[14])
        return fastas, epitopes

    return run


def test_identical_windows_are_written_once(run_start):
    fastas, _ = run_start([effects_row("T1", WT, MT_A),
                           effects_row("T2", WT, MT_A),
                           effects_row("T3", WT, MT_A)])

    assert len(fastas[("mt", EPILEN)]) == 1
    assert len(fastas[("wt", EPILEN)]) == 1


def test_rows_sharing_a_window_all_receive_its_epitopes(run_start):
    _, epitopes = run_start([effects_row("T1", WT, MT_A),
                             effects_row("T2", WT, MT_A)])

    assert epitopes["T1"] == expected_epitopes(MT_A)
    assert epitopes["T2"] == expected_epitopes(MT_A)


def test_each_row_points_at_its_own_window(run_start):
    """Shared and distinct windows interleaved: numbers must not drift.

    T3 repeats T1 after a distinct window has been written in between, so a
    reused number that were off by one would hand T3 the MT_B window.
    """
    fastas, epitopes = run_start([effects_row("T1", WT, MT_A),
                                  effects_row("T2", WT, MT_B),
                                  effects_row("T3", WT, MT_A)])

    mt = fastas[("mt", EPILEN)]
    assert sorted(mt) == [1, 2]
    assert len(set(mt.values())) == 2
    # WT is identical in all three rows
    assert len(fastas[("wt", EPILEN)]) == 1

    assert epitopes["T1"] == expected_epitopes(MT_A)
    assert epitopes["T2"] == expected_epitopes(MT_B)
    assert epitopes["T3"] == expected_epitopes(MT_A)
