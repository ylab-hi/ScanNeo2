"""Tests for VariantEffects.determine_var_bnds (#216).

The altered region used to be found by comparing wildtype and mutant position
by position. Past an in-frame insertion or deletion the mutant is the wildtype
shifted, so every later position compared unequal and the region ran to the
protein end: wildtype peptides were predicted, and reported whenever they
overlapped that region. The region is now bounded by the common prefix and the
common suffix, and an in-frame deletion becomes an empty region at the junction.
"""

import importlib
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))

# reference.py needs pyfaidx, which the CI test env does not install;
# determine_var_bnds does not touch it
sys.modules.setdefault("reference", types.ModuleType("reference"))
# another test module stubs `effects` itself; make sure this one gets the real one
if not hasattr(sys.modules.get("effects"), "VariantEffects"):
    sys.modules.pop("effects", None)
effects = importlib.import_module("effects")

import prediction  # noqa: E402


def bounds(wt, mt, var_start=0):
    ve = effects.VariantEffects.__new__(effects.VariantEffects)
    ve.data = {"var_start": var_start, "mt_seq": mt,
               "wt_seq": effects.VariantEffects.adjust_wildtype(wt, mt)}
    return ve.determine_var_bnds(wt)


@pytest.mark.parametrize("wt, mt, var_start, expected", [
    pytest.param("ABCDEFGH", "ABXDEFGH", 0, (2, 2), id="snv"),
    pytest.param("ABCDEFGH", "ABXYEFGH", 0, (2, 3), id="mnv"),
    pytest.param("ABCDEFGH", "ABXYZEFGH", 0, (2, 4), id="delins"),
    # position-by-position comparison would run to the end (index 9)
    pytest.param("ABCDEFGH", "ABCXYDEFGH", 0, (3, 4), id="inframe-insertion"),
    # DE deleted: C (2) and F (3) are now adjacent; end == start - 1
    pytest.param("ABCDEFGH", "ABCFGH", 0, (3, 2), id="inframe-deletion"),
    # the shared A's could belong to either side; prefix and suffix must not overlap
    pytest.param("XAAAAY", "XAAAY", 0, (4, 3), id="deletion-in-repeat"),
    pytest.param("XAAY", "XAAAY", 0, (3, 3), id="insertion-in-repeat"),
    pytest.param("ABCD", "ABCDXY", 0, (4, 5), id="c-terminal-extension"),
    # nothing follows the junction, so no peptide can span it
    pytest.param("ABCDEF", "ABCD", 0, (-1, -1), id="c-terminal-deletion"),
    # the mutant is a piece of the wildtype: an empty region before position 0,
    # which the overlap filter (p <= -1) can never select
    pytest.param("ABCDEFGH", "DEFGH", 0, (0, -1), id="n-terminal-deletion"),
    # the new residue L happens to equal the wildtype's last residue; the
    # region is empty, so only peptides spanning C|L are selected
    pytest.param("ABCDEFGHL", "ABCL", 3, (3, 2), id="stop-gain-matching-c-terminus"),
    pytest.param("ABCDEF", "ABCDEF", 0, (-1, -1), id="identical"),
    pytest.param("ABCDEFGHIK", "ABCWXYZ", 3, (3, 6), id="frameshift"),
    # differences before var_start are not the variant's
    pytest.param("ABCDEFGH", "QBCDEFXH", 3, (6, 6), id="scan-starts-at-var-start"),
    # fusion partner without a known wildtype: start lies beyond the wildtype,
    # so nothing can be shared and the region runs to the end
    pytest.param("", "ABCD", 1, (1, 3), id="fusion-without-wildtype"),
])
def test_bounds(wt, mt, var_start, expected):
    assert bounds(wt, mt, var_start) == expected


def test_change_entry_bounds_on_the_unpadded_wildtype(monkeypatch):
    """adjust_wildtype truncates the stored wildtype to the mutant's length,
    which removes the tail a deletion's suffix is matched on; change_entry must
    hand determine_var_bnds the wildtype as passed in."""
    monkeypatch.setattr(effects.VariantEffects, "determine_NMD", lambda self, nmd: None)
    monkeypatch.setattr(effects.VariantEffects, "get_counts", lambda self: None)
    ve = effects.VariantEffects.__new__(effects.VariantEffects)

    ve.change_entry(chrom="chr1", start=100, end=101, gene_id="G1",
                    gene_name="GENE1", transcript_id="T1", transcript=None,
                    transcript_bp=None, source="DNA", group="g1",
                    var_type="inframe_DEL", var_start=0,
                    wt_seq="ABCDEFGH", mt_seq="ABCFGH",
                    vaf=0.5, ao=10, dp=20, nmd_event=None)

    # on the truncated wildtype "ABCDEF" the region would be (3, 5)
    assert (ve.data["aa_var_start"], ve.data["aa_var_end"]) == (3, 2)


def test_deletion_yields_only_junction_spanning_epitopes(tmp_path, monkeypatch):
    """effects' bounds and subsequence fed through prediction's window and filter.

    Every k-mer is reported as a binder, so the table lists exactly the
    epitopes the region selects. Only those containing both residues flanking
    the junction are new; any other k-mer of the mutant also occurs in the
    wildtype.
    """
    wt = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWYAVRDLNSLCIDTR"
    cut = 60
    mt = wt[:cut] + wt[cut + 12:]              # 12-residue in-frame deletion

    ve = effects.VariantEffects.__new__(effects.VariantEffects)
    ve.data = {"var_start": 0, "mt_seq": mt,
               "wt_seq": effects.VariantEffects.adjust_wildtype(wt, mt)}
    ve.data["aa_var_start"], ve.data["aa_var_end"] = ve.determine_var_bnds(wt)
    ve.determine_subsequence()
    start, end = ve.data["aa_var_start"], ve.data["aa_var_end"]
    assert end == start - 1

    L = 9

    def every_kmer_binds(call, fa_file, *rest):
        allele = call[3]
        res, num = {}, None
        for line in Path(fa_file).read_text().splitlines():
            if line.startswith(">"):
                num = int(line[1:])
                continue
            res[num] = {line[i:i + L]: (allele, i, i + L - 1, 50.0, 0.1)
                        for i in range(len(line) - L + 1)}
        return res

    monkeypatch.setattr(prediction.BindingAffinities, "_run_prediction",
                        staticmethod(every_kmer_binds))

    cols = ["chr1", "100", "101", "G1", "GENE1", "T1", "DNA", "g1",
            "inframe_DEL", ve.data["wt_subseq"], ve.data["mt_subseq"], "100",
            str(start), str(end), "0.5", "10", "20", "1.0",
            "NA", "NA", "NA", "NA"]
    header = "\t".join(f"col{i}" for i in range(22))
    (tmp_path / "exitrons_variant_effects.tsv").write_text(
        header + "\n" + "\t".join(cols) + "\n")
    (tmp_path / "mhc-I.tsv").write_text("HLA-A*02:01\tHLA-A*02:01\n")

    prediction.BindingAffinities(1).start(
        str(tmp_path / "mhc-I.tsv"), str(L), str(tmp_path), "mhc-I", "exitrons")

    lines = (tmp_path / "exitrons_mhc-I_neoepitopes.txt").read_text().splitlines()
    reported = {line.split("\t")[14] for line in lines[1:]}

    junction = {mt[p:p + L] for p in range(cut - L + 1, cut)}
    assert reported == junction
    assert not any(k in wt for k in reported)
