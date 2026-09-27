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


def test_snv_is_the_single_changed_residue():
    assert bounds("ABCDEFGH", "ABXDEFGH") == (2, 2)


def test_mnv_and_delins_cover_the_changed_residues():
    assert bounds("ABCDEFGH", "ABXYEFGH") == (2, 3)
    assert bounds("ABCDEFGH", "ABXYZEFGH") == (2, 4)   # CD -> XYZ


def test_inframe_insertion_is_only_the_inserted_residues():
    # position-by-position comparison would run to the end (index 9)
    assert bounds("ABCDEFGH", "ABCXYDEFGH") == (3, 4)


def test_inframe_deletion_is_an_empty_region_at_the_junction():
    # DE deleted: C (2) and F (3) are now adjacent; end == start - 1
    assert bounds("ABCDEFGH", "ABCFGH") == (3, 2)


def test_deletion_inside_a_repeat_does_not_let_prefix_and_suffix_overlap():
    # XAAAAY -> XAAAY: the shared A's could belong to either side
    start, end = bounds("XAAAAY", "XAAAY")
    assert (start, end) == (4, 3)


def test_insertion_inside_a_repeat():
    assert bounds("XAAY", "XAAAY") == (3, 3)


def test_c_terminal_extension():
    assert bounds("ABCD", "ABCDXY") == (4, 5)


def test_c_terminal_deletion_has_no_new_peptide():
    # nothing follows the junction, so no peptide can span it
    assert bounds("ABCDEF", "ABCD") == (-1, -1)


def test_identical_sequences():
    assert bounds("ABCDEF", "ABCDEF") == (-1, -1)


def test_frameshift_runs_to_the_new_stop():
    assert bounds("ABCDEFGHIK", "ABCWXYZ", var_start=3) == (3, 6)


def test_scan_starts_at_var_start():
    # differences before var_start are not the variant's
    assert bounds("ABCDEFGH", "QBCDEFXH", var_start=3) == (6, 6)


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
