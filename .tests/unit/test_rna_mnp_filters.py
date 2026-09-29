"""Tests for the RNA filters' MNP handling (#221).

Mutect2 merges adjacent substitutions on one haplotype into an MNP. In RNA
that turns a cluster of A-to-I edits into one record (AA>GG), and a phased
pair of germline SNVs into one record (GG>AT). Both filters used to compare
whole records against single-base references and so never matched them.
"""

import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts"))
# only the pure decision functions are tested; pysam is needed by main() alone
sys.modules.setdefault("pysam", types.ModuleType("pysam"))

import annotate_rna_editing  # noqa: E402
import subtract_germline_mnps  # noqa: E402

POS = 100


def redi(*sites):
    """site_at() over a REDIportal holding exactly the given (pos, ref, ed)."""
    return lambda p, r, e: (p, r, e) if (p, r, e) in sites else None


@pytest.mark.parametrize(
    "ref, alt, sites, flagged",
    [
        pytest.param("A", "G", [(POS, "A", "G")], True, id="snv-known-edit"),
        pytest.param("A", "G", [], False, id="snv-unknown-site"),
        pytest.param("C", "T", [(POS, "C", "T")], False, id="snv-not-a-to-i"),
        pytest.param("AA", "GG", [(POS, "A", "G"), (POS + 1, "A", "G")], True,
                     id="mnp-all-bases-known-edits"),
        pytest.param("TT", "CC", [(POS, "T", "C"), (POS + 1, "T", "C")], True,
                     id="mnp-minus-strand-edits"),
        pytest.param("AGA", "GGG", [(POS, "A", "G"), (POS + 2, "A", "G")], True,
                     id="mnp-unchanged-middle-base"),
        pytest.param("AA", "GG", [(POS, "A", "G")], False, id="mnp-one-base-not-in-rediportal"),
        pytest.param("AC", "GT", [(POS, "A", "G"), (POS + 1, "C", "T")], False,
                     id="mnp-one-base-not-a-to-i"),
        pytest.param("A", "AG", [(POS, "A", "G")], False, id="indel-never-editing"),
    ],
)
def test_editing_sites(ref, alt, sites, flagged):
    rows = annotate_rna_editing.editing_sites(POS, ref, alt, redi(*sites))
    assert (rows is not None) == flagged


def germline(*snvs):
    return lambda p, r, a: (p, r, a) in snvs


@pytest.mark.parametrize(
    "ref, alts, snvs, dropped",
    [
        pytest.param("GG", ("AT",), [(POS, "G", "A"), (POS + 1, "G", "T")], True,
                     id="phased-germline-pair"),
        pytest.param("GG", ("AT",), [(POS, "G", "A")], False, id="one-base-somatic"),
        pytest.param("GG", ("AT",), [(POS, "G", "A"), (POS + 1, "G", "C")], False,
                     id="germline-other-alt"),
        pytest.param("GCG", ("ACT",), [(POS, "G", "A"), (POS + 2, "G", "T")], True,
                     id="unchanged-middle-base"),
        pytest.param("GG", ("AT", "CT"), [(POS, "G", "A"), (POS + 1, "G", "T")], False,
                     id="multiallelic-one-alt-somatic"),
        pytest.param("G", ("A",), [(POS, "G", "A")], False, id="snv-left-to-isec"),
    ],
)
def test_fully_germline(ref, alts, snvs, dropped):
    assert subtract_germline_mnps.fully_germline(POS, ref, alts, germline(*snvs)) == dropped
