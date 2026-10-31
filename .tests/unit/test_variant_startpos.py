"""Tests for variants.Variants.get_variant_startpos (#262).

VEP writes a pure insertion as Amino_acids '-/X' at Protein_position 'N-N+1',
between residues N and N+1, so an inframe insertion continues from index N.
A frameshift keeps N-1: its DownstreamProtein is not anchored on
Protein_position (#269). The fields below are VEP's own, from the TESLA
annotation.
"""

import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

# get_variant_startpos touches none of variants.py's heavy imports
for _dep in ("vcfpy", "reference", "effects"):
    sys.modules.setdefault(_dep, types.ModuleType(_dep))
sys.path.insert(0, str(PRIORITIZATION))

import variants  # noqa: E402

startpos = variants.Variants.get_variant_startpos

# HES4 ENST00000304952, residues 1-75
HES4 = "MAADTPGKPSASPMAGAPASASRTPDKPRSAAEHRKSSKPVMEKRRRARINESLAQLKTLILDALRKESSRHSK"


def test_inframe_insertion_after_residue_n():
    # intron_retention at chr1:999613; the retained intron opens after E68
    start = startpos("68-69", "-/VSRGRRPETRSLGRS*PGPPAYPRPQ", "inframe_INS")
    assert HES4[:start].endswith("LRKE")


@pytest.mark.parametrize(
    "protein_position, aa_change, csq, expected",
    [
        pytest.param("161-162", "-/X", "frameshift", 160, id="frameshift-insertion"),
        pytest.param("143", "S/FX", "frameshift", 142, id="frameshift-substitution"),
        pytest.param("212-248", "GKGAWP/X", "frameshift", 211, id="frameshift-deletion"),
        pytest.param("67", "K/KE", "inframe_INS", 66, id="insertion-written-with-anchor"),
        pytest.param("40-42", "KES/-", "inframe_DEL", 39, id="inframe-deletion"),
        pytest.param("10", "A/V", "missense", 9, id="missense"),
    ],
)
def test_startpos(protein_position, aa_change, csq, expected):
    assert startpos(protein_position, aa_change, csq) == expected
