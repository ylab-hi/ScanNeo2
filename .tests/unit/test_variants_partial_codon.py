"""Tests for variants.Variants.strip_partial_codon (#214).

VEP translates the incomplete first codon of a cds_start_NF transcript to a
leading X. Left in the wildtype, it fails a whole prediction batch; for a
frameshift past it, mt_seq = wt[:var_start] + DownstreamProtein starts with
the X and was truncated to nothing, silently dropping the row.
"""

import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

# strip_partial_codon touches none of variants.py's heavy imports
for _dep in ("vcfpy", "reference", "effects"):
    sys.modules.setdefault(_dep, types.ModuleType(_dep))
sys.path.insert(0, str(PRIORITIZATION))

import variants  # noqa: E402

strip = variants.Variants.strip_partial_codon


@pytest.mark.parametrize(
    "wt, var_start, csq, expected",
    [
        pytest.param("MKTAY", 2, "frameshift", ("MKTAY", 2), id="no-x-untouched"),
        pytest.param("XKTAY", 3, "frameshift", ("KTAY", 2), id="frameshift-after-x"),
        pytest.param("XKTAY", 3, "missense", ("KTAY", 2), id="missense-after-x"),
        pytest.param("XKTAY", 0, "frameshift", ("KTAY", 0), id="frameshift-on-x"),
        pytest.param("XKTAY", 0, "missense", (None, None), id="missense-on-x"),
        pytest.param("XKTAY", 0, "inframe_DEL", (None, None), id="inframe-on-x"),
    ],
)
def test_strip_partial_codon(wt, var_start, csq, expected):
    assert strip(wt, var_start, csq) == expected


def test_frameshift_past_x_keeps_its_mutant():
    # the construction variants.py uses for a frameshift mutant
    wt, var_start = strip("XKTAYNE", 3, "frameshift")
    mt = wt[:var_start] + "PQRS"
    assert mt == "KTPQRS"
    assert "X" not in wt and "X" not in mt
    # still points at the residue VEP's Protein_position named (A, at index 3 with the X)
    assert wt[var_start] == "A"
