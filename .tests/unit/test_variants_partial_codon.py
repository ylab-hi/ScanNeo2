"""Tests for variants.Variants.trim_at_unknown_residue (#214).

VEP translates the incomplete first codon of a cds_start_NF transcript to a
leading X. Left in the wildtype, it fails a whole prediction batch; for a
frameshift past it, mt_seq = wt[:var_start] + DownstreamProtein starts with
the X and was truncated to nothing, silently dropping the row. An X anywhere
else must still be trimmed (before the variant) or the entry skipped.
"""

import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

# trim_at_unknown_residue touches none of variants.py's heavy imports
for _dep in ("vcfpy", "reference", "effects"):
    sys.modules.setdefault(_dep, types.ModuleType(_dep))
sys.path.insert(0, str(PRIORITIZATION))

import variants  # noqa: E402

trim = variants.Variants.trim_at_unknown_residue


@pytest.mark.parametrize(
    "wt, var_start, csq, expected",
    [
        pytest.param("MKTAY", 2, "frameshift", ("MKTAY", 2), id="no-x-untouched"),
        # leading X (cds_start_NF partial codon); var_start keeps naming the same residue
        pytest.param("XKTAY", 3, "frameshift", ("KTAY", 2), id="frameshift-after-x"),
        pytest.param("XKTAY", 3, "missense", ("KTAY", 2), id="missense-after-x"),
        pytest.param("XKTAY", 3, "inframe_INS", ("KTAY", 2), id="inframe-after-x"),
        pytest.param("XKTAY", 0, "frameshift", ("KTAY", 0), id="frameshift-on-x"),
        pytest.param("XKTAY", 0, "missense", (None, None), id="missense-on-x"),
        pytest.param("XKTAY", 0, "inframe_DEL", (None, None), id="inframe-on-x"),
        # internal X
        pytest.param("MKXAYNE", 4, "missense", ("AYNE", 1), id="internal-x-before-variant"),
        pytest.param("MKTAYXE", 2, "missense", (None, None), id="internal-x-after-variant"),
        pytest.param("MKXAYNE", 2, "inframe_DEL", (None, None), id="internal-x-on-variant"),
        pytest.param("MKTAYXE", 2, "frameshift", (None, None), id="frameshift-internal-x-after"),
    ],
)
def test_trim_at_unknown_residue(wt, var_start, csq, expected):
    assert trim(wt, var_start, csq) == expected
