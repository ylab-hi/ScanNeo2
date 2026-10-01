"""Tests for filtering.SequenceSimilarity.self_similarity (#89).

The expected values are rows of a real neoepitope table (TESLA lung_patient12),
and the score is computed by the unmodified self_similarity / corr_kernel with
the workflow's BLOSUM62 matrix. When the mt epitope has no full-length wildtype
counterpart (it lies in sequence a frameshift or insertion created), the wt
epitope is empty or shorter and the score is ".".
"""

import sys
from pathlib import Path

import blosum as bl
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))

import filtering  # noqa: E402

# wt_epitope_seq, mt_epitope_seq, self-similarity as written by the workflow
ROWS = """wt_epitope_seq\tmt_epitope_seq\tself-similarity
VIGFAISQQK\tVVGFAISQQK\t0.651784218713048
MLTCPEAN\tSPVQRPTL\t2.084304252823848e-06
YH\tVPLRTVAV\t.
\tRTVAVLIRK\t.
"""


def scorer():
    # __init__ also runs the pathogen / proteome BLAST searches; only the
    # matrix it loads is needed for self_similarity
    s = object.__new__(filtering.SequenceSimilarity)
    s.matrix = bl.BLOSUM(str(REPO_ROOT / "workflow/scripts/prioritization/BLOSUM62-2.txt"))
    return s


def test_self_similarity_matches_workflow_output(tmp_path):
    table = tmp_path / "neoepitopes.txt"
    table.write_text(ROWS)
    df = pd.read_csv(table, sep="\t")  # as filtering.py reads it: empty wt -> NaN

    got = scorer().self_similarity(df["wt_epitope_seq"], df["mt_epitope_seq"])

    expected = [str(v) for v in df["self-similarity"]]
    assert [str(v) for v in got] == expected
