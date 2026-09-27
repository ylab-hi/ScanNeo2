"""collect_binding_affinities must not retain prediction results (#212).

Each unit's result is the full parsed output of one tool run. Once it is in
the prediction database it must be freed; a completed future kept alive
anywhere holds its result, and across a large transcriptomic input that
amounts to every prediction held in memory at once.

The stub returns a weak-referenceable dict per unit; the database connection
is wrapped so that, at every insert, the number of earlier results still alive
is recorded.
"""

import contextlib
import gc
import sys
import weakref
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/prioritization"))
import prediction  # noqa: E402


class Result(dict):
    """dict itself cannot be weakly referenced."""


def test_results_are_freed_once_inserted(tmp_path, monkeypatch):
    produced = []
    alive_at_insert = []

    def stub(call, fa_file, mhc_class):
        allele, L = call[3], int(call[4])
        result = Result({1: {"A" * L: (allele, 0, L - 1, 100.0, 1.0)}})
        produced.append(weakref.ref(result))
        return result

    class Counting:
        def __init__(self, db):
            self.db = db

        def __getattr__(self, name):
            return getattr(self.db, name)

        def executemany(self, sql, rows):
            if "INTO predictions" in sql:
                gc.collect()
                # results not yet inserted (this one included) are rightly
                # alive; anything beyond them was inserted and is being held
                inserted = len(alive_at_insert)
                alive = sum(r() is not None for r in produced)
                alive_at_insert.append(alive - (len(produced) - inserted))
            return self.db.executemany(sql, rows)

    orig_db = prediction.BindingAffinities.prediction_db

    @contextlib.contextmanager
    def counting_db(path):
        with orig_db(path) as db:
            yield Counting(db)

    monkeypatch.setattr(prediction.BindingAffinities, "_run_prediction",
                        staticmethod(stub))
    monkeypatch.setattr(prediction.BindingAffinities, "prediction_db",
                        staticmethod(counting_db))

    wt = "ACDEFGHIKLMNPQRSTVWYC"
    mt = wt[:10] + "W" + wt[11:]
    cols = ["chr1", "100", "101", "G1", "GENE1", "T1", "DNA", "g1", "SNV",
            wt, mt, "100", "10", "10", "0.5", "10", "20", "1.0",
            "NA", "NA", "NA", "NA"]
    header = "\t".join(f"col{i}" for i in range(22))
    (tmp_path / "exitrons_variant_effects.tsv").write_text(
        f"{header}\n" + "\t".join(cols) + "\n")
    alleles = tmp_path / "mhc-I.tsv"
    alleles.write_text("".join(f"HLA-A*0{i}:01\tHLA-A*0{i}:01\n" for i in range(1, 7)))

    prediction.BindingAffinities(1).start(
        str(alleles), "8,9,10,11", str(tmp_path), "mhc-I", "exitrons")

    # 6 alleles x 4 lengths x wt/mt
    assert len(alive_at_insert) == 48
    assert max(alive_at_insert) == 0, (
        f"earlier prediction results still alive at insert: {alive_at_insert}")
