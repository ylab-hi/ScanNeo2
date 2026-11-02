"""The mutant peptide Variants writes for a pure inframe insertion.

The record is real: a 24-residue insertion from a TESLA tumor
(chr11:46542984), reduced to the CSQ entry of ENST00000458649, with that
transcript's GTF rows. VEP writes it as Amino_acids '-/SQPV...SGI' at
Protein_position 344-345, between V344 and Q345. The expected mt_subseq occurs
verbatim in the genome-truth mutant protein: REF->ALT applied to the
transcript from the GTF and genome, spliced and translated.

The inframe path reads no genome sequence, so a stub FASTA stands in for it.
"""

import importlib
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"
DATA = Path(__file__).resolve().parent / "data/inframe_insertion"

REAL = ("vcfpy", "reference", "effects", "utility", "variants")

MT_SUBSEQ = "SLLPHQDSVPPASARATTPSFSFVSQPVRRTPLLHNFLHMLSSRSSGIQTEPFHPPEQASSTQQDQGLLNRP"


@pytest.fixture
def real_variants():
    """The real variants module, whatever stubs other tests left in sys.modules."""
    saved = {name: sys.modules.pop(name) for name in REAL if name in sys.modules}
    sys.path.insert(0, str(PRIORITIZATION))
    try:
        yield importlib.import_module("variants")
    finally:
        sys.path.remove(str(PRIORITIZATION))
        for name in REAL:
            sys.modules.pop(name, None)
        sys.modules.update(saved)


def test_insertion_follows_residue_n(real_variants, tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr11\nN\n")
    options = types.SimpleNamespace(
        reference=str(fasta),
        anno=str(DATA / "annotation.gtf"),
        counts="",
        output_dir=str(tmp_path / "out"),
    )
    real_variants.Variants(str(DATA / "record.vcf"), options, "somatic.short.indels")

    lines = (tmp_path / "out/somatic.short.indels_variant_effects.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    assert len(rows) == 1
    row = rows[0]
    assert row["transcript_id"] == "ENST00000458649"
    assert row["var_type"] == "inframe_INS"
    # the insertion follows V344: ...FSFV|SQPV...SGI|QTEP...
    assert row["mt_subseq"] == MT_SUBSEQ
    assert row["var_start"] == "344"
