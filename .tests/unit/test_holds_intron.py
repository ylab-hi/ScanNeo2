"""Tests for variants.Variants.holds_intron (#262).

An intron_retention record is only the retention spladder observed on a
transcript that splices exactly that intron. The GTF rows below are Ensembl's
own (attributes trimmed to the ones Annotation parses), loaded through the real
Annotation so that parsing, exon_map and the check run as in the workflow.
Records follow altsplc2vcf: POS is the base before the intron, END its last.
"""

import importlib.util
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

# the CI test env lacks pyfaidx; Annotation opens the fasta but holds_intron
# never reads it
try:
    import pyfaidx  # noqa: F401
except ImportError:
    sys.modules["pyfaidx"] = types.ModuleType("pyfaidx")
    sys.modules["pyfaidx"].Fasta = lambda path: None

# other tests stub `reference` in sys.modules, so load the real module by path
spec = importlib.util.spec_from_file_location(
    "reference_holds_intron", PRIORITIZATION / "reference.py")
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)

for _dep in ("vcfpy", "reference", "effects"):
    sys.modules.setdefault(_dep, types.ModuleType(_dep))
sys.path.insert(0, str(PRIORITIZATION))

import variants  # noqa: E402

GTF = """\
chr1	havana_tagene	transcript	827635	859408	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741";
chr1	havana_tagene	exon	827635	827775	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "1";
chr1	havana_tagene	exon	829003	829104	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "2";
chr1	havana_tagene	exon	847654	847806	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "3";
chr1	havana_tagene	exon	849484	849602	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "4";
chr1	havana_tagene	exon	851927	852110	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "5";
chr1	havana_tagene	exon	852671	852766	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "6";
chr1	havana_tagene	exon	853391	855121	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "7";
chr1	havana_tagene	exon	856449	859408	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000666741"; exon_number "8";
chr1	havana_tagene	transcript	827598	859446	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118";
chr1	havana_tagene	exon	827598	827775	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118"; exon_number "1";
chr1	havana_tagene	exon	829003	829104	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118"; exon_number "2";
chr1	havana_tagene	exon	851927	852110	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118"; exon_number "3";
chr1	havana_tagene	exon	852671	852766	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118"; exon_number "4";
chr1	havana_tagene	exon	853391	859446	.	+	.	gene_id "ENSG00000228794"; transcript_id "ENST00000445118"; exon_number "5";
chr1	ensembl_havana	transcript	998964	1000097	.	-	.	gene_id "ENSG00000188290"; transcript_id "ENST00000304952";
chr1	ensembl_havana	exon	999866	1000097	.	-	.	gene_id "ENSG00000188290"; transcript_id "ENST00000304952"; exon_number "1";
chr1	ensembl_havana	exon	999692	999787	.	-	.	gene_id "ENSG00000188290"; transcript_id "ENST00000304952"; exon_number "2";
chr1	ensembl_havana	exon	999526	999613	.	-	.	gene_id "ENSG00000188290"; transcript_id "ENST00000304952"; exon_number "3";
chr1	ensembl_havana	exon	998964	999432	.	-	.	gene_id "ENSG00000188290"; transcript_id "ENST00000304952"; exon_number "4";
"""


@pytest.fixture(scope="module")
def annotation(tmp_path_factory):
    d = tmp_path_factory.mktemp("annotation")
    (d / "genome.gtf").write_text(GTF)
    (d / "genome.fasta").write_text(">chr1\nN\n")
    return reference.Annotation(str(d / "genome.fasta"), str(d / "genome.gtf"))


@pytest.mark.parametrize(
    "transcript_id, pos, end, expected",
    [
        # + strand: exon 7 ends at 855121, exon 8 starts at 856449
        pytest.param("ENST00000666741", 855121, 856448, True, id="plus-splices-intron"),
        # same gene, one exon 853391-859446 already holds the intron
        pytest.param("ENST00000445118", 855121, 856448, False, id="plus-intron-is-exon"),
        # - strand HES4: exon_map walks 3'->5' genomically, adjacency is genomic
        pytest.param("ENST00000304952", 999613, 999691, True, id="minus-splices-intron"),
        # spans exon 2 (999692-999787): two introns, not one
        pytest.param("ENST00000304952", 999613, 999865, False, id="minus-exon-inside-span"),
        # flank matches on one side only
        pytest.param("ENST00000304952", 999613, 999700, False, id="minus-end-off"),
        pytest.param("ENST00000304952", 999612, 999691, False, id="minus-pos-off-by-one"),
        pytest.param("ENST00000000000", 999613, 999691, False, id="transcript-not-in-gtf"),
    ],
)
def test_holds_intron(annotation, transcript_id, pos, end, expected):
    assert variants.Variants.holds_intron(annotation, transcript_id, pos, end) is expected
