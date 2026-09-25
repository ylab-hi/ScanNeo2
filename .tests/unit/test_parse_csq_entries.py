"""Tests for variants.Variants.parse_csq_entries (PR #204).

VEP emits one CSQ entry per transcript, so a single variant record carries an
entry for every transcript it touches -- median 6, and 44.5% of records span
more than one gene. The parser returned from inside its loop, so it examined
only the *first* entry and discarded the rest: the gene a variant was
attributed to came down to VEP's ordering.

The concrete casualty was NRAS Q61H, whose record carries 47 entries -- 46 for
the overlapping CSDE1 and 1 for NRAS. CSDE1 came first, so the driver was
dropped and its validated neoepitope never reached the neoepitope table. These
tests pin the behaviour, because the failure mode is silent: no error, no
warning, just missing neoantigens.
"""

import sys
import types
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

# variants.py pulls in vcfpy and the sibling reference/effects modules, which
# drag pyfaidx/gffutils/pandas. parse_csq_entries touches none of them, so stub
# the imports rather than installing the whole runtime into the test env.
for _dep in ("vcfpy", "reference", "effects"):
    sys.modules.setdefault(_dep, types.ModuleType(_dep))
sys.path.insert(0, str(PRIORITIZATION))

import variants  # noqa: E402

CSQ_FORMAT = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|Amino_acids|Protein_position"


def entry(allele, consequence, symbol, feature, aa="", pos=""):
    return f"{allele}|{consequence}|MODERATE|{symbol}|ENSG|{feature}|{aa}|{pos}"


def parser():
    """A Variants instance without running __init__ (which reads a VCF)."""
    return variants.Variants.__new__(variants.Variants)


def test_returns_every_entry_for_the_allele():
    """All matching entries, not just the first."""
    csq = [
        entry("A", "downstream_gene_variant", "CSDE1", "ENST_A"),
        entry("A", "intron_variant", "CSDE1", "ENST_B"),
        entry("A", "missense_variant", "NRAS", "ENST_C", aa="Q/H", pos="61"),
    ]
    got = parser().parse_csq_entries(csq, CSQ_FORMAT, "A")

    assert len(got) == 3
    assert [t["Feature"] for t in got] == ["ENST_A", "ENST_B", "ENST_C"]


def test_coding_consequence_behind_non_coding_entries_is_reachable():
    """The regression: the informative entry is not first.

    Mirrors the NRAS/CSDE1 record -- many non-coding entries for an overlapping
    gene, with the driver's missense annotation last.
    """
    csq = [entry("A", "downstream_gene_variant", "CSDE1", f"ENST_{i}") for i in range(46)]
    csq.append(entry("A", "missense_variant", "NRAS", "ENST_NRAS", aa="Q/H", pos="61"))

    got = parser().parse_csq_entries(csq, CSQ_FORMAT, "A")

    assert len(got) == 47
    nras = [t for t in got if t["SYMBOL"] == "NRAS"]
    assert len(nras) == 1, "the driver's annotation was discarded"
    assert nras[0]["Amino_acids"] == "Q/H"
    assert nras[0]["Protein_position"] == "61"


def test_entries_for_other_alleles_are_excluded():
    """Allele filtering still applies -- a multiallelic record must not leak."""
    csq = [
        entry("A", "missense_variant", "GENE1", "ENST_A", aa="P/L", pos="10"),
        entry("T", "missense_variant", "GENE2", "ENST_T", aa="R/W", pos="20"),
        entry("A", "missense_variant", "GENE3", "ENST_A2", aa="G/V", pos="30"),
    ]
    got = parser().parse_csq_entries(csq, CSQ_FORMAT, "A")

    assert [t["SYMBOL"] for t in got] == ["GENE1", "GENE3"]
    assert all(t["Allele"] == "A" for t in got)


def test_no_matching_allele_returns_empty():
    csq = [entry("T", "missense_variant", "GENE1", "ENST_T", aa="P/L", pos="10")]
    assert parser().parse_csq_entries(csq, CSQ_FORMAT, "A") == []


def test_fields_are_mapped_onto_the_format_keys():
    csq = [entry("A", "missense_variant", "NRAS", "ENST_C", aa="Q/H", pos="61")]
    got = parser().parse_csq_entries(csq, CSQ_FORMAT, "A")

    assert got[0] == {
        "Allele": "A",
        "Consequence": "missense_variant",
        "IMPACT": "MODERATE",
        "SYMBOL": "NRAS",
        "Gene": "ENSG",
        "Feature": "ENST_C",
        "Amino_acids": "Q/H",
        "Protein_position": "61",
    }
