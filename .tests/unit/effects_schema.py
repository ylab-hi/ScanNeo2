"""The variant-effects table as effects.py writes it.

prediction.py looks its columns up by name, so a fixture that invents header
names exercises nothing: the lookup fails on the first row it reads. Holding the
schema in one place keeps the prediction tests from disagreeing about it, and
builds each row against the column order rather than by counting fields.

The names mirror VariantEffects.write_header. They are not imported from there
because effects.py pulls in reference.py and pyfaidx, which the test
environment does not install.
"""

from pathlib import Path

COLUMNS = [
    "chrom", "start", "end", "gene_id", "gene_name", "transcript_id",
    "source", "group", "seqtype", "var_type", "wt_subseq", "mt_subseq",
    "var_start", "aa_var_start", "aa_var_end", "vaf", "ao", "dp", "TPM",
    "NMD", "PTC_dist_ejc", "PTC_exon_number", "NMD_escape_rule",
]

HEADER = "\t".join(COLUMNS)


def effects_row(transcript_id, wt_subseq, mt_subseq, aa_var_start, aa_var_end,
                var_type="SNV", seqtype="dnaseq"):
    """One table row. Columns the prediction tests never vary are fixed here."""
    values = {
        "chrom": "chr1", "start": "100", "end": "101",
        "gene_id": "G1", "gene_name": "GENE1",
        "transcript_id": transcript_id,
        "source": "DNA", "group": "g1", "seqtype": seqtype,
        "var_type": var_type,
        "wt_subseq": wt_subseq, "mt_subseq": mt_subseq,
        "var_start": "100",
        "aa_var_start": str(aa_var_start), "aa_var_end": str(aa_var_end),
        "vaf": "0.5", "ao": "10", "dp": "20", "TPM": "1.0",
        "NMD": "NA", "PTC_dist_ejc": "NA", "PTC_exon_number": "NA",
        "NMD_escape_rule": "NA",
    }
    return "\t".join(values[name] for name in COLUMNS)


def read_table(path):
    """Rows of a written table as dicts, keyed by its own header line."""
    lines = Path(path).read_text().splitlines()
    names = lines[0].split("\t")
    return [dict(zip(names, line.split("\t"))) for line in lines[1:]]
