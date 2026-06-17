"""Process user-supplied (wildtype, mutant) protein pairs.

Reads a TSV and registers each row as a variant effect, bypassing variant
calling and VEP annotation entirely. The output table's protein-derived columns
(`wt_subseq`, `mt_subseq`, `aa_var_start`, `aa_var_end`) are populated by
`VariantEffects.change_entry` (which diffs wt against mt to locate the variant
region); the genomic / transcript / expression columns are left empty unless
the TSV supplies them.

Input format
============
A header line is required. Required columns:

    id, wildtype_protein, mutant_protein

Optional columns (each defaults to empty / -1 if absent):

    vaf, ao, dp, gene_id, gene_name, transcript_id, chrom, group, var_type

The wildtype protein is mandatory per row — without it there is no variant
region to detect, no wildtype contrast for the binding-affinity comparison, no
agretopicity, and no self-similarity check. Mutant-only rows are rejected
with a hard error.
"""

import csv
import sys
from pathlib import Path

import effects

REQUIRED_COLS = ("id", "wildtype_protein", "mutant_protein")
OPTIONAL_COLS = (
    "vaf", "ao", "dp",
    "gene_id", "gene_name", "transcript_id", "chrom", "group", "var_type",
)


class Proteins:
    def __init__(self, proteins_input, options, vartype):
        with effects.VariantEffects(options, vartype) as variant_effects:
            self.variant_effects = variant_effects

            with open(proteins_input, "r", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                self.validate_header(reader.fieldnames, proteins_input)

                for row_num, row in enumerate(reader, start=2):
                    self.process_row(row, row_num, proteins_input)

    @staticmethod
    def validate_header(fieldnames, source):
        if not fieldnames:
            sys.exit(
                f"[proteins error] {source}: empty file or missing header. "
                f"Required columns: {', '.join(REQUIRED_COLS)}."
            )
        missing = [c for c in REQUIRED_COLS if c not in fieldnames]
        if missing:
            sys.exit(
                f"[proteins error] {source}: missing required column(s): "
                f"{', '.join(missing)}. Required columns: {', '.join(REQUIRED_COLS)}."
            )

    def process_row(self, row, row_num, source):
        row_id = row.get("id", "").strip()
        wt_seq = (row.get("wildtype_protein") or "").strip()
        mt_seq = (row.get("mutant_protein") or "").strip()

        if not wt_seq or not mt_seq:
            sys.exit(
                f"[proteins error] {source}:{row_num} (id={row_id!r}): "
                f"both wildtype_protein and mutant_protein must be non-empty. "
                f"Mutant-only rows have no variant region to detect and are rejected."
            )

        # Truncate at stop-codon markers — same convention as the VEP-derived
        # path uses on `wt_seq`/`mt_seq`. Without this, an X (or *) reaching the
        # IEDB binding-affinity tool is rejected and the whole prediction batch
        # is dropped, taking neighbouring valid rows with it.
        wt_seq = self.truncate_at_stop(wt_seq)
        mt_seq = self.truncate_at_stop(mt_seq)

        if not wt_seq or not mt_seq:
            sys.exit(
                f"[proteins error] {source}:{row_num} (id={row_id!r}): "
                f"truncation at stop-codon marker ('*' or 'X') left an empty "
                f"sequence. There is no detectable variant region."
            )

        vaf = self.as_float(row.get("vaf"), default=-1.0)
        ao = self.as_int(row.get("ao"), default=-1)
        dp = self.as_int(row.get("dp"), default=-1)

        self.variant_effects.change_entry(
            chrom=row.get("chrom") or None,
            start=None,
            end=None,
            gene_id=row.get("gene_id") or None,
            gene_name=row.get("gene_name") or None,
            transcript_id=row.get("transcript_id") or None,
            transcript=None,
            transcript_bp=None,
            source="custom_protein",
            group=(row.get("group") or "custom_protein"),
            # default to "custom" so determine_NMD's "frameshift" check skips
            var_type=(row.get("var_type") or "custom"),
            # var_start=0 — determine_var_bnds scans wt vs mt from this index
            # to locate the actual variant region
            var_start=0,
            wt_seq=wt_seq,
            mt_seq=mt_seq,
            vaf=vaf,
            ao=ao,
            dp=dp,
            nmd_event=None,
        )

        if self.variant_effects.self_dissimilarity():
            self.variant_effects.write_entry()

    @staticmethod
    def truncate_at_stop(sequence):
        """Return the prefix of `sequence` before the first stop-codon marker
        (`*` or `X`). Mirrors `variants.py:Variants.scan_stop_codon` so the
        protein source produces sequences that flow cleanly into the same
        binding-affinity prediction step."""
        for marker in ("*", "X"):
            sep = sequence.find(marker)
            if sep != -1:
                sequence = sequence[:sep]
        return sequence

    @staticmethod
    def as_float(value, default):
        if value is None or value == "":
            return default
        try:
            return float(value)
        except ValueError:
            return default

    @staticmethod
    def as_int(value, default):
        if value is None or value == "":
            return default
        try:
            return int(value)
        except ValueError:
            return default
