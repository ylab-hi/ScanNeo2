"""Cohort summary of the per-source prioritization results (#227).

Writes one row per (sample, source) that ran: how many input records the
source had, how many variant effects and neoepitopes came out of them, and a
status naming the first stage that came out empty. An empty source is often
legitimate (a sample without exitrons), but can also be the symptom of an
upstream failure; this makes it visible without failing the run.

Run by the `summarize` rule:

    summarize.py --output results/summary.tsv --classes I [II] \
        --entries '[{"sample": ..., "source": ..., "dir": ..., "inputs": [...]}, ...]'
"""

import argparse
import csv
import gzip
import json
import os

def count_input_records(path):
    """Records in a source input: VCF body lines, or TSV data rows.

    VCF headers start with '#'. The fusion TSV's header does too; the custom
    protein TSV has a plain header row, which is skipped as the first line.
    """
    opener = gzip.open if path.endswith(".gz") else open
    vcf = path.endswith((".vcf", ".vcf.gz"))
    n = 0
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if not line.strip() or line.startswith("#"):
                continue
            if not vcf and i == 0:
                continue
            n += 1
    return n


def count_table(path, column=None):
    """(data rows, distinct values of `column`) of a tab-separated table."""
    if not os.path.exists(path):
        return 0, 0
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = 0
        distinct = set()
        for row in reader:
            rows += 1
            if column:
                distinct.add(row[column])
    return rows, len(distinct)


def status(input_records, variant_effects, neoepitopes):
    """The first stage that came out empty, or 'ok'."""
    if input_records == 0:
        return "no_input"
    if variant_effects == 0:
        return "no_effects"
    if sum(neoepitopes) == 0:
        return "no_neoepitopes"
    return "ok"


def summarize_source(sample, source, source_dir, inputs, classes):
    effects, _ = count_table(os.path.join(source_dir, f"{source}_variant_effects.tsv"))
    row = {
        "sample": sample,
        "source": source,
        "input_records": sum(count_input_records(p) for p in inputs),
        "variant_effects": effects,
    }
    for cls in classes:
        rows, peptides = count_table(
            os.path.join(source_dir, f"{source}_mhc-{cls}_neoepitopes.txt"),
            "mt_epitope_seq",
        )
        row[f"neoepitopes_mhc-{cls}"] = rows
        row[f"distinct_peptides_mhc-{cls}"] = peptides
    row["status"] = status(
        row["input_records"],
        effects,
        [row[f"neoepitopes_mhc-{cls}"] for cls in classes],
    )
    return row


def write_summary(entries, classes, out_path):
    """entries: dicts with sample, source, dir and inputs, in output order."""
    columns = ["sample", "source", "input_records", "variant_effects"]
    for cls in classes:
        columns += [f"neoepitopes_mhc-{cls}", f"distinct_peptides_mhc-{cls}"]
    columns.append("status")
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for e in entries:
            writer.writerow(
                summarize_source(e["sample"], e["source"], e["dir"], e["inputs"], classes)
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--output", required=True)
    p.add_argument("--classes", nargs="+", required=True)
    p.add_argument("--entries", required=True, help="JSON list of source entries")
    args = p.parse_args()
    write_summary(json.loads(args.entries), args.classes, args.output)
