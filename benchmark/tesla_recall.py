"""Score ScanNeo2's neoepitope output against the TESLA validated peptides.

Reads both TESLA supplementary tables (see benchmark/TESLA/README.md) and, for every
pMHC they tested, reports whether ScanNeo2 predicted it.

The two classes are scored in opposite directions. The immunogenic pMHC are the
positives: recall over them is the benchmark, and higher is better. The
non-immunogenic pMHC were tested and found NOT immunogenic, so they are confirmed
true negatives -- proposing one is a false positive, and a lower number is better.
Retaining few of them means they are already gone before any prioritisation filter
runs, which is the desired behaviour, not a sensitivity gap.

A pMHC is counted as found when the mutant epitope sequence AND the HLA allele both
match a row of results/<sample>/prioritization/mhc-I_neoepitopes_all.txt. Peptide-only
matches are reported too, since the two differ when HLA typing disagrees.

Usage: python benchmark/tesla_recall.py [--export <out.tsv>]
"""

import argparse
import csv
import sys
from pathlib import Path

TABLES = [
    "benchmark/TESLA/tesla_main_cohort_tableS4.tsv",
    "benchmark/TESLA/tesla_validation_cohort_tableS7.tsv",
]

# ---------------------------------------------------------------------------
# EDIT THIS FOR YOUR OWN RUN.
#
# TESLA subject id -> the directory under results/ holding that subject's
# output, i.e. where results/<sample>/prioritization/mhc-I_neoepitopes_all.txt
# lives. The values below are the names used for the run the figures in
# benchmark/TESLA/README.md were measured on.
#
# The subject ids on the left are fixed by the published tables. The names on
# the right come from whatever sample sheet configured the run and cannot be
# inferred from this repository or from results/, so they have to be set by
# hand. Delete any subject whose data you do not have: the denominators are
# derived from whatever is listed here, so a partial cohort still scores
# correctly against its own totals.
#
# Subject 10 is deliberately absent. It has 3 immunogenic pMHC, but its WES/RNA
# is not available via Synapse, so no pipeline can score it; listing it would
# only inflate the denominator with unreachable peptides.
# ---------------------------------------------------------------------------
SAMPLES = {
    "1": "melanoma_patient1",
    "2": "melanoma_patient2",
    "3": "melanoma_patient3",
    "4": "melanoma_patient4",
    "8": "melanoma_patient8",
    "9": "melanoma_patient9",
    "12": "lung_patient12",
    "16": "lung_patient16",
}

# S4 writes VALIDATED as TRUE/FALSE, S7 as 1/0.
TRUTHY = {"TRUE", "1"}


def predictions_path(sample):
    return Path(f"results/{sample}/prioritization/mhc-I_neoepitopes_all.txt")


def read_tesla(paths):
    """Assessed pMHC as {subject: [(allele, peptide, validated), ...]}."""
    out = {}
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                subject = row["PATIENT_ID"].strip()
                if not subject:
                    continue
                # the allele is the PMHC prefix in both schemas; S7 has no MHC column
                allele = row["PMHC"].split("_", 1)[0]
                out.setdefault(subject, []).append(
                    (allele, row["ALT_EPI_SEQ"].strip(),
                     row["VALIDATED"].strip().upper() in TRUTHY)
                )
    return out


def read_predicted(sample):
    """Our predictions as {(allele, peptide): [row, ...]}, allele stripped of 'HLA-'."""
    hits = {}
    with open(predictions_path(sample), newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            allele = row["allele"].removeprefix("HLA-")
            hits.setdefault((allele, row["mt_epitope_seq"]), []).append(row)
    peptides = {pep for _, pep in hits}
    return hits, peptides


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", help="write the matching prediction rows here")
    args = ap.parse_args()

    subjects = SAMPLES
    tesla = read_tesla(TABLES)

    # fail before reporting anything, so a partial cohort cannot look like a low score
    problems = [f"subject {s} is not in the TESLA tables" for s in subjects
                if s not in tesla]
    problems += [f"subject {s}: no predictions at {predictions_path(d)}"
                 for s, d in subjects.items() if not predictions_path(d).exists()]
    if problems:
        problems.append("edit SAMPLES at the top of this script for your own run")
        sys.exit("\n".join(problems))

    export_rows, fieldnames = [], None
    totals = {k: 0 for k in ("v_found", "v_tot", "n_found", "n_tot", "pep_found")}

    assessed_total = sum(len(tesla[s]) for s in subjects)
    immunogenic_total = sum(1 for s in subjects for *_, v in tesla[s] if v)
    print(f"{len(subjects)} subjects from SAMPLES: "
          f"{assessed_total} assessed pMHC, {immunogenic_total} immunogenic")
    print()
    print(f"{'subject':<9} {'sample':<20} {'immunogenic':>12} {'true neg':>12} {'assessed':>12}")
    print(f"{'':<9} {'':<20} {'(want high)':>12} {'(want low)':>12} {'':>12}")
    print("-" * 70)
    for subject in sorted(subjects, key=int):
        sample = subjects[subject]
        assessed = tesla[subject]
        hits, peptides = read_predicted(sample)

        counts = {k: 0 for k in ("v_found", "v_tot", "n_found", "n_tot", "pep_found")}
        for allele, peptide, validated in assessed:
            key = "v" if validated else "n"
            counts[f"{key}_tot"] += 1
            if peptide in peptides:
                counts["pep_found"] += 1
            if (allele, peptide) in hits:
                counts[f"{key}_found"] += 1
                for row in hits[(allele, peptide)]:
                    if fieldnames is None:
                        fieldnames = ["tesla_subject", "tesla_validated"] + list(row)
                    export_rows.append(
                        {"tesla_subject": subject,
                         "tesla_validated": "yes" if validated else "no", **row}
                    )
        for k in totals:
            totals[k] += counts[k]

        a_found = counts["v_found"] + counts["n_found"]
        a_tot = counts["v_tot"] + counts["n_tot"]
        print(
            f"{subject:<9} {sample:<20} "
            f"{counts['v_found']:>5}/{counts['v_tot']:<6} "
            f"{counts['n_found']:>5}/{counts['n_tot']:<6} "
            f"{a_found:>5}/{a_tot:<6}"
        )

    print("-" * 70)
    a_found = totals["v_found"] + totals["n_found"]
    a_tot = totals["v_tot"] + totals["n_tot"]
    print(
        f"{'TOTAL':<30} {totals['v_found']:>5}/{totals['v_tot']:<6} "
        f"{totals['n_found']:>5}/{totals['n_tot']:<6} {a_found:>5}/{a_tot:<6}"
    )
    print(f"\nimmunogenic recall      : {totals['v_found']}/{totals['v_tot']}"
          f" ({100 * totals['v_found'] / totals['v_tot']:.1f}%)   higher is better")
    print(f"true negatives retained : {totals['n_found']}/{totals['n_tot']}"
          f" ({100 * totals['n_found'] / totals['n_tot']:.1f}%)   lower is better")
    print(f"true negatives excluded : {totals['n_tot'] - totals['n_found']}"
          f" of {totals['n_tot']}")
    print(f"peptide-only match      : {totals['pep_found']}/{a_tot}"
          "   (ignores the HLA allele)")

    # every listed subject is in the tables and has predictions, both checked above,
    # so the scored denominators must equal the ones derived from the mapping
    assert a_tot == assessed_total, f"{a_tot} scored vs {assessed_total} assessed"
    assert totals["v_tot"] == immunogenic_total, (
        f"{totals['v_tot']} scored vs {immunogenic_total} immunogenic")

    if args.export:
        with open(args.export, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(export_rows)
        print(f"\nexported {len(export_rows)} prediction rows -> {args.export}")


if __name__ == "__main__":
    sys.exit(main())
