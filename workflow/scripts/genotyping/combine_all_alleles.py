import os
import sys

"""
This script combines predicted and user-provided mhc-II alleles into a
single file and also compares them with alleles from the the reference

https://help.iedb.org/hc/en-us/articles/114094151851-HLA-allele-frequencies-and-reference-sets-with-maximal-population-coverage

Usage:

python combine_all_mhcII_alleles.py '<input>' <class> <output>


class = 'mhc-I' or 'mhc-II'

"""

def parse_refset(mhc_class):
    refset = []
    with open(f"workflow/scripts/genotyping/{mhc_class}_refset.txt", "r") as fh:
        for line in fh:
            refset.append(line.rstrip())
    return refset


def main():
    infiles = sys.argv[1].split(' ')
    mhc_class = sys.argv[2]
    refset = parse_refset(mhc_class)

    alleles = {}
    per_file_stats = []  # (path, size_bytes, lines_read, matched_alleles)

    for infile in infiles:
        lines_read = 0
        matched = 0
        size_bytes = os.path.getsize(infile)

        with open(infile, "r") as fh_in:
            for line in fh_in:
                lines_read += 1
                cols = line.rstrip().split("\t")
                if len(cols) != 2:
                    print(f"Invalid input file: {infile}")
                    sys.exit(1)

                source = cols[0]
                mhc = cols[1]

                # chop down the alleles to first two fields
                # e.g., HLA-A*02:01:01 becomes HLA-A*02:01
                if len(mhc.split(':')) > 2:
                    mhc = mhc.split(':')[0] + ':' + mhc.split(':')[1]


                if mhc in refset:
                    matched += 1
                    if mhc not in alleles:
                        alleles[mhc] = source
                    else:
                        alleles[mhc] += ',' + source

        per_file_stats.append((infile, size_bytes, lines_read, matched))

    outfile = sys.argv[3]
    with open(outfile, 'w') as fh_out:
        for allele in alleles:
            fh_out.write(f'{alleles[allele]}\t{allele}\n')

    if len(alleles) == 0:
        print("ERROR: No valid alleles collected from any input source.\n")
        print("Per-source breakdown:")
        for path, size, lines, matched in per_file_stats:
            if size == 0:
                state = "EMPTY (0 bytes) — upstream rule produced no output"
            elif lines == 0:
                state = "no content rows"
            elif matched == 0:
                state = (f"{lines} row(s) read, 0 matched the {mhc_class} "
                         f"reference set ({len(refset)} alleles)")
            else:
                # Unreachable here (len(alleles) == 0 implies no file matched),
                # but keep a sensible fallback.
                state = f"{lines} row(s), {matched} matched allele(s)"
            print(f"  {path} — {state}")

        print("\nLikely upstream culprit:")
        print("  For predicted alleles: check the HLA-filtered BAMs")
        print(f"    results/<sample>/hla/{mhc_class}/reads/<group>_<nartype>_flt_{{SE,PE}}.bam")
        print("    - If near-empty: the yara HLA-panel filter discarded most reads")
        print("      (low HLA-region coverage in input FASTQ/BAM, or off-target sequencing)")
        print("    - If non-empty but per-batch OptiType results are 0 bytes: the")
        print("      optitype wrapper hit its '<10 reads per batch' empty-output branch")
        print(f"      (see results/<sample>/hla/{mhc_class}/genotyping/<group>_<nartype>_flt_*/)")
        print("  For user-provided alleles: check the sample sheet's custom_hla_I / custom_hla_II columns")
        sys.exit(1)

main()
