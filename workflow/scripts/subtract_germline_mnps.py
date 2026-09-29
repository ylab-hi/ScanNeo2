"""Drop RNA MNPs whose every changed base is a matched-normal germline SNV.

Tumor-only RNA Mutect2 merges adjacent substitutions on one haplotype into an
MNP (GG>AT), while the germline reference comes from HaplotypeCaller, which
emits the same pair as two SNV records. The exact-allele subtraction
(bcftools isec) therefore never matches such an MNP, and a phased pair of
germline SNVs would pass as a somatic call. An MNP is dropped here when every
base it changes is a germline SNV with the same alternate base; one base that
is not keeps it as a candidate somatic call.

Usage: subtract_germline_mnps.py <in.vcf.gz> <germline.vcf.gz> <out.vcf.gz>
"""

import sys

import pysam


def fully_germline(pos, ref, alts, is_germline):
    """True if every ALT is an MNP whose changed bases are all germline SNVs.

    is_germline(pos, ref, alt) answers for one base. A multiallelic record is
    dropped only if all of its ALTs are, so a somatic allele is never lost
    with a germline one.
    """
    if len(ref) < 2 or not alts:
        return False
    for alt in alts:
        if len(alt) != len(ref):
            return False
        changed = [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]
        if not changed or not all(is_germline(p, r, a) for p, r, a in changed):
            return False
    return True


def main():
    in_vcf, germline_vcf, out_vcf = sys.argv[1], sys.argv[2], sys.argv[3]

    germline = pysam.VariantFile(germline_vcf)
    germline_contigs = set(germline.header.contigs)
    vcf_in = pysam.VariantFile(in_vcf)
    # "wz" = BGZF-compressed, matching the .vcf.gz output
    vcf_out = pysam.VariantFile(out_vcf, "wz", header=vcf_in.header)

    dropped = 0
    for record in vcf_in:
        if record.contig in germline_contigs:

            def is_germline(pos, ref, alt):
                return any(
                    g.pos == pos and g.ref == ref and alt in (g.alts or ())
                    for g in germline.fetch(record.contig, pos - 1, pos)
                )

            if fully_germline(record.pos, record.ref, record.alts, is_germline):
                dropped += 1
                continue
        vcf_out.write(record)

    vcf_out.close()
    vcf_in.close()
    sys.stderr.write(f"dropped {dropped} MNPs made entirely of germline SNVs\n")


if __name__ == "__main__":
    main()
