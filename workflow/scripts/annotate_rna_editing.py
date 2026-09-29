"""Annotate RNA somatic SNVs with REDIportal A-to-I editing evidence.

For each SNV whose substitution matches the A-to-I signature (genomic A>G on
plus-strand genes, T>C on minus-strand genes), look up the position in the
tabix-indexed REDIportal TABLE1 atlas. If a REDIportal site at that position
carries the same Ref>Ed substitution, add INFO flags describing the edit. An
MNP is checked base by base and flagged only if every changed base is one:

  RE          known A-to-I editing site (REDIportal)
  RE_nTissues number of normal GTEx tissues edited there (constitutiveness)
  RE_exonic   the edit is exonic/nonsynonymous (REDIportal ANNOVAR annotation)
  RE_TCGA     number of TCGA tumor samples with the edit

The downstream split (issue #186, phase 3) uses these to separate point
mutations from editing and to drop constitutive (self) edits.

Usage: annotate_rna_editing.py <in.vcf.gz> <rediportal.txt.gz> <out.vcf.gz>
"""

import sys

import pysam

# REDIportal TABLE1 columns (0-based) after tabix-skipped header.
REGION, POSITION, REF, ED = 1, 2, 3, 4
EXONICFUNC, NTISSUES, NTCGA = 12, 24, 34

# genomic representation of an A-to-I edit: A>G (+ strand), T>C (- strand)
AI_SIGNATURE = {("A", "G"), ("T", "C")}


def editing_sites(pos, ref, alt, site_at):
    """REDIportal rows for every base an SNV or MNP changes, or None.

    Mutect2 merges adjacent substitutions on one haplotype into an MNP, and
    hyper-edited regions edit neighbouring adenosines on the same reads, so a
    cluster of known edits arrives as one record (AA>GG). It counts as editing
    only if every changed base is an A-to-I substitution that REDIportal lists
    at that position; one base that is not makes the record a candidate
    mutation. site_at(pos, ref, ed) returns the matching row or None.
    """
    if len(alt) != len(ref):
        return None
    changed = [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]
    if not changed or any((r, a) not in AI_SIGNATURE for _, r, a in changed):
        return None
    rows = [site_at(p, r, a) for p, r, a in changed]
    return None if None in rows else rows


def main():
    in_vcf, redi_path, out_vcf = sys.argv[1], sys.argv[2], sys.argv[3]

    redi = pysam.TabixFile(redi_path)
    redi_contigs = set(redi.contigs)

    vcf_in = pysam.VariantFile(in_vcf)
    header = vcf_in.header
    header.info.add("RE", 0, "Flag", "Known A-to-I RNA editing site (REDIportal)")
    header.info.add(
        "RE_nTissues",
        1,
        "Integer",
        "Normal GTEx tissues edited at this site (REDIportal nTissues)",
    )
    header.info.add(
        "RE_exonic", 0, "Flag", "Editing site is exonic/nonsynonymous (REDIportal)"
    )
    header.info.add(
        "RE_TCGA",
        1,
        "Integer",
        "TCGA tumor samples with the edit (REDIportal nTCGASamples)",
    )
    # "wz" = BGZF-compressed, matching the .vcf.gz output ("w" would write plain text)
    vcf_out = pysam.VariantFile(out_vcf, "wz", header=header)

    annotated = 0
    for record in vcf_in:
        if record.contig in redi_contigs:

            def site_at(pos, ref, ed):
                for line in redi.fetch(record.contig, pos - 1, pos):
                    fields = line.split("\t")
                    if int(fields[POSITION]) == pos and fields[REF] == ref and fields[ED] == ed:
                        return fields
                return None

            # Mutect2 SNVs are effectively biallelic, but check every ALT so a
            # multiallelic site with an A-to-I allele is still annotated.
            for alt in record.alts or ():
                rows = editing_sites(record.pos, record.ref, alt, site_at)
                if rows is None:
                    continue
                record.info["RE"] = True
                # an MNP is as constitutive, and as widespread in TCGA, as its
                # least-edited base; exonic if any of its bases is
                ntissues = [int(f[NTISSUES]) for f in rows if f[NTISSUES].isdigit()]
                if ntissues:
                    record.info["RE_nTissues"] = min(ntissues)
                if any(len(f) > EXONICFUNC and "nonsynonymous" in f[EXONICFUNC] for f in rows):
                    record.info["RE_exonic"] = True
                ntcga = [int(f[NTCGA]) for f in rows if len(f) > NTCGA and f[NTCGA].isdigit()]
                if ntcga:
                    record.info["RE_TCGA"] = min(ntcga)
                annotated += 1
                break
        vcf_out.write(record)

    vcf_out.close()
    vcf_in.close()
    sys.stderr.write(f"annotated {annotated} SNVs/MNPs as known A-to-I editing sites\n")


if __name__ == "__main__":
    main()
