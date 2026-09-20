"""Annotate RNA somatic SNVs with REDIportal A-to-I editing evidence.

For each SNV whose substitution matches the A-to-I signature (genomic A>G on
plus-strand genes, T>C on minus-strand genes), look up the position in the
tabix-indexed REDIportal TABLE1 atlas. If a REDIportal site at that position
carries the same Ref>Ed substitution, add INFO flags describing the edit:

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
        # Mutect2 SNVs are effectively biallelic, but check every ALT so a
        # multiallelic site with an A-to-I allele is still annotated.
        if len(record.ref) == 1 and record.contig in redi_contigs:
            matched = False
            for alt in record.alts or ():
                if len(alt) != 1 or (record.ref, alt) not in AI_SIGNATURE:
                    continue
                for line in redi.fetch(record.contig, record.pos - 1, record.pos):
                    fields = line.split("\t")
                    if (
                        int(fields[POSITION]) == record.pos
                        and fields[REF] == record.ref
                        and fields[ED] == alt
                    ):
                        record.info["RE"] = True
                        if fields[NTISSUES].isdigit():
                            record.info["RE_nTissues"] = int(fields[NTISSUES])
                        if len(fields) > EXONICFUNC and "nonsynonymous" in fields[EXONICFUNC]:
                            record.info["RE_exonic"] = True
                        if len(fields) > NTCGA and fields[NTCGA].isdigit():
                            record.info["RE_TCGA"] = int(fields[NTCGA])
                        annotated += 1
                        matched = True
                        break
                if matched:
                    break
        vcf_out.write(record)

    vcf_out.close()
    vcf_in.close()
    sys.stderr.write(f"annotated {annotated} SNVs as known A-to-I editing sites\n")


if __name__ == "__main__":
    main()
