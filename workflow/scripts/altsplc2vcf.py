import sys
import os
import re
import vcfpy
import configargparse
from pyfaidx import Fasta
import gzip 
from pathlib import Path

"""
    python altsplc_to_variants.py '<input_txts>' <genome_fasta>
"""

def add_info_field(header, key, number, dtype, desc):
    info = vcfpy.OrderedDict([
        ('ID', key),
        ('Number', number),
        ('Type', dtype),
        ('Description', desc)
    ])
    header.add_info_line(info)
    return header

def create_header(samplename, ref):
    samples = vcfpy.SamplesInfos([samplename])
    # extract samples from filename
    header = vcfpy.Header(
            lines = [
                    vcfpy.HeaderLine(key="fileformat", value="VCFv4.2"),
            ],
            samples=samples)
    
    header = add_info_field(header, 'SVTYPE', '1', 'String', 'Type of structural variant')
    header = add_info_field(header, 'STRAND', '1', 'String', 'Junction Strand')
    header = add_info_field(header, 'END', '1', 'Integer', 'End position of the structural variant')
    header = add_info_field(header, 'SVLEN', '1', 'Integer', 'Length of the structural variant')
    header = add_info_field(header, 'AO', '1' , 'Integer', 'Number of reads supporting the alt splicing')
    header = add_info_field(header, 'PSI', '1', 'Float', 'Percentage spliced-in')
    header = add_info_field(header, 'GeneName', '1', 'String', 'Gene name')
    header = add_info_field(header, 'GRP', '1', 'String', 'Name of the group')
    header = add_info_field(header, 'SRC', '1', 'String', 'Source of the variant')
    # SplAdder runs on RNA only, so this is constant here; it exists so every
    # source writes the same provenance fields (see issue #244)
    header = add_info_field(header, 'SEQ', '1', 'String', 'Sequencing type the variant was called from')
            
    header.add_format_line(vcfpy.OrderedDict(
        [('ID','GT'),
         ('Number',1),
         ('Type', 'String'),
         ('Description', 'Genotype')
         ]))
    
    # add contig lines
    for chrom in ref.keys():
        header.add_contig_line(vcfpy.OrderedDict(
            [('ID', chrom),
             ('length', len(ref[chrom]))]
        ))


    return header


def main():
    options = parse_arguments()
    seq = Fasta(options.reference)
    
    header = create_header(options.group, seq)
    writer = vcfpy.Writer.from_path(options.output, header)

    # different spladder events can imply the same genomic change -- an exon
    # skipped on its own and again as part of a run, say -- and each occurrence
    # would otherwise be written again and counted as a further candidate
    seen = set()

    def write_once(rec):
        key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0].value)
        if key in seen:
            return
        seen.add(key)
        writer.write_record(rec)

    # sorted so that which event type keeps a shared allele is fixed, rather
    # than following the order the filesystem happens to return
    for root, dirs, files in sorted(os.walk(options.input)):
        for file in sorted(files):
            if file.endswith(".txt.gz"):
                print(file)
                # determine type of splicing from file
                pattern = re.compile(r'merge_graphs_(\w+)_C[1-3].confirmed.txt.gz')
                match = pattern.search(file)
                if match:
                    event_type = match.group(1)

                    fh = gzip.open(os.path.join(root, file), 'r')
                    next(fh)
                    for line in fh:
                        column = line.decode().rstrip().split('\t')

                        if (event_type == "alt_5prime" or
                            event_type == "alt_3prime" or
                            event_type == "exon_skip" or
                            event_type == "intron_retention" or
                            event_type == "mult_exon_skip"):

                            try:
                                chrom = column[0]
                                strand = column[1]
                                event_id = column[2]

                                # https://github.com/ratschlab/spladder/issues/168
                                annotated = column[3]

                                gene_name = column[4]

                                if event_type == "mult_exon_skip":
                                    # spladder lists every skipped exon here,
                                    # colon separated, so these two columns hold
                                    # a list rather than a single coordinate
                                    starts = [int(v) for v in column[7].split(':')]
                                    ends = [int(v) for v in column[8].split(':')]
                                    start = min(starts)
                                    end = max(ends)
                                else:
                                    start = int(column[7])
                                    end = int(column[8])

                            # in some cases lines are incomplete (faulty spladder output)
                            except (IndexError, ValueError):
                                print("malformed event - skipping event")
                                continue

                            # every record anchors on the base before the event,
                            # so an event starting at the first base of a
                            # chromosome has no anchor to place it on. The
                            # negative index would wrap to the chromosome's end
                            # and build the allele from unrelated sequence.
                            if start < 2:
                                print("event at chromosome start - skipping event")
                                continue

                            if (event_type == "alt_5prime" or
                                event_type == "alt_3prime" or
                                event_type == "intron_retention"):

                                sv_type = "INS"
                                sv_len = end - start + 1

                                # POS is start-1, the base before the event, so
                                # the alleles must begin there too. pyfaidx
                                # indexes from 0, which puts that base at
                                # start-2 and the event's own first base at
                                # start-1.
                                ref = str(seq[chrom][start-2:start-1]).upper()
                                alt = (str(seq[chrom][start-2:start-1]) + str(seq[chrom][start-1:end])).upper()
                            
                                ao = column[12]

                                if event_type == "intron_retention":
                                    psi = column[15]
                                else:
                                    psi = column[16]

                                rec = vcfpy.Record(
                                        CHROM=chrom,
                                        POS=start-1,
                                        ID=[event_id],
                                        REF=ref,
                                        ALT=[vcfpy.Substitution(type_=sv_type, value=alt)],
                                        QUAL='.',
                                        FILTER=['PASS'],
                                        INFO = vcfpy.OrderedDict(
                                            [('SVTYPE', sv_type),
                                             ('END', end),
                                             ('AO', ao),
                                             ('STRAND', strand),
                                             ('GeneName', gene_name),
                                             ('SVLEN', sv_len),
                                             ('PSI', psi),
                                             ('GRP', options.group),
                                             ('SRC', event_type),
                                             ('SEQ', 'rnaseq')]),

                                        FORMAT=['GT'],
                                        calls = [vcfpy.Call(sample=options.group, data={'GT': '0/1'})]
                                )
                                write_once(rec)

                            elif (event_type == "exon_skip" or
                                  event_type == "mult_exon_skip"):
                                sv_type = "DEL"
                                sv_len = -(end - start + 1)

                                """ for some reason spladder reports exon skipping
                                events in which the start and end position of the
                                skipped exon are the same... skip those events """
                                if start == end:
                                    continue

                                # anchored on start-1 like the INS types above
                                ref = (str(seq[chrom][start-2:start-1]) + str(seq[chrom][start-1:end])).upper()
                                alt = str(seq[chrom][start-2:start-1]).upper()

                                if event_type == "mult_exon_skip":
                                    # the deleted span reaches from the first
                                    # skipped exon to the last, so it also covers
                                    # the introns between them. Those are spliced
                                    # out either way, and both flanking introns
                                    # survive, so the transcript still joins e1 to
                                    # e3 -- the isoform spladder observed.
                                    ao = column[18]
                                    psi = column[19]
                                else:
                                    ao = column[16]
                                    psi = column[17]

                                rec = vcfpy.Record(
                                        CHROM=chrom,
                                        POS=start-1,
                                        ID=[event_id],
                                        REF=ref,
                                        ALT=[vcfpy.Substitution(type_=sv_type, value=alt)],
                                        QUAL='.',
                                        FILTER=['PASS'],
                                        INFO = vcfpy.OrderedDict(
                                            [('SVTYPE', sv_type),
                                             ('END', end),
                                             ('AO', ao),
                                             ('STRAND', strand),
                                             ('GeneName', gene_name),
                                             ('SVLEN', sv_len),
                                             ('PSI', psi),
                                             ('GRP', options.group),
                                             ('SRC', event_type),
                                             ('SEQ', 'rnaseq')]),

                                        FORMAT=['GT'],
                                        calls = [vcfpy.Call(sample=options.group, data={'GT': '0/1'})]
                                )
                                write_once(rec)
    writer.close()



def parse_arguments():
    p = configargparse.ArgParser()
    p.add('-i', '--input', required=True, help='Input folder with spladder results')
    p.add('-r', '--reference', required=True, help='Genome fasta file') 
    p.add('-g', '--group', required=True, help='Group name') # needed for vcf output
    p.add('-o', '--output', required=True, help='Output vcf file')
    p.add('-l', '--log', required=False, help='Log file')

    options = p.parse_args()

    return options


if __name__ == "__main__":
    main()

