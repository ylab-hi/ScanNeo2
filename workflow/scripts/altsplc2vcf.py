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

    # search for input files (in spladder output)
    for root, dirs, files in os.walk(options.input):
        for file in files:
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
                            event_type == "mutex_exon"):

                            try:
                                chrom = column[0]
                                strand = column[1]
                                event_id = column[2]

                                # https://github.com/ratschlab/spladder/issues/168
                                annotated = column[3]

                                gene_name = column[4]

                                start = int(column[7])
                                end = int(column[8])
                            
                            # in some cases lines are incomplete (faulty spladder output)
                            except IndexError:
                                print("IndexError - skipping event")
                                continue

                            if (event_type == "alt_5prime" or
                                event_type == "alt_3prime" or
                                event_type == "intron_retention"):

                                sv_type = "INS"
                                sv_len = end - start + 1

                                ref = str(seq[chrom][start-1:start]).upper()
                                alt = (str(seq[chrom][start-1:start]) + str(seq[chrom][start:end+1])).upper()
                            
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
                                             ('SRC', event_type)]),

                                        FORMAT=['GT'],
                                        calls = [vcfpy.Call(sample=options.group, data={'GT': '0/1'})]
                                )
                                writer.write_record(rec)

                            elif (event_type == "exon_skip"):
                                sv_type = "DEL"
                                sv_len = -(end - start + 1)
                                
                                e2_start = column[7]
                                e2_end = column[8]

                                """ for some reason spladder reports exon skipping
                                events in which the start and end position of the 
                                skipped exon are the same... skip those events """
                                if e2_start == e2_end:
                                    continue

                                ref = (str(seq[chrom][start-1:start]) + str(seq[chrom][start:end+1])).upper()
                                alt = str(seq[chrom][start-1:start]).upper()
                                
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
                                             ('SRC', event_type)]),

                                        FORMAT=['GT'],
                                        calls = [vcfpy.Call(sample=options.group, data={'GT': '0/1'})]
                                )
                                writer.write_record(rec)

                            elif (event_type == "mutex_exons"):
                                sv_type = "DEL"

                                e2_start = column[7]
                                e2_end = column[8]

                                e3_start = column[9]
                                e3_end = column[10]

                                sv_len = []
                                sv_len.append(-(e2_end - e2_start + 1))
                                sv_len.append(-(e3_end - e3_start + 1))

                                ref = []
                                ref.append((str(seq[chrom][e2_start-1:e2_start]) + str(seq[chrom][e2_start:e2_end-1])).upper())
                                ref.append((str(seq[chrom][e3_start-1:e3_start]) + str(seq[chrom][e3_start:e3_end-1])).upper())

                                alt = []
                                alt.append(str(seq[chrom][e2_start-1:e2_start]).upper())
                                alt.append(str(seq[chrom][e3_start-1:e3_start]).upper())

                                ao = []
                                e1e2_conf = column[17]
                                e3e4_conf = column[20]
                                ao.append(e1e2_conf)
                                ao.append(e3e4_conf)

                                psi = column[21]

                                for i in range(0,2):

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
                                                 ('SRC', event_type)]),

                                            FORMAT=['GT'],
                                            calls = [vcfpy.Call(sample=options.group, data={'GT': '0/1'})]
                                    )
                                    writer.write_record(rec)
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

