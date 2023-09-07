#!/usr/bin/env python
#-*- coding: utf-8 -*-
#===============================================================================
# convert exitron results table file to VCF file
#===============================================================================
import os
import argparse
from pyfaidx import Fasta
import configparser


def config_getter():
    this_dir = os.path.dirname(os.path.realpath(__file__))
    config_default = os.path.join(os.path.dirname(this_dir), '../resources/scanexitron_config.ini')
    config = configparser.ConfigParser(os.environ)
    config.read(config_default)
    hg38_ref   = config.get("fasta", "hg38")
    hg19_ref   = config.get("fasta", "hg19")
    return {'hg38_ref': hg38_ref, 'hg19_ref': hg19_ref}


def exitron2vcf(in_file, out_vcf, ref='hg38', config=config_getter()):
    if ref == 'hg19':
        fasta = config['hg19_ref']
    elif ref == 'hg38':
        fasta = config['hg38_ref']

    genome_seq = Fasta(fasta, sequence_always_upper=True, as_raw=True)
    sample = os.path.splitext(os.path.basename(in_file))[0]
    vcf_header ='''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">
##INFO=<ID=DP,Number=2,Type=Integer,Description="Total depth of junction">
##INFO=<ID=SpliecedSite,Number=2,Type=String,Description="Splieced site">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Junction strand">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=AO,Number=1,Type=Integer,Description="Reads support of the exitron">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=PSO,Number=1,Type=Float,Description="Percent spliced-out">
##INFO=<ID=PSI,Number=1,Type=Float,Description="Percent spliced-in">
##INFO=<ID=GeneName,Number=1,Type=String,Description="Gene name">
##INFO=<ID=GeneID,Number=1,Type=String,Description="Gene ID">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}
'''.format(sample)
    outfile = open(out_vcf, 'w')
    outfile.write(vcf_header)
    with open(in_file) as f:
        f.readline()
        for line in f:
            l = line.rstrip().split('\t')
            chrom = l[0]
            pos = int(l[1])
            end = int(l[2])
            _alt = genome_seq[chrom][pos-1:pos]
            _ref = genome_seq[chrom][pos-1:end-1]
            idx = l[3]
            ao = int(l[4])
            strand = l[5]
            gene_name = l[6]
            length = int(l[7])
            #end = pos + length
            spliced_site = l[8]
            gene_id = l[9].split('.')[0]
            pso = float(l[10])
            psi = float(l[11])
            dp = int(l[12])
            outfile.write('{}\t{}\t{}\t{}\t{}\t.\t.\tSVTYPE=DEL;END={};AO={};DP={};STRAND={};SpliecedSite={};GeneName={};GeneID={};SVLEN=-{};PSO={};PSI={}\tGT\t0/1\n'.format(chrom,pos,idx,_ref,_alt,end-1,ao,dp,strand,spliced_site,gene_name,gene_id,length,pso,psi))
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', action='store', dest='input', help="Input Exitron results table", required=True)
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="output type (default: %(default)s)", choices=['hg19', 'hg38'], default='hg38')
    parser.add_argument('-o', '--output', action='store', dest='output', help="(default: %(default)s)",default='output.vcf')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    exitron2vcf(in_file=args.input, out_vcf=args.output, ref=args.ref)
