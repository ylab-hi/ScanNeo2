import sys
import subprocess
import math
from pathlib import Path

"""
This script is used to prepare the input for the HLA typing (classII) using
HLA-HD. It requires paired-end reads. In this script the single-end reads in 
.fastq or .bam format are splitted into two files (if necessary)

Usage:
    python prep_input_mhcII_hlayping.py <input> <output_fwd> <output_rev>
"""

def calc_revcomp(seq):
    """
    Calculate the reverse complement of a given sequence
    """
    revcomp = ""
    for i in range(len(seq)):
        if seq[i] == "A":
            revcomp = "T" + revcomp
        elif seq[i] == "T":
            revcomp = "A" + revcomp
        elif seq[i] == "C":
            revcomp = "G" + revcomp
        elif seq[i] == "G":
            revcomp = "C" + revcomp
        else:
            revcomp = "N" + revcomp
    return revcomp


def split(infile, out_fwd, out_rev):
    print(f"Input file: {infile}")
    print(f"Reverse reads: {out_rev}")
    print(f"Forward reads: {out_fwd}")
    call = "samtools fastq"
    call += " -1 " + out_fwd + " -2 " + out_rev
    call += " -0 /dev/null" # discard unpaired reads
    call += " -s /dev/null" # discard singletons
    call += " " + infile
    print(call)
    subprocess.call(call, shell=True)
    print("Done.")

def main():
    inputfile = sys.argv[1]
    out_fwd = sys.argv[2]
    out_rev = sys.argv[3]

    fext = Path(inputfile).suffix # determine file extension
    print("Detected File extension: ", fext)
    if fext == ".fq":
        in_fh = open(inputfile, "r")
        first_line = in_fh.readline()

        # check if fq is interleaved (paired-end)
        if "/" in first_line:
            print(f"Interleaved fastq file detected. Splitting into two files...")
            split(inputfile, out_fwd, out_rev)

        else:
            print(f"Single-end fastq file detected. Splitting into two files...")
            print(f"Forward reads: {out_fwd}")
            print(f"Reverse reads: {out_rev}")
            in_fh = open(inputfile[0], "r")
            out_fwd_fh = open(out_fwd, "w")
            out_rev_fh = open(out_rev, "w")

            for i, val in enumerate(in_fh):
                line = val.strip()
                if i % 4 == 0: # identifier
                    out_fwd_fh.write(line+"/1\n")
                    out_rev_fh.write(line+"/2\n")
                if i % 4 == 1: # sequence
                    seqmid = math.ceil(len(line)/2)
                    out_fwd_fh.write(line[:seqmid+1] + "\n")
                    if seqmid > 10: 
                        seqmid -= 10
                    # small overlap for reverse reads
                    out_rev_fh.write(calc_revcomp(line[seqmid-10:]))
                if i % 4 == 2: # +
                    out_fwd_fh.write("+\n")
                    out_rev_fh.write("+\n")
                if i % 4 == 3: # quality
                    out_fwd_fh.write(line[:seqmid+1] + "\n")
                    out_rev_fh.write(line[seqmid-10:] + "\n")

            in_fh.close()
            out_fwd_fh.close()
            out_rev_fh.close()


    elif fext == ".bam": # bamfile is paired-end (at least in ScanNeo2)
        print(f"Paired-end bam file detected. Splitting into two files...")
        split(inputfile, out_fwd, out_rev)

main()
