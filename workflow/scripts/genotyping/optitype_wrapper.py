import sys
import os
import subprocess
import tempfile
from pathlib import Path

"""
    Usage:
        python3 optitype_wrapper.py <nartype> <prefix> <outpath> <bam1> [bam2] [...]
"""

# HLA typing saturates well below this depth, but OptiType's hit matrix grows
# with read count -- both memory AND runtime. One TESLA RNA sample had ~455k HLA
# reads: uncapped it OOM'd at 64 GB, and even capped at 100k OptiType spent ~1.5h
# building the matrix and timed out. 50k keeps thousands-fold coverage per HLA
# locus (no typing accuracy lost) while bounding memory and runtime.
READ_CAP = 50000


def count_reads(bam):
    out = subprocess.run(
        ["samtools", "view", "-c", bam],
        stdout=subprocess.PIPE, universal_newlines=True, check=True)
    return int(out.stdout.strip())


def main():
    if len(sys.argv) < 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <nartype> <prefix> <outpath> "
            f"<bam1> [bam2] [...]")

    nartype = "dna" if sys.argv[1] == "DNA" else "rna"
    prefix = sys.argv[2]
    outpath = sys.argv[3]
    inbams = sys.argv[4:]

    try:
        with tempfile.TemporaryDirectory() as tmp:
            # Decide one subsample fraction from the first input and apply the
            # SAME seed+fraction to every input, so samtools keeps identical
            # read names across R1/R2 and mates stay paired. Then coordinate-sort
            # and index each (razers3/OptiType read BAM, and samtools index needs
            # coordinate order; the split chunks arrive name-sorted).
            n = count_reads(inbams[0])
            subsample = None
            if n > READ_CAP:
                subsample = "{:.6f}".format(42 + READ_CAP / n)  # 42 = seed
                print(f"Subsampling {n} -> ~{READ_CAP} reads/file (samtools -s {subsample})")

            prepped = []
            for filename in inbams:
                dst = os.path.join(tmp, os.path.basename(filename))
                if subsample:
                    sub = dst + ".sub.bam"
                    subprocess.run(
                        ["samtools", "view", "-s", subsample, "-b", filename, "-o", sub],
                        check=True)
                    subprocess.run(["samtools", "sort", "-o", dst, sub], check=True)
                else:
                    subprocess.run(["samtools", "sort", "-o", dst, filename], check=True)
                subprocess.run(["samtools", "index", dst], check=True)
                prepped.append(dst)

            if count_reads(prepped[0]) < 10:
                print("Input BAM file: " + prepped[0] + " is empty")
                Path(outpath + prefix + "_coverage_plot.pdf").touch()
                Path(outpath + prefix + "_result.tsv").touch()
            else:
                optitype_cmd = [
                    "OptiTypePipeline.py",
                    "--input", *prepped,
                    "--outdir", outpath,
                    "--prefix", prefix,
                    "--" + nartype,
                    "-v",
                ]
                print(" ".join(optitype_cmd))
                subprocess.run(optitype_cmd, check=True)

    except subprocess.CalledProcessError as e:
        print("samtools failed: " + str(e))
        raise

main()
