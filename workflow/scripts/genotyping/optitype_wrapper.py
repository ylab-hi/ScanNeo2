import sys
import os
import subprocess
import tempfile
from pathlib import Path

"""
    Usage:
        python3 optitype_wrapper.py <nartype> <prefix> <outpath> <bam1> [bam2] [...]
"""

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
            # GATK SplitSamByNumberOfReads emits name-grouped (SO:unsorted)
            # chunks so read mates stay paired across the split; samtools index
            # requires coordinate order, so sort each chunk before indexing and
            # hand the sorted copies to OptiType.
            sorted_bams = []
            for filename in inbams:
                sorted_bam = os.path.join(tmp, os.path.basename(filename))
                subprocess.run(
                    ["samtools", "sort", "-o", sorted_bam, filename], check=True)
                subprocess.run(["samtools", "index", sorted_bam], check=True)
                sorted_bams.append(sorted_bam)

            result = subprocess.run(
                ["samtools", "view", "-c", sorted_bams[0]],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                universal_newlines=True, check=True)

            if int(result.stdout.strip()) < 10:
                print("Input BAM file: " + sorted_bams[0] + " is empty")
                Path(outpath + prefix + "_coverage_plot.pdf").touch()
                Path(outpath + prefix + "_result.tsv").touch()
            else:
                optitype_cmd = [
                    "OptiTypePipeline.py",
                    "--input", *sorted_bams,
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
