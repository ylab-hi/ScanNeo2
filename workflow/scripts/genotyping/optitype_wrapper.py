import sys
import os
import subprocess
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

    for filename in inbams:
        if not os.path.exists(filename + '.bai'):
            print("Indexing BAM file: " + filename)
            subprocess.run(["samtools", "index", filename], check=True)

    try:
        result = subprocess.run(
            ["samtools", "view", "-c", inbams[0]],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            universal_newlines=True, check=True)

        if int(result.stdout.strip()) < 10:
            print("Input BAM file: " + inbams[0] + " is empty")
            Path(outpath + prefix + "_coverage_plot.pdf").touch()
            Path(outpath + prefix + "_result.tsv").touch()
        else:
            optitype_cmd = [
                "OptiTypePipeline.py",
                "--input", *inbams,
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
