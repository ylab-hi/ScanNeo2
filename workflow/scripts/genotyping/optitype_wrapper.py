import sys
import os
import subprocess
from pathlib import Path

"""
    Usage:
        python3 optitype_wrapper.py <input> <nartype> <prefix> <outpath>
"""

def main():
    inbams_arg = sys.argv[1]
    nartype = "dna" if sys.argv[2] == "DNA" else "rna"
    prefix = sys.argv[3]
    outpath = sys.argv[4]

    inbams = inbams_arg.split()

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
            print("Input BAM file: " + inbams_arg + " is empty")
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
