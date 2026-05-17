import os
import sys
import tempfile
import subprocess

"""
    This scripts is a wrapper for featureCounts. It checks if the input is in paired-end or 
    single-end mode and calls the appropriate function.

    Usage:
        featureCounts_wrapper.py <input_bam> <output> <annotation> <mapq> <threads>
"""
def main():
    # determine 
    inputbam = sys.argv[1]
    outputfile = sys.argv[2]
    annofile = sys.argv[3]
    mapq = sys.argv[4]
    threads = sys.argv[5]

    mode = get_mode(inputbam)
    featureCounts(mode, inputbam, outputfile, annofile, mapq, threads)


def get_mode(inputbam):
    """ This function checks if the input is in paired-end or single-end mode. """
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run(
        ["samtools", "view", "-f", "1", inputbam, "-o", tmpfile.name],
        check=True)

    if os.stat(tmpfile.name).st_size == 0:
        return "SE"
    else:
        return "PE"

def featureCounts(mode, inputfile, outputfile, annofile, mapq, threads):
    """ This function calls featureCounts"""
    cmd = ["featureCounts"]
    if mode == "PE":
        cmd.append("-p")
    cmd.extend([
        "-F", "GTF",
        "-t", "gene",
        "-g", "gene_id",
        "--fracOverlap", "0.2",
        "-Q", str(mapq),
        "-T", str(threads),
        "-a", annofile,
        "-o", outputfile,
        inputfile,
    ])

    # run featureCounts
    subprocess.run(cmd, check=True)

main()
