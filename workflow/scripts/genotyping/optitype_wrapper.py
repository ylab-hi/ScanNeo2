import sys
import os
import subprocess

"""
    Usage:
        python3 optitype_wrapper.py <input> <nartype> <prefix> <outpath>
"""

def main():
    inbams = sys.argv[1]
    nartype = "dna" if sys.argv[2] == "DNA" else "rna"
    prefix = sys.argv[3]
    outpath = sys.argv[4]

    for filename in inbams.split(' '):
        if not os.path.exists(filename + '.bai'):
            print("Indexing BAM file: " + filename)
            samtools_idx = subprocess.Popen("samtools index " + 
                                            filename, shell=True)

    # check if input file is not empty (using samtools)
    try:
        samtools = subprocess.Popen("samtools view -c " + inbams.split(' ')[0] + " 2>/dev/null", 
                                    stdout=subprocess.PIPE, shell=True)

        if int(samtools.communicate()[0]) < 10:
            print("Input BAM file: " + inbams + " is empty")
            # create an empty output file
            pdf = subprocess.Popen("touch " + outpath + prefix + "_coverage_plot.pdf", shell=True)
            tsv = subprocess.Popen("touch " + outpath + prefix + "_result.tsv", shell=True)
        else:
            # call optitype 
            optitype_call = "OptiTypePipeline.py --input " + inbams + " --outdir " + outpath + " --prefix " + prefix + " --" + nartype + " -v"
            print(optitype_call)
            os.system(optitype_call)
#            optitype = subprocess.Popen(optitype_call, shell=True)


    except subprocess.CalledProcessError as e:
        print("samtools failed")

main()
