import os
import sys
import pysam

def split(input_bamfile, outdir):
    # Open the input BAM file
    bam = pysam.AlignmentFile(input_bamfile, "rb")
    chromosomes = list(bam.references)

    os.makedirs(outdir, exist_ok=True)

    # create a dictionary of output files
    for chrom in chromosomes:
        # define output path
        outbam = os.path.join(outdir, f"{chrom}.bam")

        with pysam.AlignmentFile(outbam, "wb", template=bam) as outbam_fh:
            for read in bam.fetch(chrom):
                outbam_fh.write(read)

        # Index each per-chromosome BAM in place so the .bai lives inside the
        # checkpoint's directory() output. A standalone index rule emits a .bai
        # nested under that directory, which snakemake rejects as a
        # ChildIOException when the consuming rule runs as a remote SLURM jobstep.
        pysam.index(outbam)

    bam.close()

    print(f"Splitting {input_bamfile} into {outdir} is done.")


def main():
    bamfile = sys.argv[1]
    outdir = sys.argv[2]
    split(bamfile, outdir)


if __name__ == '__main__':
    main()
