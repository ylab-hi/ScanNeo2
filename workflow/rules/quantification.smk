rule feature_counts:
    input:
        # list of sam or bam files
        samples=get_aligned_reads,
        annotation="resources/refs/genome.gtf"
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/{sample}/{seqtype}/featurecounts/{group}",
            ".txt",
            ".txt.summary",
            ".txt.jcounts",
        ),
    threads: 2
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra=f"""--fracOverlap 0.2 -J -f -p -t exon -g transcript_id \
          -Q {config['mapq']}""",
    log:
        "logs/{sample}/featurecounts/{seqtype}_{group}.log",
    wrapper:
        "v2.7.0/bio/subread/featurecounts"

