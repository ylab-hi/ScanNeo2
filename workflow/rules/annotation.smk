rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    wrapper:
        "v1.31.1/bio/vep/plugins"

# annotation of long indels (as determined by transindel)
rule vep_longindels:
    input:
        calls="results/{sample}/rnaseq/indel/m2.indel.vcf",
        plugins="resources/vep/plugins",
        fasta="resources/refs/genome.fasta",
        fai="resources/refs/genome.fasta.fai", # fasta index
        gff="resources/refs/genome.gtf.gz",
        csi="resources/refs/genome.gtf.gz.csi",
    output:
        calls="results/{sample}/annotation/long.indel.vcf",
        stats="results/{sample}/annotation/long.indel.html"
    params:
        plugins=["NMD"],
        extra="--everything",  # optional: extra arguments
    log:
        "logs/vep/{sample}_annotate.log",
    threads: 4
    wrapper:
        "v1.31.1/bio/vep/annotate"


use rule vep_longindels as vep_shortindels with:
    input:
        calls="results/{sample}/rnaseq/indel/m2.indel.vcf",
        plugins="resources/vep/plugins"
    output:
        calls="results/{sample}/annotation/short.indel.vcf",
        stats="results/{sample}/annotation/short.indel.html"

use rule vep_longindels as vep_snvs with:
    input:
        calls="results/{sample}/rnaseq/indel/m2.snvs.vcf",
        plugins="resources/vep/plugins",
    output:
        calls="results/{sample}/rnaseq/annotation/somatic.snvs.vcf",
        stats="results/{sample}/rnaseq/annotation/somatic.snvs.html"











