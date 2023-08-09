rule download_prediction_tools:
  output:
    directory("workflow/scripts/mhc_i/")
  message:
    "Downloading MHC I prediction tools"
  log:
    "logs/download_mhc_i_tools.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o - https://downloads.iedb.org/tools/mhci/3.1.4/IEDB_MHC_I-3.1.4.tar.gz \
          | tar xz -C workflow/scripts/
    """

rule variants_to_peptides:
  input:
    var=get_variants,
    pep="resources/refs/peptides.fasta",
    alleles="results/{sample}/hla/alleles.tsv"

  output:
    "results/{sample}/neoantigens/results.tsv"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/variants_to_peptides.yml"
  params:
    length="-l {config["priorization"]["mhc_i"]["len"]}"
  shell:
    """
      python3 workflow/scripts/variants_to_peptide.py \
          -p {input.pep} \
          -v {input.var} \
          -o {output} \
          -a {input.alleles} \
          {params.length} > {log}
    """

rule all:
  input:
    "results/{sample}/neoantigens/results.tsv"
