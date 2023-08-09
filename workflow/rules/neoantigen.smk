#rule download_prediction_tools:
  #output:
    #""

  #shell:
    #"""
      #curl -L -o - https://downloads.iedb.org/tools/mhci/LATEST/IEDB_MHC_I-3.1.4.tar.gz | 
    #"""

rule variants_to_peptides:
  input:
    get_variants,
    pep="resources/refs/peptides.fasta",
    alleles="results/{sample}/hla/alleles.tsv"
  output:
    "results/{sample}/neoantigens/results.tsv"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/variants_to_peptides.yml"
  shell:
    """
      python3 workflow/scripts/variants_to_peptide.py \
          -p {input.pep} \
          -v {input[0]} \
          -o {output} \
          -a {input.alleles}
          -l {config['priorization']['len_class_i']}
    """

rule all:
  input:
    "results/{sample}/neoantigens/results.tsv"
