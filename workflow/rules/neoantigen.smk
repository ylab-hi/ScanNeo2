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
    pep = "resources/refs/peptides.fasta"
  output:
    "results/{sample}/neoantigens/results.tsv"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/variants_to_peptides.yml"
  shell:
    """
      python3 workflow/scripts/variants_to_peptide.py \
          resources/refs/Homo_sapiens/GRCh38.pep.all.fa \
          {input} {output}
    """

rule all:
  output:
    "results/{sample}/neoantigens/results.tsv"
