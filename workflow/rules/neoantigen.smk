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
      curl -L -o - https://downloads.iedb.org/tools/immunogenicity/3.0/IEDB_Immunogenicity-3.0.tar.gz \
          | tar xz -C workflow/scripts/
    """

rule variants_to_peptides:
  input:
    var=get_variants,
    alleles="results/{sample}/hla/alleles.tsv"
  output:
    "results/{sample}/neoantigens/results.tsv"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/variants_to_peptides.yml"
  params:
    length="{config['priorization']['mhc_i']['len']}"
  shell:
    """
      python3 workflow/scripts/variants_to_peptides.py \
          -v '{input.var}' -o {output} \
          -a {input.alleles} \
          -o {output} \
          -p 8-11 > {log}
    """
