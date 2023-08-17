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

rule compile_peptides_from_variants:
  input:
    var=get_variants,
    alleles="results/{sample}/hla/alleles.tsv"
  output:
    "results/{sample}/neoantigens/peptides.tsv"
  message:
    "Compile peptides from variants on sample:{wildcards.sample}"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/priorization.yml"
  params:
  shell:
    """
      python3 workflow/scripts/compile_peptides_from_variants.py \
          -v '{input.var}' -o {output} > {log}
    """

rule priorization:
  input:
    peptides="results/{sample}/neoantigens/peptides.tsv",
    alleles="results/{sample}/hla/alleles.tsv"
  output:
    "results/{sample}/neoantigens/final.tsv"
  message:
    "Predicting affinities on sample:{wildcards.sample}"
  log:
    "logs/mhci/{sample}_affinities.log"
  conda:
    "../envs/priorization.yml""
  shell:
    """
      python workflow/scripts/predict_affinities.py \
          {input.peptides} {input.alleles} {output} > {log}
      python workflow/scripts/predict_immunogenicity.py \
          {output} >> {log}
    """







