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
    "results/{sample}/priorization/peptides.tsv"
  message:
    "Compile peptides from variants on sample:{wildcards.sample}"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/priorization.yml"
  params:
  shell:
    """
      python workflow/scripts/compile_peptides_from_variants.py \
          -i '{input.var}' -o {output} > {log}
    """

rule priorization:
  input:
    peptides="results/{sample}/priorization/peptides.tsv",
    alleles="results/{sample}/hla/alleles.tsv"
  output:
    "results/{sample}/priorization/neoantigens.tsv"
  message:
    "Predicting affinities on sample:{wildcards.sample}"
  log:
    "logs/mhci/{sample}_affinities.log"
  threads: config["threads"]
  params:
    epitope_length=config["priorization"]["mhc_i"]["length"]
  conda:
    "../envs/priorization.yml"
  shell:
    """
      python workflow/scripts/predict_affinities.py \
          -i {input.peptides} -a {input.alleles} \
          -e {params.epitope_length} -t {threads} -o {output} 
      python workflow/scripts/predict_immunogenicity.py \
          {output}
    """
