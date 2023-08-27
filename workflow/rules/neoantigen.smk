rule download_prediction_affinity_tools:
  output:
    directory("workflow/scripts/mhc_i/")
  message:
    "Downloading MHC I prediction tools"
  log:
    "logs/download_prediction_affinity_tools.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o - https://downloads.iedb.org/tools/mhci/3.1.4/IEDB_MHC_I-3.1.4.tar.gz \
          | tar xz -C workflow/scripts/
    """

rule download_prediction_binding_affinity_tools:
  output:
    directory("workflow/scripts/immunogenicity/")
  message:
    "Downloading immunogenicity prediction tools"
  log:
    "logs/download_immunogenicity_tools.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
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
    pred_aff="workflow/scripts/mhc_i/",
    pred_imm="workflow/scripts/immunogenicity/",
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
    mhc_i=config["priorization"]["lengths"]["MHC-I"]
  conda:
    "../envs/priorization.yml"
  shell:
    """
      python workflow/scripts/predict_affinities.py \
          -i {input.peptides} -a {input.alleles} \
          -e {params.mhc_i} -t {threads} -o {output} 
      python workflow/scripts/predict_immunogenicity.py \
          {output}
    """
