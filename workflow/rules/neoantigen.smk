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

rule priorization:
  input:
    var=get_variants,
    fus=get_fusions,
    mhcI="results/{sample}/hla/mhc-I.tsv",
    peptide="resources/refs/peptide.fasta",
    annotation="resources/refs/genome_tmp.gtf",
  output:
    effects = "results/{sample}/priorization/variant_effects.tsv",
    affinities = "results/{sample}/priorization/binding_affinities.tsv",
  message:
    "Predicting binding affinities on sample:{wildcards.sample}"
  conda:
    "../envs/priorization.yml"
  shell:
    """
      python workflow/scripts/priorization/compile.py \
          -i '{input.var}' -f {input.fus} \
          --output_dir results/{wildcards.sample}/priorization/ \
          -p {input.peptide} -a {input.annotation} \
          --confidence low \
          --mhcI {input.mhcI} \
          --mhc_len {config["priorization"]["lengths"]["MHC-I"]}
    """
