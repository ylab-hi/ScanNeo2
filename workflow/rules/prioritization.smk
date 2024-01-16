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
    mhcI=get_mhcI,
    mhcII=get_mhcII,
    peptide="resources/refs/peptide.fasta",
    annotation="resources/refs/genome_tmp.gtf",
    counts="results/{sample}/quantification/allcounts.txt"
  output:
    "results/{sample}/prioritization/",
  message:
    "Predicting binding affinities on sample:{wildcards.sample}"
  conda:
    "../envs/prioritization.yml"
  log:
    "logs/prioritization/{sample}.log"
  threads: config["threads"]
  params:
    mhc_class = f"""config["priorization"]["class"]""",
    mhcI_len = f"""config["priorization"]["lengths"]["MHC-I"]""",
    mhcII_len = f"""config["priorization"]["lengths"]["MHC-II"]""",
  shell:
    """
      python workflow/scripts/priorization/compile.py \
          -i '{input.var}' -f {input.fus} \
          --output_dir results/{wildcards.sample}/priorization/ \
          -p {input.peptide} -a {input.annotation} \
          --confidence medium \
          --mhc_class {params.mhc_class} \
          --mhcI {input.mhcI} \
          --mhc_len {params.mhcI_len} \
          --mhcII {input.mhcII} \
          --mhcII_len {params.mhcII_len} \
          --counts {input.counts} \
          --threads {threads} \
          --output_dir {output}
    """
