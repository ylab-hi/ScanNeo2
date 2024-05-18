rule download_mhcI_ba_tools:
  output:
    directory("workflow/scripts/mhc_i/")
  message:
    "Downloading MHC I prediction tools"
  log:
    "logs/download_mhc-I_ba_tools.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o - https://downloads.iedb.org/tools/mhci/3.1.4/IEDB_MHC_I-3.1.4.tar.gz \
          | tar xz -C workflow/scripts/
    """

rule download_mhcII_ba_tools:
  output:
    directory("workflow/scripts/mhc_ii/")
  message:
    "Downloading MHC II prediction tools"
  log:
    "logs/download_mhc-II_ba_tools.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o - https://downloads.iedb.org/tools/mhcii/3.1/IEDB_MHC_II-3.1.tar.gz \
          | tar xz -C workflow/scripts/ > {log} 2>&1
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

rule prioritization:
  input:
    snv=get_prioritization_snvs,
    indels=get_prioritization_indels,
    long_indels=get_prioritization_long_indels,
    altsplicing=get_prioritization_altsplicing,
    exitrons=get_prioritization_exitrons,
    fusions=get_fusions,
    custom=get_prioritization_custom,
    mhcI=get_prioritization_mhcI,
    mhcII=get_prioritization_mhcII,
    refgenome="resources/refs/genome.fasta",
    peptide="resources/refs/peptide.fasta",
    annotation="resources/refs/genome_tmp.gtf",
    counts=get_prioritization_counts,
    mhcI_ba="workflow/scripts/mhc_i/",
    mhcI_im="workflow/scripts/immunogenicity/",
    mhcII_ba="workflow/scripts/mhc_ii/"
  output:
    directory("results/{sample}/prioritization/"),
  message:
    "Prioritize on sample:{wildcards.sample}"
  conda:
    "../envs/prioritization.yml"
  log:
    "logs/prioritization/{sample}.log"
  threads: config["threads"]
  params:
    mhc_class = f"""{config["prioritization"]["class"]}""",
    mhcI_len = f"""{config["prioritization"]["lengths"]["MHC-I"]}""",
    mhcII_len = f"""{config["prioritization"]["lengths"]["MHC-II"]}""",
  shell:
    """
      python workflow/scripts/prioritization/compile.py \
          --SNV "{input.snv}" \
          --indels "{input.indels}" \
          --long_indels "{input.long_indels}" \
          --exitrons "{input.exitrons}" \
          --altsplicing "{input.altsplicing}" \
          --fusions "{input.fusions}" \
          --custom "{input.custom}" \
          --proteome {input.peptide} \
          --anno {input.annotation} \
          --confidence medium \
          --mhc_class {params.mhc_class} \
          --mhcI "{input.mhcI}" \
          --mhcI_len "{params.mhcI_len}" \
          --mhcII "{input.mhcII}" \
          --mhcII_len "{params.mhcII_len}" \
          --counts "{input.counts}" \
          --threads {threads} \
          --output_dir {output} \
          --reference {input.refgenome} > {log} 2>&1
    """
