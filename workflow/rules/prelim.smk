# preliminary rules
rule create_tmp_folder:
  output:
    directory("tmp/")
  message:
    "Creating temporary folder"
  conda:
    "../envs/basic.yml"
  log:
    "logs/prelim.log"
  shell:
    """
      mkdir -p {output}
    """

