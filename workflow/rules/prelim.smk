# preliminary rules
rule create_tmp_folder:
  output:
    directory("tmp/")
  message:
    "Creating temporary folder"
  conda:
    "../envs/basic.yml"
  shell:
    """
      mkdir -p {output}
    """

