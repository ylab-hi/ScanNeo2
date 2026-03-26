# preliminary rules
rule create_tmp_folder:
  output:
    directory("tmp/")
  message:
    "Creating temporary folder"
  conda:
    "../envs/basic.yml"
  log:
    "logs/ref/create_tmp_folder.log"
  shell:
    """
      mkdir -p {output} > {log} 2>&1
    """

