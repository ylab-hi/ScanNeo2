rule augment_custom_variants:
  input:
    get_custom_variants
  output:
    "results/{sample}/custom/custom.variants_augmented.vcf"
  message:
    "Augmenting custom variants on sample:{wildcards.sample}"
  log:
    "logs/{sample}/custom/augment_custom_variants.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/add_infos_to_vcf.py \
          {input} \
          custom \
          custom \
          {output} > {log} 2>&1
    """

rule sort_custom_variants:
  input:
    "results/{sample}/custom/custom.variants_augmented.vcf"
  output:
    "results/{sample}/variants/custom.vcf.gz"
  message:
    "Sorting and compressing exitrons on sample:{wildcards.sample}"
  log:
    "logs/{sample}/custom/sort_custom_variants.log"
  conda:
    "../envs/bcftools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} > {log} 2>&1
    """
