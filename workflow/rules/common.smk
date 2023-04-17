def get_fqs(wildcards):
    return config["samples"][wildcards.sample]

#rule merge_vcf:
#    input:
#        "results/exitron/final.vcf"
#        "results/splicing/final.vcf"
#        "results/genefusion/final.vcf"
#    output:
#        "results/variants.vcf"
#
#    shell:
#        "workflow/scripts/merge_vcf.py {config[refgen]} {input} {output} 2> logs/error.err"
#        
        
