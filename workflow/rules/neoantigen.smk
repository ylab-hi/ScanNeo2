rule candidate_filtering:
    input:
        "results/vep/output.tab"
    output:
        "results/candidates.tab"
    shell:
        """
        """

rule binding_prediction:
    input:
        cad="results/candidates.tab",
        tsv="results/hla/all_result.tsv",
    output:
        "results/candidates_binding.tab"
    shell:
        """
        """

rule immunogenecity:
    input:
        "results/candidates_binding.tab"
    output:
        "results/neoantigens_unfiltered.txt"
    shell:
        """
        """


rule final_filtering:
    input:
        "results/neoantigens_unfiltered.txt"
    output:
        "results/neoantigens.txt"
    shell:
        """
        """
        
