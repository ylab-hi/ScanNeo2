# The BLAST protein database files that `blastp` requires. Measured against the
# pinned blast: removing .pdb makes blastp abort with "File ....pdb not found",
# while the other v5 metadata files it writes (.pot, .ptf, .pto, .pjs) serve
# blastdbcmd/taxonomy lookups and are individually removable with identical
# results, so they are deliberately not tracked. Declared once and shared by the
# rule that builds the database and the rule that reads it, so the two sets
# cannot drift apart.
PROTEOME_BLASTDB = multiext(
    "resources/refs/proteome_blastdb", ".phr", ".pin", ".psq", ".pdb"
)


rule download_mhcI_ba_tools:
    output:
        directory("workflow/scripts/mhc_i/"),
    log:
        "logs/download/mhcI_ba_tools.log",
    conda:
        "../envs/basic.yml"
    message:
        "Downloading MHC I prediction tools"
    shell:
        """
        (curl --fail -L -o - https://downloads.iedb.org/tools/mhci/3.1.4/IEDB_MHC_I-3.1.4.tar.gz \
            | tar xz -C workflow/scripts/) >{log} 2>&1
        """


rule download_mhcII_ba_tools:
    output:
        directory("workflow/scripts/mhc_ii/"),
    log:
        "logs/download/mhcII_ba_tools.log",
    conda:
        "../envs/basic.yml"
    message:
        "Downloading MHC II prediction tools"
    shell:
        """
        (
            curl --fail -L -o - https://downloads.iedb.org/tools/mhcii/3.1.12/IEDB_MHC_II-3.1.12.tar.gz \
                | tar xz -C workflow/scripts/
            cp misc/mhc_ii/methods/netmhciipan-3.2-executable/netmhciipan_3_2_executable/netmhciipan_python_interface.py \
                workflow/scripts/mhc_ii/methods/netmhciipan-3.2-executable/netmhciipan_3_2_executable/
            cp misc/mhc_ii/methods/netmhciipan-4.1-executable/netmhciipan_4_1_executable/netmhciipan_python_interface.py \
                workflow/scripts/mhc_ii/methods/netmhciipan-4.1-executable/netmhciipan_4_1_executable/
            cp misc/mhc_ii/methods/netmhciipan-4.2-executable/netmhciipan_4_2_executable/netmhciipan_python_interface.py \
                workflow/scripts/mhc_ii/methods/netmhciipan-4.2-executable/netmhciipan_4_2_executable/
            cp misc/mhc_ii/methods/netmhciipan-4.3-executable/netmhciipan_4_3_executable/netmhciipan_python_interface.py \
                workflow/scripts/mhc_ii/methods/netmhciipan-4.3-executable/netmhciipan_4_3_executable/
        ) >{log} 2>&1
        """


rule download_prediction_binding_affinity_tools:
    output:
        directory("workflow/scripts/immunogenicity/"),
    log:
        "logs/download/immunogenicity_tools.log",
    conda:
        "../envs/basic.yml"
    message:
        "Downloading immunogenicity prediction tools"
    shell:
        """
        (curl --fail -L -o - https://downloads.iedb.org/tools/immunogenicity/3.0/IEDB_Immunogenicity-3.0.tar.gz \
            | tar xz -C workflow/scripts/) >{log} 2>&1
        """


rule make_proteome_blastdb:
    input:
        peptide="resources/refs/peptide.fasta",
    output:
        PROTEOME_BLASTDB,
    log:
        "logs/ref/make_proteome_blastdb.log",
    conda:
        "../envs/prioritization.yml"
    params:
        # derived from the output rather than hardcoded, so the prefix still
        # resolves when the outputs are staged (no shared filesystem)
        prefix=lambda w, output: os.path.splitext(output[0])[0],
    message:
        "Building BLAST database for the reference proteome"
    shell:
        """
        makeblastdb -in {input.peptide} -dbtype prot -out {params.prefix} \
            >{log} 2>&1
        """


# the sources whose prediction is large enough to use a full node; the rest
# (a few to a few thousand records) finish in about a minute on 8
PRIORITIZATION_LARGE_SOURCES = {"somatic.snvs", "somatic.short.indels", "altsplicing"}
PRIORITIZATION_CLASSES = {"I": ["I"], "II": ["II"], "BOTH": ["I", "II"]}[
    config["prioritization"]["class"]
]


# One job per (sample, source), so a sample's sources are predicted
# concurrently and its wall-clock is its slowest source rather than their sum.
rule prioritize_source:
    input:
        variants=get_prioritization_source,
        mhcI=get_prioritization_mhcI,
        mhcII=get_prioritization_mhcII,
        refgenome="resources/refs/genome.fasta",
        # built by its own rule, so concurrent source jobs never race to create
        # it on first open
        refgenome_idx="resources/refs/genome.fasta.fai",
        peptide="resources/refs/peptide.fasta",
        annotation="resources/refs/genome_tmp.gtf",
        counts=get_prioritization_counts,
        mhcI_ba=get_mhcI_ba_tools,
        mhcII_ba=get_mhcII_ba_tools,
        mhcI_im=get_mhcI_immunogenicity_tools,
        proteome_db=PROTEOME_BLASTDB,
    output:
        directory("results/{sample}/prioritization/{source}/"),
    log:
        "logs/{sample}/prioritization/{source}.log",
    wildcard_constraints:
        source="|".join(re.escape(s) for s in PRIORITIZATION_SOURCES),
    conda:
        "../envs/prioritization.yml"
    # The binding-affinity pool is one set of (allele x epitope length x wt/mt x
    # FASTA batch) units, so it uses as many threads as it is given, unlike the
    # rest of the workflow; 48 fills a 52-core node. Snakemake clamps this to
    # --cores, so a smaller local run is unaffected.
    threads: lambda wildcards: 48 if wildcards.source in PRIORITIZATION_LARGE_SOURCES else 8
    params:
        flag=lambda wildcards: PRIORITIZATION_SOURCES[wildcards.source][1],
        mhc_class=f"""{config["prioritization"]["class"]}""",
        mhcI_len=f"""{config["prioritization"]["lengths"]["MHC-I"]}""",
        mhcII_len=f"""{config["prioritization"]["lengths"]["MHC-II"]}""",
    message:
        "Prioritize {wildcards.source} on sample:{wildcards.sample}"
    # The combined tables are removed first: they are not this job's output, so
    # a failed run would otherwise leave an earlier sample-wide table in place,
    # and the report takes a present table as a finished sample.
    shell:
        """
        rm -f results/{wildcards.sample}/prioritization/mhc-*_neoepitopes_all.txt
        python workflow/scripts/prioritization/compile.py \
            {params.flag} "{input.variants}" \
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
            --reference {input.refgenome} >{log} 2>&1
        """


rule combine_neoepitopes:
    input:
        get_prioritization_source_dirs,
    output:
        "results/{sample}/prioritization/mhc-{cls}_neoepitopes_all.txt",
    log:
        "logs/{sample}/prioritization/combine_neoepitopes_mhc-{cls}.log",
    wildcard_constraints:
        cls="I|II",
    # a concatenation, so it runs on the controller instead of queueing a job
    localrule: True
    conda:
        "../envs/basic.yml"
    params:
        # each source directory holds <source>_mhc-<cls>_neoepitopes.txt
        tables=lambda wildcards, input: [
            os.path.join(
                d,
                f"{os.path.basename(os.path.normpath(d))}_mhc-{wildcards.cls}_neoepitopes.txt",
            )
            for d in input
        ],
    message:
        "Combining the per-source MHC-{wildcards.cls} neoepitope tables on sample:{wildcards.sample}"
    # header from the first table, the others without theirs, in source order;
    # /dev/null keeps awk off stdin when a sample has no source
    shell:
        """
        awk 'FNR > 1 || NR == 1' {params.tables} /dev/null >{output} 2>{log}
        """


# One row per (sample, source) that ran, so a source that came out empty is
# visible without opening its tables; the status report lists those rows.
rule summarize:
    input:
        dirs=[e["dir"] for e in summary_entries()],
        sources=[p for e in summary_entries() for p in e["inputs"]],
    output:
        "results/summary.tsv",
    log:
        "logs/summarize.log",
    localrule: True
    conda:
        "../envs/basic.yml"
    params:
        # one SAMPLE|SOURCE|DIR|INPUT[,INPUT...] token per source that runs
        entries=[
            "|".join([e["sample"], e["source"], e["dir"], ",".join(e["inputs"])])
            for e in summary_entries()
        ],
        classes=PRIORITIZATION_CLASSES,
    message:
        "Summarizing the per-source prioritization results of all samples"
    shell:
        """
        python workflow/scripts/summarize.py \\
        --output {output} \\
        --classes {params.classes} \\
        --entries {params.entries:q} >{log} 2>&1
        """
