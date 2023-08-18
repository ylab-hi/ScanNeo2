FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="59e93cf8f17044dd4c67a8f4249f98cb59359406b22460eaee3da067704d29e3"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.26.0/bio/bwa/index/environment.yaml
#   prefix: /conda-envs/c2273db56372af81b557728cddb6ccd7
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bwa =0.7.17
RUN mkdir -p /conda-envs/c2273db56372af81b557728cddb6ccd7
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.26.0/bio/bwa/index/environment.yaml /conda-envs/c2273db56372af81b557728cddb6ccd7/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.26.0/bio/star/index/environment.yaml
#   prefix: /conda-envs/62ead550e24a99a2d873f1fc4ab2e9d3
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - star =2.7.10b
RUN mkdir -p /conda-envs/62ead550e24a99a2d873f1fc4ab2e9d3
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.26.0/bio/star/index/environment.yaml /conda-envs/62ead550e24a99a2d873f1fc4ab2e9d3/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.28.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/f49ecd58ca2f64505be52944e62793eb
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.17
#     - snakemake-wrapper-utils =0.5.3
RUN mkdir -p /conda-envs/f49ecd58ca2f64505be52944e62793eb
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.28.0/bio/samtools/faidx/environment.yaml /conda-envs/f49ecd58ca2f64505be52944e62793eb/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.31.1/bio/vep/annotate/environment.yaml
#   prefix: /conda-envs/fc332fe00ffc9fb8703741177c361800
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - ensembl-vep =109.3
#     - bcftools =1.16
#     - perl-encode-locale =1.05
#     - perl =5.32.1
RUN mkdir -p /conda-envs/fc332fe00ffc9fb8703741177c361800
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.31.1/bio/vep/annotate/environment.yaml /conda-envs/fc332fe00ffc9fb8703741177c361800/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.31.1/bio/vep/plugins/environment.yaml
#   prefix: /conda-envs/b823871fbdb958824831a5d19ef346f7
#   channels:
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - python =3.11.3
RUN mkdir -p /conda-envs/b823871fbdb958824831a5d19ef346f7
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.31.1/bio/vep/plugins/environment.yaml /conda-envs/b823871fbdb958824831a5d19ef346f7/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.2.1/bio/fastp/environment.yaml
#   prefix: /conda-envs/371fd38cc0bae1ae0eff1fdbcd1fb191
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastp =0.23.4
RUN mkdir -p /conda-envs/371fd38cc0bae1ae0eff1fdbcd1fb191
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.2.1/bio/fastp/environment.yaml /conda-envs/371fd38cc0bae1ae0eff1fdbcd1fb191/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.2.1/bio/fastqc/environment.yaml
#   prefix: /conda-envs/4c5e83029c4c995ec347593201305c58
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.5.3
RUN mkdir -p /conda-envs/4c5e83029c4c995ec347593201305c58
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.2.1/bio/fastqc/environment.yaml /conda-envs/4c5e83029c4c995ec347593201305c58/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.3.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/07acb90021639eb8e01afad917437175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.17
RUN mkdir -p /conda-envs/07acb90021639eb8e01afad917437175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.3.0/bio/samtools/index/environment.yaml /conda-envs/07acb90021639eb8e01afad917437175/environment.yaml

# Conda environment:
#   source: workflow/envs/basic.yml
#   prefix: /conda-envs/433868d00d0dc4be52f375e4825c470c
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.8
#     - bwa=0.7.17
#     - samtools=1.16.1
RUN mkdir -p /conda-envs/433868d00d0dc4be52f375e4825c470c
COPY workflow/envs/basic.yml /conda-envs/433868d00d0dc4be52f375e4825c470c/environment.yaml

# Conda environment:
#   source: workflow/envs/manipulate_vcf.yml
#   prefix: /conda-envs/6eaaf5e74784e6f8c26fdcd742e36151
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - pip=23.1.2
#     - pip:
#         - vcfpy==0.13.6
#         - intervaltree==3.1.0
RUN mkdir -p /conda-envs/6eaaf5e74784e6f8c26fdcd742e36151
COPY workflow/envs/manipulate_vcf.yml /conda-envs/6eaaf5e74784e6f8c26fdcd742e36151/environment.yaml

# Conda environment:
#   source: workflow/envs/optitype.yml
#   prefix: /conda-envs/fb6854ad4d0c4f83583dc02e8514c04a
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - optitype =1.3.5
#     - pysam
RUN mkdir -p /conda-envs/fb6854ad4d0c4f83583dc02e8514c04a
COPY workflow/envs/optitype.yml /conda-envs/fb6854ad4d0c4f83583dc02e8514c04a/environment.yaml

# Conda environment:
#   source: workflow/envs/priorization.yml
#   prefix: /conda-envs/68a7987b5f06060b54044245991e2e18
#   channels:
#   - bioconda
#   - conda-forge
#   dependencies:
#   - python>=3.6
#   - vcfpy
#   - pyfaidx
#   - configargparse
RUN mkdir -p /conda-envs/68a7987b5f06060b54044245991e2e18
COPY workflow/envs/priorization.yml /conda-envs/68a7987b5f06060b54044245991e2e18/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools.yml
#   prefix: /conda-envs/73c25277031fb2557ef8f6ce52038e12
#   channels:
#     - bioconda
#   dependencies:
#    - samtools
#    - bcftools
RUN mkdir -p /conda-envs/73c25277031fb2557ef8f6ce52038e12
COPY workflow/envs/samtools.yml /conda-envs/73c25277031fb2557ef8f6ce52038e12/environment.yaml

# Conda environment:
#   source: workflow/envs/transindel.yml
#   prefix: /conda-envs/1f213f31e028357cabefb1af5246e324
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - samtools
#     - pysam
#     - htseq
#     - pyfaidx
#     - pip
#     - pip:
#       - vcfpy==0.13.6
RUN mkdir -p /conda-envs/1f213f31e028357cabefb1af5246e324
COPY workflow/envs/transindel.yml /conda-envs/1f213f31e028357cabefb1af5246e324/environment.yaml

# Conda environment:
#   source: workflow/envs/yara.yml
#   prefix: /conda-envs/b92240a0a4d42d7bd8bbe57c8f338de6
#   channels:
#     - bioconda
#   dependencies:
#     - yara
#     - samtools
RUN mkdir -p /conda-envs/b92240a0a4d42d7bd8bbe57c8f338de6
COPY workflow/envs/yara.yml /conda-envs/b92240a0a4d42d7bd8bbe57c8f338de6/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/c2273db56372af81b557728cddb6ccd7 --file /conda-envs/c2273db56372af81b557728cddb6ccd7/environment.yaml && \
    mamba env create --prefix /conda-envs/62ead550e24a99a2d873f1fc4ab2e9d3 --file /conda-envs/62ead550e24a99a2d873f1fc4ab2e9d3/environment.yaml && \
    mamba env create --prefix /conda-envs/f49ecd58ca2f64505be52944e62793eb --file /conda-envs/f49ecd58ca2f64505be52944e62793eb/environment.yaml && \
    mamba env create --prefix /conda-envs/fc332fe00ffc9fb8703741177c361800 --file /conda-envs/fc332fe00ffc9fb8703741177c361800/environment.yaml && \
    mamba env create --prefix /conda-envs/b823871fbdb958824831a5d19ef346f7 --file /conda-envs/b823871fbdb958824831a5d19ef346f7/environment.yaml && \
    mamba env create --prefix /conda-envs/371fd38cc0bae1ae0eff1fdbcd1fb191 --file /conda-envs/371fd38cc0bae1ae0eff1fdbcd1fb191/environment.yaml && \
    mamba env create --prefix /conda-envs/4c5e83029c4c995ec347593201305c58 --file /conda-envs/4c5e83029c4c995ec347593201305c58/environment.yaml && \
    mamba env create --prefix /conda-envs/07acb90021639eb8e01afad917437175 --file /conda-envs/07acb90021639eb8e01afad917437175/environment.yaml && \
    mamba env create --prefix /conda-envs/433868d00d0dc4be52f375e4825c470c --file /conda-envs/433868d00d0dc4be52f375e4825c470c/environment.yaml && \
    mamba env create --prefix /conda-envs/6eaaf5e74784e6f8c26fdcd742e36151 --file /conda-envs/6eaaf5e74784e6f8c26fdcd742e36151/environment.yaml && \
    mamba env create --prefix /conda-envs/fb6854ad4d0c4f83583dc02e8514c04a --file /conda-envs/fb6854ad4d0c4f83583dc02e8514c04a/environment.yaml && \
    mamba env create --prefix /conda-envs/68a7987b5f06060b54044245991e2e18 --file /conda-envs/68a7987b5f06060b54044245991e2e18/environment.yaml && \
    mamba env create --prefix /conda-envs/73c25277031fb2557ef8f6ce52038e12 --file /conda-envs/73c25277031fb2557ef8f6ce52038e12/environment.yaml && \
    mamba env create --prefix /conda-envs/1f213f31e028357cabefb1af5246e324 --file /conda-envs/1f213f31e028357cabefb1af5246e324/environment.yaml && \
    mamba env create --prefix /conda-envs/b92240a0a4d42d7bd8bbe57c8f338de6 --file /conda-envs/b92240a0a4d42d7bd8bbe57c8f338de6/environment.yaml && \
    mamba clean --all -y
