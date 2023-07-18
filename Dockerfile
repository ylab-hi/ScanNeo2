yes
FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="cc0971e5f0db6491b8c7ea4fcb1b8c181f666c0a7975d5a224d603fec0e6440a"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.1.1/bio/trimmomatic/pe/environment.yaml
#   prefix: /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trimmomatic ==0.36
#     - pigz ==2.3.4
#     - snakemake-wrapper-utils ==0.1.3
RUN mkdir -p /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.1.1/bio/trimmomatic/pe/environment.yaml /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe/environment.yaml

# Conda environment:
#   source: workflow/envs/basic.yml
#   prefix: /conda-envs/f05d91be109b4d9b1e6f8f8371699cb6
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.8
#     - bwa=0.7.17
#     - samtools
RUN mkdir -p /conda-envs/f05d91be109b4d9b1e6f8f8371699cb6
COPY workflow/envs/basic.yml /conda-envs/f05d91be109b4d9b1e6f8f8371699cb6/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe --file /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe/environment.yaml && \
    mamba env create --prefix /conda-envs/f05d91be109b4d9b1e6f8f8371699cb6 --file /conda-envs/f05d91be109b4d9b1e6f8f8371699cb6/environment.yaml && \
    mamba clean --all -y
