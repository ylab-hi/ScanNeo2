FROM condaforge/mambaforge:latest
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    git 

# # copy scanneo2 into container
# RUN git clone https://github.com/ylab-hi/ScanNeo2/
# # check out branch
# RUN cd ScanNeo2 && git checkout mhcIItyping

# install snakemake in conda environment
RUN conda create -n scanneo2 -c bioconda snakemake snakemake-wrapper-utils && \
    echo "source activate scanneo2" >> ~/.bashrc

# enable strict channels
RUN conda config --set channel_priority strict

# WORKDIR /ScanNeo2

SHELL ["/bin/bash", "-c"]

