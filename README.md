<div align="left">
    <h1>ScanNeo2</h1>
    <img src="https://img.shields.io/github/v/release/ylab-hi/ScanNeo2">
    <img src="https://github.com/ylab-hi/ScanNeo2/actions/workflows/linting.yml/badge.svg" alt="Workflow status badge">
    <img src="https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg">
    <img alt="Citation Badge" src="https://api.juleskreuer.eu/citation-badge.php?doi=10.1093/bioinformatics/btad659">
    <img src="https://img.shields.io/github/downloads/ylab-hi/ScanNeo2/total.svg">
    <img src="https://img.shields.io/github/contributors/ylab-hi/ScanNeo2">
    <img src="https://img.shields.io/github/last-commit/ylab-hi/ScanNeo2">
    <img src="https://img.shields.io/github/commits-since/ylab-hi/ScanNeo2/latest">
    <img src="https://img.shields.io/github/stars/ylab-hi/ScanNeo2?style=social">
    <img src="https://img.shields.io/github/forks/ylab-hi/ScanNeo2?style=social">
</div>

## What is ScanNeo2
`Scanneo2` is a snakemake workflow for the prediction of neoantigens from multiple sources. In its current state, 
this includes canonical-splicing, exitron-splicing, gene fusion, indels and snvs.

## Getting Started

In principle, Scanneo2 aims to resolve its dependencies automatically and requires only snakemake and snakedeploy.

### Prerequisites

Before installing and using ScanNeo2, make sure you have the following software installed:

1. [Mamba](https://github.com/conda-forge/miniforge#mambaforge): An open-source package manager. 
Mamba should be installed independently by the user.

### Installation

To get started with ScanNeo2, follow the steps below:

1. Create and activate a new environment with Mamba using the environment file from the ScanNeo2 repository:

    ```bash
    mamba env create --file https://raw.githubusercontent.com/ylab-hi/ScanNeo2/main/environment.yml
    mamba activate scanneo2
    ```

    Note: ScanNeo2 requires Snakemake >= 8.x.x is not compatible with Snakemake <= 8.x.x. 

2. Deploy ScanNeo2:

```
mkdir -p /path/to/your/working/directory/
cd /path/to/your/working/directory/
git clone --recurse-submodules https://github.com/ylab-hi/ScanNeo2.git
cd ScanNeo2
```

3. (Optional) Install HLA-HD

    ScanNeo2 employs HLA-HD for HLA Class II genotyping which is required when ScanNeo2 has been configured to predict class II neoantigens. 
    Due to license reasons it has to be installed manually (download request). Please follow the instructions on the official 
    [website](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/). ScanNeo2 has been tested using HLA-HD v1.7.1


<!--

    ```bash
    mkdir -p /path/to/your/working/directory/
    cd /path/to/your/working/directory/
    snakedeploy deploy-workflow https://github.com/ylab-hi/ScanNeo2 . --tag v0.1.0
    ```
-->

3. Configure ScanNeo2 by editing the `config/config.yml` file. Make sure to adjust parameters to suit your needs and data. Note that the paths in the config file need to be relative to the directory from which you run snakemake.

### Running the Workflow

To run the workflow, use the following command:

```bash
cd /path/to/your/working/directory/
snakemake --cores all --software-deployment-method conda 
```

Note that snakemake currently does not support Mamba v2. We recommend to either use Mamba v1, or provide the option `--conda-frontend conda`.

In addition, custom configfiles can be configured using `--configfile <path/to/configfile>`. In principle, this merely 
overwrites the default config, and should include all key/value pairs of the valid config file.

For more detailed instructions and explanations on how to use ScanNeo2, please consult the [wiki](https://github.com/ylab-hi/ScanNeo2/wiki).

## Docker Support

For added convenience, we also provide a ready-to-use [Docker](https://hub.docker.com/r/yanglabinfo/scanneo2)
Container for ScanNeo2. This container encapsulates the environment required to run ScanNeo2, making it even 
easier to get started. 

## Test data

In additon, we provided test data in `.tests/integration` with configuration and resulting file that can be used 
to test the installation. 

For example, the following will run the the workflow with custom provided variants:
```
snakemake --cores all --sdm conda --configfile .tests/integration/custom-test/config/config.yaml
```
Also, the following command will run the workflow with a small raw RNA-seq dataset:
```
snakemake --cores all --sdm conda --configfile .tests/integration/indel-test/config/config.yaml
```

## System requirements

We recommend to run ScanNeo on a system with at least 64GB


## Troubleshooting

Please don't hesitate to reach out with any questions or feedback

## Citation

```bibtex
@article{Schafer2023Nov,
	author = {Sch{\ifmmode\ddot{a}\else\"{a}\fi}fer, Richard A. and Guo, Qingxiang and Yang, Rendong},
	title = {{ScanNeo2: a comprehensive workflow for neoantigen detection and immunogenicity prediction from diverse genomic and transcriptomic alterations}},
	journal = {Bioinformatics},
	volume = {39},
	number = {11},
	pages = {btad659},
	year = {2023},
	month = nov,
	issn = {1367-4811},
	publisher = {Oxford Academic},
	doi = {10.1093/bioinformatics/btad659}
}
```

## License

ScanNeo2 is licensed under MIT License.

## Contact

If you have any issues or queries about ScanNeo2, please [raise an issue on GitHub](https://github.com/ylab-hi/ScanNeo2/issues/new).
