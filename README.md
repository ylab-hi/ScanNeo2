<div align="left">
    <h1>ScanNeo2</h1>
    <img src="https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg">
    <img src="https://github.com/ylab-hi/ScanNeo2/actions/workflows/linting.yml/badge.svg" alt="Workflow status badge">
</div>

## What is ScanNeo2
`Scanneo2` is a snakemake workflow for the prediction of neoantigens from multiple sources. In its current state, 
this includes canonical-splicing, exitron-splicing, gene fusion, indels and snvs.

## Getting Started

In principle, Scanneo2 aims to resolve its dependencies automatically and requires only snakemake and snakedeploy.

### Prerequisites

Before installing and using ScanNeo2, make sure you have the following software installed:

1. [Mamba](https://github.com/conda-forge/miniforge#mambaforge): An open-source package manager. Mamba should be installed independently by the user.

### Installation

To get started with ScanNeo2, follow the steps below:

1. Create and activate a new environment with Mamba using the environment file from the ScanNeo2 repository:

    ```bash
    mamba env create --file https://raw.githubusercontent.com/ylab-hi/ScanNeo2/main/environment.yml
    mamba activate scanneo2
    ```

    Note: This installs Snakemake v7.32.x. In its current form, ScanNeo2 is not comptabile with Snakemake >= 8.x.x. 
    If ScanNeo2 is configured to use the exitron module, singularity needs to be installed. For that, the 
    `environment_singularity.yml' can be used. However, most HPC servers provide their own module installation.

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
    [website](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/). ScanNeo2 has been tested using HLA-HD v1.7.0


<!--

    ```bash
    mkdir -p /path/to/your/working/directory/
    cd /path/to/your/working/directory/
    snakedeploy deploy-workflow https://github.com/ylab-hi/ScanNeo2 . --tag v0.1.0
    ```
-->

3. Configure ScanNeo2 by editing the `config/config.yml` file. Make sure to adjust parameters to suit your needs and data.

### Running the Workflow

To run the workflow, use the following command:

```bash
cd /path/to/your/working/directory/
snakemake --cores all --use-conda
```

As mentioned above, when exitron detection is activated the singularity option `--use-singularity` has to be used as well.

```bash
snakemake --cores all --software-deployment-method conda apptainer 
```

In addition, custom configfiles can be configured using `--configfile <path/to/configfile>`. In principle, this merely 
overwrites the default config, and should include all key/value pairs of the valid config file.

For more detailed instructions and explanations on how to use ScanNeo2, please consult the [wiki](https://github.com/ylab-hi/ScanNeo2/wiki).

## Docker Support

For added convenience, we also provide a ready-to-use [Docker](https://hub.docker.com/r/yanglabinfo/scanneo2)
Container for ScanNeo2. This container encapsulates the environment required to run ScanNeo2, making it even 
easier to get started. 

## Test data

In additon, we provided test data in `.tests/integration` with configuration and resulting file that can be used to test the installation

## System requirements

We recommend to run ScanNeo on a system with at least 64GB


## Conclusion

ScanNeo2 provides an accessible, efficient method for predicting neoantigens. Its comprehensive support for multiple sources of neoantigens, along with its ease of installation and use, make it a powerful tool for researchers in the field. Please don't hesitate to reach out with any questions or feedback - we're always looking to improve ScanNeo2.

## Citation

If ScanNeo2 has proven useful in your work please cite it using the linked publication.

## License

ScanNeo2 is licensed under MIT License.

## Contact

If you have any issues or queries about ScanNeo2, please [raise an issue on GitHub](https://github.com/ylab-hi/ScanNeo2/issues/new).
