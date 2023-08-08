-<div align="left">
    <h1>ScanNeo2</h1>
    <img src="https://img.shields.io/badge/snakemake-≥6.4.1-brightgreen.svg">
    <img src="https://github.com/ylab-hi/ScanNeo2/actions/workflows/linting.yml/badge.svg" alt="Workflow status badge">
</div>

## What is ScanNeo2
`Scanneo2` is a snakemake workflow for the prediction of neoantigens from multiple sources. In its current state, 
this includes canonical-splicing, exitron-splicing, gene fusion, indels and snvs.

## Getting Started

In principle, Scanneo2 aims to resolve its dependencies automatically and requires only snakemake and snakedeploy.

## Quickstart

Install `snakemake` and `snakedeploy`
```
mamba env create --file https://github.com/ylab-hi/ScanNeo2/blob/devel/environment.yml
mamba activate scanneo2
```
Deploy Scanneo2
```
mkdir -p /path/to/working/directory/
cd /path/to/working/directory/
snakedeploy deploy-workflow https://github.com/ylab-hi/scanneo2 . --tag v0.1.0
```
Configure ScanNeo2 by modifying `config/config.yml`

Run the workflow
```
cd scanneo
snakemake --cores all --use-conda
```
Please consult the [wiki](https://github.com/ylab-hi/ScanNeo2/wiki) for detailed instruction and explanations.

### Docker

We also provide a ready-to-use [Docker Container](https://hub.docker.com/r/yanglabinfo/scanneo2) 
that can be used to use Scanneo2.

## Introduction

ScanNeo2 is a comprehensive, Snakemake-based workflow designed to predict neoantigens from a variety of sources. It currently supports canonical-splicing, exitron-splicing, gene fusion, indels, and SNVs. This powerful tool provides a streamlined approach to neoantigen prediction, simplifying the process and enabling more efficient research.

## Getting Started

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

2. Deploy ScanNeo2:

    ```bash
    mkdir -p /path/to/your/working/directory/
    cd /path/to/your/working/directory/
    snakedeploy deploy-workflow https://github.com/ylab-hi/scanneo2 . --tag v0.1.0
    ```

3. Configure ScanNeo2 by editing the `config/config.yml` file. Make sure to adjust parameters to suit your needs and data.

### Running the Workflow

To run the workflow, use the following command:

```bash
cd /path/to/your/working/directory/
snakemake --cores all --use-conda
```

For more detailed instructions and explanations on how to use ScanNeo2, please consult the [wiki](https://github.com/ylab-hi/ScanNeo2/wiki).

## Docker Support

For added convenience, we also provide a ready-to-use Docker Container for ScanNeo2. This container encapsulates the environment required to run ScanNeo2, making it even easier to get started.

## Conclusion

ScanNeo2 provides an accessible, efficient method for predicting neoantigens. Its comprehensive support for multiple sources of neoantigens, along with its ease of installation and use, make it a powerful tool for researchers in the field. Please don't hesitate to reach out with any questions or feedback - we're always looking to improve ScanNeo2.

## License

ScanNeo2 is licensed under MIT License.

## Contact

If you have any issues or queries about ScanNeo2, please [raise an issue on GitHub](https://github.com/ylab-hi/ScanNeo2/issues/new) or contact us at [richard.schaefer@northwestern.edu](mailto:richard.schaefer@northwestern.edu).
