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




