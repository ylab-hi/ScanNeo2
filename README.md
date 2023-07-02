<div align="center">
    <h1>ScanNeo2 (working title)</h1>
    <img src="https://github.com/ylab-hi/ScanNeo2/actions/workflows/linting.yml/badge.svg" alt="Workflow status badge">
</div>


## What is ScanNeo2
`Scanneo2` is a snakemake workflow for the prediction of neoantigens from 
multiple sources. In its current state, this includes 

## Getting Started

In principle, Scanneo2 aims to 


## Quickstart

Install `snakemake` and `snakedeploy`
```
mamba env create --file https://
mamba activate scanneo2
```
Deploy Scanneo2
```
mkdir -p /path/to/working/directory/
cd /path/to/working/directory/
snakedeploy deploy-workflow https://github.com/ylab-hi/scanneo2 . --tag v0.1.0
```
Run the workflow
```
snakemake --cores all --use-conda
```

Please consult the wiki for detailed instruction and explanations.


### Docker

We also provide a ready to use Docker Container that can be used 




