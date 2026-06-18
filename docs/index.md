# What is Scanneo2

Scanneo2 is a snakemake workflow for detecting neoantigens from multiple sources. These include alternative splicing, exitron splicing, gene fusion, indels, and SNVs. User-supplied VCFs and `(wildtype, mutant)` protein-pair TSVs are also accepted as input, bypassing variant calling and VEP. See the [Configuration](https://github.com/ylab-hi/ScanNeo2/wiki/Configuration#data) page for details. 


# Installation & Usage


## Install & Deploy Workflow

Scanneo2 requires snakemake and snakedeploy which are most easily installed using the Mamba Package Manager (or by conda itself). If neither mamba nor conda is present, it is best to install it from [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html). We provide an `environment.yml` which allows creating a mamba directory that includes the above-defined dependencies. 

Simply, run
```
mamba env create --file https://raw.githubusercontent.com/ylab-hi/ScanNeo2/main/environment.yml
mamba activate scanneo2
```
Note: The exitron module in ScanNeo2 requires `apptainer` (formerly `singularity`) which needs to be installed. If this is not provided (as usually on HPC systems), the `environment_apptainer.yml` can be used for that. If the exitron module is deactivated (activate: false), `apptainer` is not required. 

In the following Scanneo2 needs to be deployed. For that, create a working directory and enter it, e.g.

```
mkdir -p /path/to/working/directory/
cd /path/to/working/directory/
```

Now the workflow needs to be deployed which can be done using

```
git clone --recurse-submodules https://github.com/ylab-hi/ScanNeo2
```

## Configuration

To configure this workflow, modify [config/config.yaml](https://github.com/ylab-hi/ScanNeo2/blob/main/config/config.yaml) according to your needs, following the explanations provided in the file. It is to be noted that each parameter can also be specified over the command line. For a detailed description of the parameters, please refer to section about the [configuration](https://github.com/ylab-hi/ScanNeo2/wiki/Configuration) in the wiki. Different locations of configuration files can be specified using `--configfile config.yaml`. 

In addition, the parameters can also be specified using the command line such as `--config parameter=value`. This is not recommended as this can be done on the top-level. 

## Usage

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda (using the [Mamba package manager](https://github.com/mamba-org/mamba) by default), run Snakemake with

```
snakemake --cores all --software-deployment-method conda
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module that has been defined by the deployment

## Testdata

Can be found at `.tests/integration`








