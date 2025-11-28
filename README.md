# fixed-tree-experiments

## Data
The BEAST XML configuration files are located in [ssstree_datasets](ssstree_datasets/).
Due to data sharing limitations, sequences for the SARS-CoV-2 datasets need to be downloaded from GISAID. These folders contain a `gisaid.txt` file containing the GISAID accession IDs.

## Dependencies

This project requires the following phylogenetic inference tools:

| Program | Version |
| ------------ | -------- |
| [BEAGLE]   | [4.0.1](https://github.com/beagle-dev/beagle-lib/releases/tag/v4.0.1) |
| [BEAST]    | [1.10.5pre_thorney_0.1.2](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre_thorney_v0.1.2) |
| [iqtree]  | [2.2.2.6](https://github.com/iqtree/iqtree2/releases/tag/v2.2.2.6) |
| [lsd]     |[2.3](https://github.com/tothuhien/lsd2/releases/tag/v.2.3) |
| [treetime] | [0.11.4](https://github.com/neherlab/treetime/releases/tag/v0.11.4) |


## Pipeline with docker or singularity
There is no need to install dependencies with docker or singularity.

### Running the pipeline with docker

    nextflow run main.nf -profile docker

### Running the pipeline with singularity and PBS

    nextflow -C configs/pbs.config run main.nf -profile singularity

Since the pipeline will take weeks to run to completion one should use a high performance computer. An example of configuration file for PBS Pro can be found in the [configs](configs/) folder.


## Installation

```
git clone 4ment/fixed-tree-experiments.git
```

## Install phylogenetic programs

Install **BEAST**, **IQ-TREE** and **LSD** in the `bin` directory.

## conda setup

Other Python dependencies, including **TreeTime**, are managed with Anaconda:

```
conda env create -f environment.yml
conda activate fixed-tree
```

You may need to set the `NXF_CONDA_CACHEDIR` variable for nextflow to locate the conda environment. You need to do something like this:
```
export NXF_CONDA_CACHEDIR=$HOME/miniconda3/envs
```

## Running the pipeline
```
nextflow run main.nf
```

This command creates a `results/` directory in the pipeline working directory. The directory contains all intermediate and final outputs needed to generate the figures.

## Generate figures
The results can be visualized using the included R Markdown report. Required R packages are managed with `renv`. Run once:
```
Rscript -e 'renv::restore()'
```

Then render the results to a PDF:
```
Rscript -e 'rmarkdown::render("summary.Rmd")'
```

### Generate network graphs
The network figure is generated separately using a python script. The python script dependencies are included in the conda environment *fixed-tree*.

```
python scripts/build-net.py --multi --results results --out networks.pdf
```

This script will create a pdf file containing the graphs and it will also generate a file containing z-scores contained in the manuscript.

[BEAGLE]: https://github.com/beagle-dev/beagle-lib
[BEAST]: https://github.com/beast-dev/beast-mcmc
[iqtree]: https://github.com/iqtree/iqtree2/
[lsd]: https://github.com/tothuhien/lsd2
[treetime]: https://github.com/neherlab/treetime/
