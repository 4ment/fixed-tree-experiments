# fixed-tree-experiments

## Dependencies

This project requires the following phylogenetic inference tools:

| Program | Version |
| ------------ | -------- |
| [BEAGLE]   | [4.0.1](https://github.com/beagle-dev/beagle-lib/releases/tag/v4.0.1) |
| [BEAST]    | [1.10.5pre_thorney_0.1.2](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre_thorney_v0.1.2) |
| [iqtree]  | [2.2.2.6](https://github.com/iqtree/iqtree2/releases/tag/v2.2.2.6) |
| [lsd]     |[2.3](https://github.com/tothuhien/lsd2/releases/tag/v.2.3) |
| [treetime] | [0.11.4](https://github.com/neherlab/treetime/releases/tag/v0.11.4) |



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

You may need to set the `NXF_CONDA_CACHEDIR` variable for nextflow to locate the conda environment
```
export NXF_CONDA_CACHEDIR=$HOME/miniconda3/envs
```

## Running the pipeline
```
nextflow run main.nf
```

## Generate figures
The results can be visualized using the included R Markdown report. Required R packages are managed with `renv`. Run once:
```
Rscript -e 'renv::restore()'
```

Then render the results to a PDF:
```
Rscript -e 'rmarkdown::render("plot.Rmd")'
```


[BEAGLE]: https://github.com/beagle-dev/beagle-lib
[BEAST]: https://github.com/beast-dev/beast-mcmc
[iqtree]: https://github.com/iqtree/iqtree2/
[lsd]: https://github.com/tothuhien/lsd2
[treetime]: https://github.com/neherlab/treetime/
