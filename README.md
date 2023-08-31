# fixed-tree-experiments

## Dependencies

- iqtree2
- lsd2
- nextflow
- BEAGLE
- BEAST

## Installation

```
git clone 4ment/fixed-tree-experiments.git
```

## Install dependencies with anaconda

```
conda env create -f environment.yml
conda activate fixed-tree
```

You may need to set the `NXF_CONDA_CACHEDIR` variable for nextflow to locate the conda envirnoment
```
export NXF_CONDA_CACHEDIR=$HOME/miniconda3/envs
```

## Running the pipeline
```
nextflow run main.nf --datasets ssstree_datasets
```