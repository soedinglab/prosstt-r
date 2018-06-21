## PROSSTT and prosstt

Single-cell RNAseq is revolutionizing cellular biology, and many algorithms are developed for the analysis of scRNAseq data. Simulations provide an easy way to test the performance of trajectory inference methods on realistic data with a known "gold standard".

PROSSTT (PRObabilistic Simulations of ScRNA-seq Tree-like Topologies) is a package with code for the simulation of scRNAseq data for complex dynamic processes such as cell differentiation. PROSSTT is [open source GPL-licensed software implemented in Python](https://github.com/soedinglab/prosstt).

`prosstt` is a lightweight R package with the methods necessary to evaluate predicted lineage trees of simulated data. Additionally, for simulations produced by PROSSTT, this package contains the necessary I/O code to read parameter files and visualize simulations and predictions.

## Installation

You can install the latest release version from CRAN (pending):

`install.packages("prosstt", dependencies = TRUE)`

Alternatively, you can install directly from github using `devtools`:

```
library(devtools)
install_github("soedinglab/prosstt-r")
```

## Dependencies

`prosstt` imports the following packages:

* `igraph`, for convenient graph representations and graph utility functions,
* `infotheo`, for the calculation of mutual information,
* `readr`, for some convenient I/O operations,
* `viridis`, for the sequential colormaps needed for pseudotime,
* `LSD`, for the discrete, divergint colormaps needed for branches.

## How to use

A notebook with a detailed tutorial is [provided](https://soedinglab.github.io/prosstt-r/example/).