# SEMgraph

## Overview

**SEMgraph** estimates networks and causal relations in complex systems through
Structural Equation Modeling (SEM). **SEMgraph** runs with a directed graph
(adjacency matrix) that encodes the hypothesized or data-driven causal links
among variables. It comes with the following functionalities:

- Interchangeable model representation as either an **igraph** object 
or the corresponding SEM in **lavaan** syntax. Model management functions 
include graph-to-SEM conversion, automated covariance matrix regularization, 
graph conversion to DAG, and tree (arborescence) from correlation matrices.

- Heuristic filtering, node and edge weighting, resampling and 
parallelization settings for fast fitting in case of very large models.

- Automated data-driven model building and improvement, through causal 
structure learning and bow-free interaction search and latent variable 
confounding adjustment.

- Perturbed paths finding, community searching and sample scoring, 
together with graph plotting utilities, tracing model architecture 
modifications and perturbation (i.e., activation or repression) routes.

## Installation

The latest stable version can be installed from CRAN:

``` r
install.packages("SEMgraph")
```

The latest development version can be installed from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("fernandoPalluzzi/SEMgraph")
```

## References

Grassi M, Palluzzi F, Tarantino B. **SEMgraph: an R package for causal network inference of high-throughput data with structural equation models**. Bioinformatics, 2022 Aug 30; 38(20):btac567. https://doi.org/10.1093/bioinformatics/btac567

Grassi M, Tarantino B. **SEMgsa: topology-based pathway enrichment analysis with structural equation models**. BMC Bioinformatics, 2022 Aug 17; 23(1):344. https://doi.org/10.1186/s12859-022-04884-8

Grassi M, Tarantino B. **SEMtree: tree-based structure learning methods with structural equation models**. Bioinformatics, 2023 June 09; 39(6):btad377. https://doi.org/10.1093/bioinformatics/btad377

Grassi M, Tarantino B. **SEMbap: Bow-free covariance search and data de-correlation**. PLoS Comput Biol, 2024 Sep 11; 20(9):e1012448. https://doi.org/10.1371/journal.pcbi.1012448

Grassi M, Tarantino B. **SEMdag: Fast learning of Directed Acyclic Graphs via node or layer ordering**. PLoS ONE. 2025 Jan 08; 20(1): e0317283. https://doi.org/10.1371/journal.pone.0317283
