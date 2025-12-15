# SEMgraph

## Overview

**SEMgraph**  Estimate causal relations in network or in complex systems with **Structural Equation Modeling** (SEM) using as input a **directed graph** that encodes the hypothesized or data-driven causal relationships among variables, a **data matrix** with *n* samples and *p* variables, and (optional) a  binary **group** vector of experimental conditions for the *n* samples.
**SEMgraph** comes with the following functionalities (see Figure):
- Interchangeable model representation as either an [igraph](https://igraph.org/) object or the corresponding SEM in [lavaan](https://lavaan.ugent.be/) syntax. Model management functions include graph-to-SEM conversion, automated covariance matrix regularization, graph conversion to DAG (Directed Acyclic Graph), and tree (arborescence) from correlation matrices.
- Heuristic filtering, node and edge weighting, resampling and parallelization settings for fast fitting in case of very large models.
- Automated data-driven model building and improvement, through causal structure learning and bow-free interaction search and latent variable confounding adjustment.
- Perturbed paths finding, community searching and sample scoring, together with graph plotting utilities, tracing model architecture modifications and perturbation (i.e., activation or repression) routes.

<p align="center">
   <img src="README.png" width=85% height=85%>
</p>
&nbsp;

## Installation

The latest stable version can be installed from [**CRAN**](https://CRAN.R-project.org/package=SEMgraph):
``` r
install.packages("SEMgraph")
```

The latest development version can be installed from GitHub:
``` r
devtools::install_github("fernandoPalluzzi/SEMgraph")
```
## Getting Started

The full list of **SEMgraph** functions with examples and a tutorial is available [**HERE**](https://fernandopalluzzi.github.io/SEMgraph/).
&nbsp;

## Companion packages

(Optional) install also:

[**SEMdata**](https://github.com/fernandoPalluzzi/SEMdata) It contains useful high-throughput sequencing data, reference networks, and pathways for SEMgraph training:

``` r
devtools::install_github("fernandoPalluzzi/SEMdata")
```

[**SEMdeep**](https://github.com/BarbaraTarantino/SEMdeep). It provides a SEM-based framework using **machine learning (ML)** and **deep neural network (DNN)** algorithms. 

``` r
install.packages("SEMdeep")
```

## Associated projects with SEMgraph functions:

- [**SEMgsa**](https://github.com/fernandoPalluzzi/SEMgraph/tree/master/SEMgsa). SEM-based gene set analysis tool enabling perturbed pathway and gene finding by exploiting the causal structure of a graph.
- [**SEMtree**](https://github.com/fernandoPalluzzi/SEMgraph/tree/master/SEMtree). Tree structure learning methods implemented with graph and data-driven SEM-based algorithms.
- [**SEMbap**](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/SEMbap). Bow-free covariance search and data de-correlation.
- [**SEMdag**](https://github.com/fernandoPalluzzi/SEMgraph/tree/master/SEMdag). Fast learning of Directed Acyclic Graphs via node or layer ordering.
&nbsp;

## References

### SEMgraph

Grassi M, Palluzzi F, Tarantino B. **SEMgraph: an R package for causal network inference of high-throughput data with structural equation models**. Bioinformatics, 2022 Aug 30; 38(20):btac567. https://doi.org/10.1093/bioinformatics/btac567

### Associated projects

Grassi M, Tarantino B. **SEMgsa: topology-based pathway enrichment analysis with structural equation models**. BMC Bioinformatics, 2022 Aug 17; 23(1):344. https://doi.org/10.1186/s12859-022-04884-8

Grassi M, Tarantino B. **SEMtree: tree-based structure learning methods with structural equation models**. Bioinformatics, 2023 June 09; 39(6):btad377. https://doi.org/10.1093/bioinformatics/btad377

Grassi M, Tarantino B. **SEMbap: Bow-free covariance search and data de-correlation**. PLoS Comput Biol, 2024 Sep 11; 20(9):e1012448. https://doi.org/10.1371/journal.pcbi.1012448

Grassi M, Tarantino B. **SEMdag: Fast learning of Directed Acyclic Graphs via node or layer ordering**. PLoS ONE. 2025 Jan 08; 20(1): e0317283. https://doi.org/10.1371/journal.pone.0317283
