# SEMgraph
**Causal network inference and discovery** with **Structural Equation Modeling**

**SEMgraph**  Estimates networks and causal relations in complex systems through
Structural Equation Modeling (SEM). **SEMgraph** comes with the following functionalities:

- Interchangeable model representation as either an **igraph** object 
or the corresponding SEM in **lavaan** syntax. Model management functions 
include graph-to-SEM conversion, automated covariance matrix regularization, 
graph conversion to DAG, and graph creation from correlation matrices.

- Heuristic filtering, node and edge weighting, resampling and 
parallelization settings for fast fitting in case of very large models.

- Automated data-driven model building and improvement, through causal 
structure learning and bow-free interaction search and latent variable 
confounding adjustment.

- Perturbed paths finding, community searching and sample scoring, 
together with graph plotting utilities, tracing model architecture 
modifications and perturbation (i.e., activation or repression) routes.

## Installation

The latest stable version can be installed from [**CRAN**](https://CRAN.R-project.org/package=SEMgraph):

``` r
install.packages("SEMgraph")
```

The latest development version can be installed from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("fernandoPalluzzi/SEMgraph")
```

Do not forget to install the [**SEMdata**](https://github.com/fernandoPalluzzi/SEMdata) 
package too! It contains useful high-throughput sequencing data, reference networks, 
and pathways for SEMgraph training:

``` r
devtools::install_github("fernandoPalluzzi/SEMdata")
```

## Getting help

See our website [**HERE**](https://fernandopalluzzi.github.io/SEMgraph/) for help and examples.

## Available datasets

**SEMgraph** offers several datasets to work with, for both training and research. They include (**\*\*** available with the [**SEMdata**](https://github.com/fernandoPalluzzi/SEMdata) expansion):

- [**KEGG**](https://www.genome.jp/kegg/) directed reference network of 5934 nodes and 77158 edges, derived from the **KEGG** database.

- **KEGG pathways**. A comprehensive list of 225 KEGG pathways (last update: November 2021).

- [**Reactome**](https://reactome.org) directed reference network of 9762 nodes and 416128 edges, derived from **Reactome** DB. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Reactome pathways**. A comprehensive list of 1641 pathways (last update: April 2020). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- [**STRING**](https://string-db.org/) interactome (version 10.5) of 9725 nodes and 170987 edges. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Amyotrophic Lateral Sclerosis** (ALS) RNA-seq dataset of 139 cases and 21 healthy controls, from Tam O.H. *et al.*, 2019 (GEO accession: [**GSE124439**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Frontotemporal Dementia** (FTD) DNA methylation dataset 150 cases and 150 healthy controls, from Li Y. *et al.*, 2014 (GEO accession: [**GSE53740**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53740)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **COVID-19** RNA-seq dataset of 46 critical and 23 non-critical COVID-19 cases in young patients, from Carapito R. *et al.*, 2022 (GEO accession: [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Flow cytometry** data and causal model from [Sachs *et al.*, 2005](https://www.science.org/lookup/doi/10.1126/science.1105809).

&nbsp;

# SEMgraph tutorial

The following section offers an overview of **SEMgraph** functionalities. Starting from model fitting, it will gently introduce functions for model learning, weighting, clustering, and evaluation of causal effects and model perturbation. This section includes:

1. **Manual dependencies installation**
2. **Causal effects estimation, model learning, extension, and clusterinng (Amyotrophic Lateral Sclerosis dataset)**
3. **Gene Set Analysis and perturbed subnetwork/module extraction (Frontotemporal Dementia dataset)**

## 1. Manual dependencies installation

```r
# Defining dependency installers

installFromCran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

installFromBioconductor <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# SEMgraph package
devtools::install_github("fernandoPalluzzi/SEMgraph")

# SEMdata package
devtools::install_github("fernandoPalluzzi/SEMdata")

# Required packages

pkgs.cran <- c("BiocManager", "boot", "cate", "corpcor", "dagitty",
               "flip", "gdata", "ggm", "GGMncv", "glmnet",
	       "leaps", "mgcv", "pbapply", "protoclust")

pkgs.bioc <- c("graph", "RBGL", "Rgraphviz")

install_1 <- sapply(pkgs.cran, installFromCran)
install_2 <- sapply(pkgs.cran, installFromBioconductor)

# Installing suggested packages
install.packages("org.Hs.eg.db")
install.packages("huge")

options(warn = -1)

# Library loading
library(SEMgraph)
library(SEMdata)
```

&nbsp;

# References

Palluzzi F, Grassi M. **SEMgraph: An R Package for Causal Network Analysis of High-Throughput Data with Structural Equation Models**. 3 Jan 2022; arXiv:2103.08332.
