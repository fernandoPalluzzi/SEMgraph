# SEMgraph
**Causal network inference and discovery** with **Structural Equation Modeling**

**SEMgraph**  Estimate networks and causal relations in complex systems through
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

- **KEGG pathways**. A comprehensive list of 306 KEGG pathways (last update: April 2020). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- [**Reactome**](https://reactome.org) directed reference network of 9762 nodes and 416128 edges, derived from the **Reactome** database. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Reactome pathways**. A comprehensive list of 1641 pathways (last update: April 2020). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- [**STRING**](https://string-db.org/) interactome (version 10.5) of 9725 nodes and 170987 edges. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Amyotrophic Lateral Sclerosis** (ALS) RNA-seq dataset of 139 cases and 21 healthy controls, from Tam O.H. *et al.*, 2019 (GEO accession: [**GSE124439**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Frontotemporal Dementia** (FTD) DNA methylation dataset 150 cases and 150 healthy controls, from Li Y. *et al.*, 2014 (GEO accession: [**GSE53740**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53740)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Flow cytometry** data and causal model from [Sachs *et al.*, 2005](https://www.science.org/lookup/doi/10.1126/science.1105809).

## Coming soon

Next versions of SEMgraph will include new updated databases, new functionalities for *de novo* (data-driven) causal model learning, and new inference methods.

## References

Palluzzi F, Grassi M. **SEMgraph: An R Package for Causal Network Analysis of High-Throughput Data with Structural Equation Models**. Sep 2021; arXiv:2103.08332.
