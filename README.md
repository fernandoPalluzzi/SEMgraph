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

<p align="center">
    <img src="https://github.com/fernandoPalluzzi/SEMgraph/blob/master/docs/figures/SEMgraph_workflow.png" width=85% height=85%>
</p>

&nbsp;

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

&nbsp;

## Getting help

A gentle introduction to **SEMgraph** functionalities is available at our [**DOCs page**](https://github.com/fernandoPalluzzi/SEMgraph/edit/master/docs).

The full list of **SEMgraph** functions with examples is available at our website [**HERE**](https://fernandopalluzzi.github.io/SEMgraph/).

&nbsp;

## SEMgraph-related projects

- [**SEMgsa**](https://github.com/fernandoPalluzzi/SEMgraph/tree/master/SEMgsa_replication). SEM-based gene set analysis tool enabling perturbed pathway and gene finding by exploiting the causal structure of a graph.

- [**SEMtree**](https://github.com/fernandoPalluzzi/SEMgraph/tree/master/SEMtree). Tree structure learning methods implemented with graph and data-driven SEM-based algorithms.

- [**SEMbap**](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/SEMbap). Bow-free covariance search and data de-correlation.

- [**SEMdag**](). Fast learning of Directed Acyclic Graphs via node or layer ordering.

&nbsp;

## Available datasets

### Create updated pathway and reference network versions

**SEMgraph** and **SEMdata** reference datasets are freezed to benchmarked versions. If you would like to get the latest version of your favourite database, you can use either the R package 
[**graphite**](https://bioconductor.org/packages/release/bioc/html/graphite.html) ([**graphite tutorial**](https://bioconductor.org/packages/release/bioc/vignettes/graphite/inst/doc/graphite.pdf)), 
or our simple wrapper function, contained in the R script [**loadPathwayData.R**](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/loadPathwayData.R). The script comes with descriptions and examples.

### Latest stable datasets

**SEMgraph** offers several verified datasets to work with, for both training and research. They include (**\*\*** available with the [**SEMdata**](https://github.com/fernandoPalluzzi/SEMdata) expansion):

- [**KEGG**](https://www.genome.jp/kegg/) directed reference network of 5934 nodes and 77158 edges, derived from the **KEGG** database.

- **KEGG pathways**. A comprehensive list of 227 KEGG pathways (last update: February 2024).

- [**Reactome**](https://reactome.org) directed reference network of 9762 nodes and 416128 edges, derived from **Reactome** DB. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Reactome pathways**. A comprehensive list of 1641 pathways (last update: April 2020). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- [**STRING**](https://string-db.org/) interactome (version 10.5) of 9725 nodes and 170987 edges. [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Amyotrophic Lateral Sclerosis** (ALS) RNA-seq dataset of 139 cases and 21 healthy controls, from Tam O.H. *et al.*, 2019 (GEO accession: [**GSE124439**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Frontotemporal Dementia** (FTD) DNA methylation dataset 150 cases and 150 healthy controls, from Li Y. *et al.*, 2014 (GEO accession: [**GSE53740**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53740)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **COVID-19** RNA-seq dataset of 46 critical and 23 non-critical COVID-19 cases in young patients, from Carapito R. *et al.*, 2022 (GEO accession: [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114)). [[**\*\***](https://github.com/fernandoPalluzzi/SEMdata)]

- **Flow cytometry** data and causal model from [Sachs *et al.*, 2005](https://www.science.org/lookup/doi/10.1126/science.1105809).

&nbsp;

# References

### SEMgraph

Grassi M, Palluzzi F, Tarantino B. **SEMgraph: an R package for causal network inference of high-throughput data with structural equation models**. Bioinformatics, 2022 Aug 30; 38(20):btac567. https://doi.org/10.1093/bioinformatics/btac567

### Associated projects

Grassi M, Tarantino B. **SEMgsa: topology-based pathway enrichment analysis with structural equation models**. BMC Bioinformatics, 2022 Aug 17; 23(1):344. https://doi.org/10.1186/s12859-022-04884-8

Grassi M, Tarantino B. **SEMtree: tree-based structure learning methods with structural equation models**. Bioinformatics, 2023 June 09; 39(6):btad377. https://doi.org/10.1093/bioinformatics/btad377
