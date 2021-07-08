# SEMgraph
**Network Analysis** and **Causal Learning** with **Structural Equation Modeling**

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

The latest stable version can be installed from [**CRAN**](https://cran.r-project.org/web/packages/SEMgraph/index.html):

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

## Coming soon

Next versions of SEMgraph will include new functionalities for *de novo* (data-driven) model building. A draft version of this functionality is implemented in the unstable release 1.0-beta2 "Parallel Universes".

## References

Palluzzi F, Grassi M. **SEMgraph: An R Package for Causal Network Analysis of High-Throughput Data with Structural Equation Models**. Mar 2021; arXiv:2103.08332.
