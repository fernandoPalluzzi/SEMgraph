# SEMgraph
**Network Analysis** and **Causal Learning** with **Structural Equation Modeling**.

# Overview

**SEMgraph** is an **R** package for network analysis and causal inference and causal structure learning through Structural Equation Modeling (SEM).
**SEMgraph** comes with the following functionalities:

- Interchangeable model representation as either an **igraph** object 
or the corresponding SEM in **lavaan** syntax. Model management functions 
include graph-to-SEM conversion, automated covariance matrix regularization, 
graph conversion to DAG, and graph creation from correlation matrices.

- Heuristic filtering, node and edge weighting, resampling and 
parallelization settings for fast fitting in case of very large models.

- Automated data-driven model building and improvement, through causal 
structure learning, bow-free interaction search, and latent variable 
confounding adjustment.

- Perturbed paths finding, community searching and sample scoring, 
together with graph plotting utilities, tracing model architecture 
modifications and perturbation (i.e., activation or repression) routes.

## Installation

The development version of **SEMgraph** can be installed in **R** with the following line:

```{r, echo = FALSE}
devtools::install_github("fernandoPalluzzi/SEMgraph")
```
Do not forget to install the [**SEMdata**](https://github.com/fernandoPalluzzi/SEMdata) package too! It contains useful high-throughput 
sequencing data, reference networks, and pathways for **SEMgraph** training:

```{r, echo = FALSE}
devtools::install_github("fernandoPalluzzi/SEMdata")
```
