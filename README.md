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

## 2. Amyotrophic Lateral Sclerosis (ALS) data analysis.

### 2.1. The ALS dataset.

**SEMdata** provides the ALS RNA-seq dataset of 139 cases and 21 healthy controls, from Tam O.H. *et al.*, 2019 (GEO accession: GSE124439). Raw data were pre-processed applying batch effect correction, using the sva R package (Leek et al., 2012), to remove data production center and brain area biases. Using multidimensional scaling-based clustering, ALS-specific and HC-specific clusters were generated. Misclassified samples were blacklisted and removed from the dataset. Since the expression of many genes is significantly different from Gaussian, we apply a nonparanormal transform with the **huge** package, to relax the normality assumption.

```r
# ALS input graph
summary(alsData$graph)

# ALS RNA-seq expression data
dim(alsData$exprs)

# group = {1: case, 0: control} vector
table(alsData$group)

# Nonparanormal transform
library(huge)
data.npn <- huge.npn(alsData$exprs)
```

### 2.2. Model fitting.

**SEMgraph** offers three main modes of model fitting: (i) common model fitting, (ii) node perturbation evaluation, and (iii) edge perturbation evaluation. **SEMgraph** will automatically take care of applying shrinkage methods in case of high dimensionality (#variables >> #subjects), heuristics and parallelization settings for large graphs. In addition, perturbation evaluation enables the extraction of differentially regulated nodes (DRNs) and edges (DREs).

```r
# ALS model fitting (sem0: common model, no groups)
# The whole dataset is used to fit the model and perturbation is not evaluated (group = NULL)
sem0 <- SEMrun(graph = alsData$graph, data = data.npn)

est0 <- parameterEstimates(sem0$fit)
head(est0)

# Other possible output (not shown in the manuscript)
summary(sem0$fit)


# ALS model fitting
# sem1: common model, group influence on nodes (node perturbation)

sem1 <- SEMrun(graph = alsData$graph, data = data.npn, group = alsData$group)

est1 <- parameterEstimates(sem1$fit)
head(est1)

# Other possible output (not shown in the manuscript).
summary(sem1$fit)


## Figure 1. Estimated group effects on nodes and direct effects. ----##

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(sem1$graph)$label <- mapIds(org.Hs.eg.db, V(sem1$graph)$name,
                              column = 'SYMBOL',
                              keytype = 'ENTREZID')
# Graph plot
pdf("Figure1.pdf", width = 14, height = 9)
gplot(sem1$graph, cex.main = 2.5, fontsize = 22)
dev.off()

##--------------------------------------------------------------------##


# RICF fitting
# Fast fitting for large graphs (SE estimation disabled)

ricf1 <- SEMrun(alsData$graph, data.npn, alsData$group, algo = "ricf")

# Other possible output (not shown in the manuscript)
summary(ricf1$fit)
print(ricf1$gest)


# ALS model fitting
# sem2: one model for eaach group, group influence on edges (edge perturbation)

sem2 <- SEMrun(alsData$graph, data.npn, alsData$group, fit = 2)

# Other possible output (not shown in the manuscript).
summary(sem2$fit)


# Perturbed graph elements

# Differentially Regualted Nodes (DRNs)
DRN <- sem1$gest[sem1$gest$pvalue < 0.05,]
nrow(DRN)
head(DRN)

# Differentially Regulated Edges (DREs)
DRE <- sem2$dest[sem2$dest$pvalue < 0.05,]
nrow(DRE)
head(DRE)


# CGGM fitting
# Fast edge perturbation calculation for large graphs
cggm2 <- SEMrun(alsData$graph, data.npn, alsData$group, fit = 2, algo = "cggm")

# Other possible output (not shown in the manuscript)
summary(cggm2$fit$Group_0)
summary(cggm2$fit$Group_1)
print(cggm2$dest)
```

### 2.3. Total effect estimation as Average Causal Effect (ACE).

In its general definition, the total effect (TE) between two nodes x and y in a path P = x -> ... -> y is the sum of the DE x -> y (i.e., the "beta" coefficient) and the IE = b(1,2)<sub>*;</sub>b(2,3)<sub>*;</sub>...<sub>*;</sub>b(k-1,k).
One convenient way of estimating the TE is through the definition of ACE by [Pearl J, 1998](https://doi.org/10.1177/0049124198027002004). The simplest estimation of the TE as ACE is possible in DAGs, thorugh linear regression. The parent set pa(X) of X blocks all backdoor (i.e., confounding) paths from X to Y, and the ACE is equal to the b(Y,X|Z) coefficient in a multiple regression Y ~ X + pa(X).

```r
# Average Causal Effect (ACE)
# optimal adjustement set (type = "optimal")
# source -> sink ACE (effect = "source2sink");
# Benjamini-Hochberg multiple test adjustment (method = "BH");
# 5% significance level (alpha = 0.05)

ace <- SEMace(graph = alsData$graph, data = data.npn,
			  type = "optimal", effect = "source2sink",
			  method = "BH", alpha = 0.05)


# Sort by decreasing abs(z) value
ace <- ace[order(abs(ace$z), decreasing = TRUE),]
nrow(ace)
head(ace)


# SOD1-CASP3 path extraction and fitting

source <- as.character(ace$source[6])
sink <- as.character(ace$sink[6])
path <- SEMpath(alsData$graph, data.npn, alsData$group,
                from = source, to = sink, 
                path = "directed",
                verbose = TRUE)


# All directed paths extraction, fit, and perturbation evaluation

paths <- pathFinder(alsData$graph, data.npn, alsData$group, ace = ace)
print(paths$dfp)
```

&nbsp;

# References

Palluzzi F, Grassi M. **SEMgraph: An R Package for Causal Network Analysis of High-Throughput Data with Structural Equation Models**. 3 Jan 2022; arXiv:2103.08332.
