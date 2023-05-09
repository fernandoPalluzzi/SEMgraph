# SEMgraph tutorial

The following section offers an overview of **SEMgraph** functionalities. Starting from model fitting, it will gently introduce functions for model learning, weighting, clustering, and evaluation of causal effects and model perturbation. This section includes:

1. **Supplementary packages installation**
2. **Causal effects estimation, model learning, extension, and clusterinng (Amyotrophic Lateral Sclerosis dataset)**
3. **Gene Set Analysis and perturbed subnetwork/module extraction (Frontotemporal Dementia dataset)**

Please, also visit our website [**HERE**](https://fernandopalluzzi.github.io/SEMgraph/) for a complete list of SEMgraph functions and related examples.

&nbsp;

## 1. Supplementary packages.

Besides the required packages, SEMgraph suggests the use of org.Hs.eg.db, for gene ID conversion.
SEMgraph uses entrez IDs to avoid special chatacters (such as hyphens or slashes), but it can use official gene symbols as labels.

```r
install.packages("org.Hs.eg.db")
```

Sometimes, it could be useful to relax Gaussianity constraints when fitting a SEM.
The huge package does it by applying a nonparanormal transform (PMID: 26834510).

```r
install.packages("huge")
```

&nbsp;

## 2. Amyotrophic Lateral Sclerosis (ALS) data analysis.

### 2.1. The ALS dataset.

**SEMdata** provides the ALS RNA-seq dataset of 139 cases and 21 healthy controls, from Tam O.H. *et al.*, 2019 (GEO accession: [GSE124439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439)). Raw data were pre-processed applying batch effect correction, using the sva R package (Leek et al., 2012), to remove data production center and brain area biases. Using multidimensional scaling-based clustering, ALS-specific and HC-specific clusters were generated. Misclassified samples were blacklisted and removed from the dataset. Since the expression of many genes is significantly different from Gaussian, we apply a nonparanormal transform with the **huge** package, to relax the normality assumption.

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

**Figure 1. ALS model fitting.** Estimated group effects on nodes and direct effects (red/pink: activation, blue/lightblue: repression).

![alt text](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/docs/figures/Figure1.png)

### 2.3. Total effect estimation as Average Causal Effect (ACE).

Suppose having a path P = X -> M1 -> ... -> Mk -> Y between two nodes X and Y, connected through k mediators. The total effect (TE = DE + IE) is the sum of the direct effect X -> Y, DE = b(X,Y), and the indirect effect through the mediators, IE = b(X,M1)\*b(M1,M2)\*...\*b(Mk,Y).
One convenient way of estimating the TE is through the definition of ACE by [Pearl J, 1998](https://doi.org/10.1177/0049124198027002004). The simplest estimation of the TE as ACE is possible in directed acyclic graphs (DAGs), thorugh linear regression. The parent set pa(X) of X blocks all backdoor (i.e., confounding) paths from X to Y, and the ACE is equal to the coefficient b(X,Y) in a multiple regression Y ~ X + pa(X).

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

#  The pathFinder function extracts all directed paths, fits them, and evaluates their perturbation

paths <- pathFinder(alsData$graph, data.npn, alsData$group, ace = ace)
print(paths$dfp)
```

### 2.4. Model estimation strategies.

The input graph is a picture of the current knowledge. Besides the evaluation of known relationships, data can be used to infer new ones. Four model search strategies are implemented in the `modelSearch` function. The *basic* one is completely data-driven and requires an input graph only to establish the initial topological order. In the *direct* strategy, The input graph structure is improved through direct (i.e., adjacent) link search, followed by reference-based ("gnet" argument) interaction validation and import from the reference network, with no mediators (d = 1). In the *outer* strategy, new interactions and connectors (i.e., mediators) will be searched and imported from the reference network. This strategy is analogous to the "outer" one, but disables external mediator search. In other words, new indirect paths are generated by adding new interactions between nodes only from the input model.

```r
# Model Search with search = "basic" and beta = 0.1
# The beta argument determines the minimum absolute LASSO coefficient
# for an edge to be included in the output model.
# Reducing beta (up to 0) will increase model complexity.

# This strategy is data-driven (gnet = NULL) and 
# only direct connections are added (d = 0; i.e., no new mediators).

model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = NULL, d = 0,
		     search = "basic",
		     beta = 0.1,
		     method = "BH",
		     alpha = 0.05,
		     verbose = TRUE)


## Figure 2. ALS improved model (basic strategy). --------------------##

V(model$graph)$label <- mapIds(org.Hs.eg.db, V(model$graph)$name,
                               column = 'SYMBOL',
                               keytype = 'ENTREZID')

pert <- SEMrun(model$graph, model$data, alsData$group)

path <- SEMpath(model$graph, model$data, alsData$group,
                from = "6647",
                to = "4741",
                path = "directed",
                verbose = TRUE)

pdf("Figure2.pdf", width = 14, height = 9)
par(mfrow=c(3,1), mar=c(1,1,1,1))
gplot(model$graph, main = "\nA)  ALS model structure", cex.main = 3, fontsize = 40)
gplot(pert$graph, main = "\nB)  ALS model perturbation", cex.main = 3, fontsize = 40)
gplot(path$map, main = "\nC)  SOD1-NEFM path", cex.main = 3, fontsize = 40)
dev.off()

##--------------------------------------------------------------------##


# Other possible strategies (not shown in the manuscript).

# Direct strategy.
# - Knowledge-based estimation (gnet should be a directed reference network).
# - Only direct interactions are inferred (d is fixed to 1).
# - New interactions are inferred from data, during the execution.
# - A smaller starting beta value is suggested (beta = 0.05)

model1 <- modelSearch(graph = alsData$graph, data = data.npn,
                      gnet = kegg, d = 1,
		      search = "direct",
		      beta = 0.05,
		      method = "BH",
		      alpha = 0.05,
		      verbose = FALSE)

# Inner strategy.
# - The reference network is used to validate new interactions and mediators.
# - Inferred mediators must already belong to the input graph.
# - Larger d values increase model complexity (suggested: d = 2).

model2 <- modelSearch(graph = alsData$graph, data = data.npn,
                      gnet = kegg, d = 2,
		      search = "inner",
		      beta = 0.05,
		      method = "BH",
		      alpha = 0.05,
		      verbose = FALSE)
					 
# Outer strategy.
# - Knowledge-based estimation (gnet should be a directed reference network).
# - Up to d - 1 mediators can be imported from the reference network.
# - Larger d values increase model complexity (suggested: d = 2).

model3 <- modelSearch(graph = alsData$graph, data = data.npn,
                      gnet = kegg, d = 2,
		      search = "outer",
		      beta = 0.05,
		      method = "BH",
		      alpha = 0.05,
		      verbose = FALSE)

# graph models visualization:	 
par(mfrow=c(2,2), mar=rep(2,4))
gplot(model$graph, main="basic")
gplot(model1$graph, main="direct")
gplot(model2$graph, main="inner")
gplot(model3$graph, main="outer")

```

**Figure 2. ALS improved model (basic strategy).** Data-driven model improvement of the input ALS model: (A) added connectors in green; (B) improved model fitting and perturbation evaluation (red/pink: activation, blue/lightblue: repression); (C) SOD1-NEFM directed path.

![alt text](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/docs/figures/Figure2.png)

### 2.5. Communities and factor scores.

The modular structure of biological networks could often reveal local effects and perturbed routes and communities hidden within a larger and more complex context. **SEMgraph** allows to detect and estimate these local properties as follows.

```r
# Improved ALS model (model$graph) clustering and scoring, using alatent variable
# "hidden" model (LV), edge betweeness clustering (EBC) algorithm, and a minimum 
# cluster size of 5 nodes.
# Other clustering algorithms can be exploited (e.g., the walktrap community 
# detection algorithm, WTC) to improve the interpretation of results.

LV <- clusterScore(model$graph, model$data, alsData$group,
                   type = "ebc",
		   HM = "LV",
		   size = 5)
				  
table(LV$membership)

head(parameterEstimates(LV$fit))


# Clustering only (no scores calculation)

C <- clusterGraph(model$graph, type = "ebc", HM = "LV", size = 5, verbose = TRUE)

# Cluster plot utility

cg <- cplot(graph = model$graph, membership = LV$membership, verbose = TRUE)
list(cg)
gplot(cg$graph)


# Cluster extraction, fitting, and perturbation evaluation

cls <- extractClusters(graph = model$graph, data = model$data, group = alsData$group,
                       membership = LV$membership)
print(cls$dfc)


## Figure 3. Edge betweenness clusters (EBC) mapped over the improved model. --##

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(cg$graph)$label <- mapIds(org.Hs.eg.db, V(cg$graph)$name, 'SYMBOL', 'ENTREZID')

# Set node colors
V(cg$graph)$color[V(cg$graph)$color == 2] <- "lightsalmon"
V(cg$graph)$color[V(cg$graph)$color == 3] <- "lightgreen"
V(cg$graph)$color[V(cg$graph)$color == 4] <- "lightyellow"

# Graph plot 3
pdf("Figure3.pdf", width = 16, height = 8)
gplot(cg$graph, fontsize = 30)
dev.off()

##-----------------------------------------------------------------------------##
```

&nbsp;

## 3. Frontotemporal Dementia (FTD) data analysis.

The FTD dataset coming with **SEMdata** is a data matrix of 256 rows (subjects; 150 FTD patients and 150 healthy controls) and 16560 columns (genes) containing the value of the first principal component of DNAme levels, obtained applying a principal component analysis to methylated CpG sites within the promoter region, for each gene (genes with unmethylated CpGs in both conditions were discarded). This dataset was derived from the study by Li Y. *et al.*, 2014 (GEO accession: [GSE53740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53740)).

```r
# Library loading

library(SEMgraph)
library(SEMdata)
library(huge)

# Nonparanormal transform of DNAme PC1 data

pc1.npn <- huge.npn(ftdDNAme$pc1)

dim(pc1.npn)

# Defining groups

group <- ftdDNAme$group

table(group)
```

### 3.1. Gene Set Analysis (GSA).

In absence of a conceptual model (e.e., one built from an expert's indication), the input graph can be inferred from data, literature or both of them. In the following example, we will take advantage of our FTD dataset and known FTD-associated pathways from the KEGG database. To this end, GSA can be used to assess the actual perturbation of known pathways, given data, and extract those genes (i.e., seeds) underlying perturbed routes. Seeds can also be defined using knowledge, such as importing a list of known disease-associated genes, or from other data sources (e.g. mutational, transcriptional, or epigenetic data).

```r
# Known FTD-related pathway selection from KEGG

ftd.pathways <- c("MAPK signaling pathway",
                  "Protein processing in endoplasmic reticulum",
                  "Endocytosis",
                  "Wnt signaling pathway",
                  "Notch signaling pathway",
                  "Neurotrophin signaling pathway")

# Pathway extraction from the kegg.pathways object

j <- which(names(kegg.pathways) %in% ftd.pathways)

# Gene set analysis (GSA) with 5000 permutations

ftd.gsa <- SEMgsa(kegg.pathways[j], pc1.npn, group, n_rep = 5000)
print(ftd.gsa$gsa)

# Seed extraction

seed <- unique(unlist(ftd.gsa$DEG))
length(seed)
```

### 3.2. Network weighting and perturbed backbone extraction.

While seeds pinpoint nodes with specific properties, edges can be weighted on the base of quantitative and phenotype data to define preferential ways of information exchange (e.g., perturbation propagation) through the network. **SEMgraph** offers different alternatives to generate data-driven weights, but the fastest of them is based on the Fisher's "r-to-z" method, to test the group difference between correlation coefficients of pairs of interacting nodes (Fisher, 1915). Both seeds and weights can be used to extract the perturbed core(s) of a network to highlight its disease-associated components.

```r
# KEGG interactome weighting

W <- weightGraph(kegg, pc1.npn, group, method = "r2z")
summary(W)

# Perturbed backbone as a Steiner tree (Kou's algorithm)
# A Steiner tree is the minimum cost (distance) graph
# including all the specified seeds.

ST <- SEMtree(W, data = pc1.npn, seed = seed, type = "ST", eweight = "pvalue")
summary(ST)

# Perturbation evaluation with raw data

sem1 <- SEMrun(ST, pc1.npn, group)

# Data deconfounding and perturbation evaluation

# When a DAG is used for causal network inference, missing edges 
# are often masked by unmeasured confounding variables.
# This step might reduce this effect through Bow-free Acyclic Path (BAP) 
# search and data deconfounding.
adj.pc1 <- SEMbap(ST, pc1.npn, method = "bonferroni", alpha = 5E-06)$data

adj.sem1 <- SEMrun(ST, adj.pc1, group)

# Tree agglomerative hierarchical clustering (TAHC)
# This optional step allow us to detect possible communities
# within our tree.

C <- clusterGraph(ST, type = "tahc")
cg <- cplot(ST, membership = C)
list(cg)


## Figure 4. FTD perturbed backbone. ---------------------------------##

# Set node color and sze

cg2 <- cg$HM2

V(cg2)$color <- ifelse(V(cg2)$name %in% seed, "green", "white")
#gplot(cg2, l = "neato")

V(cg2)$size <- 3*degree(cg2, mode="total")
#gplot(cg2, l = "neato")

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(cg2)$label <- mapIds(org.Hs.eg.db, V(cg2)$name, 'SYMBOL', 'ENTREZID')
#gplot(cg2, l = "neato")

pdf("Figure4.pdf", width = 16, height = 12)
gplot(cg2, l = "neato", fontsize = 30)
dev.off()

##--------------------------------------------------------------------##
```

### 3.3. Locating differentially connected genes.

The SEMgraph differential causal inference (DCI) module enables the detection of perturbed nodes and edges for large graphs. This module is useful when the aim is to find a perturbed backbone of essential disease-associated nodes nodes.

```r
# Input graph as the union of FTD KEGG pathways

gU <- graph.union(kegg.pathways[j])
gU <- properties(gU)[[1]]
summary(gU)

# SEM-based differential causal inference (DCI) with the EBC clustering algorithm

gD <- SEMdci(gU, pc1.npn, group, type = "ebc", method = "BH", alpha = 0.05)
summary(gD)

# Perturbation evaluation of the 3th component of gD

gC <- properties(gD)
gC3 <- gC[[3]]
sem1 <- SEMrun(gC3, pc1.npn, group, fit = 1, algo = "cggm")
sem2 <- SEMrun(sem1$graph, pc1.npn, group, fit = 2, algo = "cggm")


## Figure 5. DCI perturbed backbone. ---------------------------------##

# Set node and edge colors
gC3 <- sem2$graph

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(gC3)$label <- mapIds(org.Hs.eg.db, V(gC3)$name, "SYMBOL", "ENTREZID")

png("Figure5.png", width = 16, height = 12, units = 'in', res = 400)
gplot(gC3, l = "fdp", fontsize = 30)
dev.off()

##--------------------------------------------------------------------##
```

**Figure 5. DCI perturbed backbone.** Figure showing the DCI perturbed backbone of the original FTD model (red/pink: activation, blue/lightblue: repression).

![alt text](https://github.com/fernandoPalluzzi/SEMgraph/blob/master/docs/figures/Figure5.png)

&nbsp;

# References

Palluzzi F, Grassi M. **SEMgraph: An R Package for Causal Network Analysis of High-Throughput Data with Structural Equation Models**. 3 Jan 2022; arXiv:2103.08332.
