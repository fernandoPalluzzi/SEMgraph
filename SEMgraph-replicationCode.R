### SEMgraph: An R Package for Causal Network Analysis of High-Throughput 
### Data with Structural Equation Models

                        # Replication code #


### Installation ###

# SEMgraph URL: https://github.com/fernandoPalluzzi/SEMgraph

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

pkgs.cran <- c("BiocManager", "cate", "corpcor", "dagitty", "diffusr",
               "flip", "gdata", "GGMncv", "glmnet", "pbapply", "pcalg",
               "protoclust", "RcppEigen")

pkgs.bioc <- c("graph", "ggm", "Rgraphviz")

install_1 <- sapply(pkgs.cran, installFromCran)
install_2 <- sapply(pkgs.cran, installFromBioconductor)

# Installing suggested packages
install.packages("org.Hs.eg.db")
install.packages("huge")



###  Section 3. The SEMgraph package.  ###


## 3.1. Getting started with SEMgraph: SEM fitting functions.

# Environment preparation

library(SEMgraph)
library(SEMdata)

graph <- properties(kegg.pathways$"Amyotrophic lateral sclerosis (ALS)")[[1]]

library(huge)

data.npn <- huge.npn(alsData$exprs)

# ALS model fitting (sem0: common model, no groups)

sem0 <- SEMrun(graph = alsData$graph, data = data.npn)

est <- parameterEstimates(sem0$fit)
head(est)

# ALS model fitting (sem1: common model, group influence on nodes)

sem1 <- SEMrun(graph = alsData$graph, data = data.npn, group = alsData$group)

est1 <- parameterEstimates(sem1$fit)
head(est1)


## Figure 2. Estimated group effects on nodes and direct effects. ----##

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(sem1$graph)$label <- mapIds(org.Hs.eg.db, V(sem1$graph)$name,
                              column = 'SYMBOL',
                              keytype = 'ENTREZID')
# Graph plot
pdf("Figure2.pdf", width = 14, height = 9)
gplot(sem1$graph, cex.main = 2.5, fontsize = 22)
dev.off()

##--------------------------------------------------------------------##


# RICF approximation

ricf1 <- SEMrun(alsData$graph, data.npn, alsData$group, algo = "ricf")
summary(ricf1$fit)

# Global perturbation effects (D, A, E)

head(ricf1$gest)

# ALS model fitting (sem2: two models, group influence on edges)

sem2 <- SEMrun(alsData$graph, data.npn, alsData$group, fit = 2)

est1 <- parameterEstimates(sem1$fit)
head(est1)

# Perturbed graph elements

# Differentially Regualted Nodes (DRNs)
DRN <- sem1$gest[sem1$gest$pvalue < 0.05,]
nrow(DRN)
head(DRN)

# Differentially Regulated Edges (DREs)
DRE <- sem2$dest[sem2$dest$pvalue < 0.05,]
nrow(DRE)
head(DRE)


## 3.2. Total effect estimation.

# Average Causal Effect (ACE) with Benjamini-Hochberg adjustment (BH)

ace <- SEMace(graph = alsData$graph, data = data.npn, method = "BH")

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


## 3.3. Gene set analysis (GSA).

# Remove too small or too large pathways

n <- unlist(lapply(1:length(kegg.pathways),
            function(x) vcount(kegg.pathways[[x]])))
blacklist <- which(n < 5 | n > 500)
length(blacklist)
pathways <- kegg.pathways[-blacklist]

# GSA with Benjamini-Hochberg adjustment (BH)

GSA <- SEMgsa(pathways, data.npn, alsData$group, method = "BH", alpha = 0.05)



### Section 4. Causal structure learning. ###


## 4.1. Directed Acyclic Graph (DAG) estimation.
DAG <- SEMdag(graph = alsData$graph, data = data.npn, gnet = kegg,
              d = 2,
              beta = 0,
              lambdas = NA,
              verbose = FALSE)


## 4.2. Bow-free Acyclic Path (BAP) search.

BAP <- SEMbap(graph = alsData$graph, data = data.npn,
              method = "BH", 
              alpha = 0.05,
              limit = 30000,
              verbose = FALSE)


## 4.3. Graph extension.

ext <- extendGraph(g = list(DAG$dag, DAG$dag.red), data = data.npn,
                   gnet = kegg,
                   verbose = FALSE)


## 4.4. Model estimation strategies.

# Basic strategy (manuscript example).

# - Data-driven model estimation (fixed gnet = NULL and d = 0).
# - Decrease beta to improve model fitting.
model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = NULL,
                     d = 0,
                     search = "basic",
                     beta = 0.1, 
                     alpha = 0.05,
                     pstop = TRUE,
                     verbose = FALSE)


## Figure 3. ALS improved model (basic strategy). --------------------##

V(model$graph)$label <- mapIds(org.Hs.eg.db, V(sem1$graph)$name,
                               column = 'SYMBOL',
                               keytype = 'ENTREZID')

pert <- SEMrun(model$graph, model$data, alsData$group)

path <- SEMpath(model$graph, model$data, alsData$group,
                from = "6647",
                to = "4741",
                path = "directed",
                verbose = TRUE)

pdf("Figure3.pdf", width = 14, height = 19)
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
# - A smaller starting beta value is suggested (beta = 0.05).
model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = kegg,
                     d = 1,
                     search = "direct",
                     beta = 0.05,
                     alpha = 0.05,
                     pstop = TRUE,
                     verbose = FALSE)

# Outer strategy.
# - Knowledge-based estimation (gnet should be a directed reference network).
# - Up to d - 1 mediators can be imported from the reference network.
# - Larger d values increase model complexity (suggested: d = 2).
model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = kegg,
                     d = 2,
                     search = "outer",
                     beta = 0.1, 
                     alpha = 0.05,
                     pstop = TRUE,
                     verbose = FALSE)

# Inner strategy.
# - The reference network is used to validate new interactions and mediators.
# - Inferred mediators must already belong to the input graph.
# - Larger d values increase model complexity (suggested: d = 2).
model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = kegg,
                     d = 2,
                     search = "inner",
                     beta = 0.1, 
                     alpha = 0.05,
                     pstop = TRUE,
                     verbose = FALSE)



### Section 5. Network clustering and scoring. ###


# Improved ALS model clustering and scoring, using a latent variable 
# "hidden" model (LV), edge betweeness clustering (EBC) algorithm, 
# and a minimum cluster size of 5 nodes.

U <- clusterScore(model$graph, model$data, alsData$group,
                  HM = "LV",
                  type = "ebc",
                  size = 5)

scores <- parameterEstimates(U$fit)
head(scores)

# Clustering only (no scores calculation)

C <- clusterGraph(model$graph, type = "wtc", size = 5, verbose = FALSE)

# Cluster plot utility

G <- cplot(graph = model$graph, membership = U$membership, map = TRUE)

# Cluster extraction, fitting, and perturbation evaluation

G <- extractClusters(model$graph, model$data, alsData$group,
                     membership = U$membership)


## Figure 4. Edge betweenness clusters (EBC) mapped over the improved model. --##

# Graph plot
pdf("E1_ALS_clusters.pdf", width = 16, height = 8)
gplot(G$graph, fontsize = 30)
dev.off()

##-----------------------------------------------------------------------------##


### Section 6. Network weighting and filtering. ###


# Network weighting
# weighting method: r2z;
# seed thresholds (alpha, h, q): c(0.05, 0.5, 0.5).

W <- weightGraph(alsData$graph, alsData$exprs, alsData$group,
                 method = "r2z",
                 seed = c(0.05, 0.5, 0.5))

# Active module(s) finding (using the random walk with restart algorithm)

R <- activeModule(W, type = "rwr", seed = "pvlm", eweight = "pvalue")



### Section 7. Graph conversion utility. ###


# SEM (lavaan syntax) from an igraph object

als.sem <- graph2lavaan(model$graph)

# igraph network object from SEM (lavaan syntax)

als.graph <- lavaan2graph(als.sem, directed = TRUE, psi = TRUE)

# Extract a DAG from a network (igraph format)

DAG <- graph2dag(model$graph, model$data)

# Extract an undirected network from a correlation matrix R

R <- cor(model$data)
U <- corr2graph(R, n = nrow(model$data), type = "marg",
                method = "BH", 
                aplha = 0.05)



### Section 8. Disease module detection. ###


# FTD-related pathway selection

ftd.pathways <- c("MAPK signaling pathway",
                  "Protein processing in endoplasmic reticulum",
                  "Endocytosis",
                  "Wnt signaling pathway",
                  "Notch signaling pathway",
                  "Neurotrophin signaling pathway")
j <- which(names(kegg.pathways) %in% ftd.pathways)

# Nonparanormal transform of DNAme PC1 data

pc1.npn <- huge.npn(ftdDNAme$pc1)

# Gene set analysis (GSA)

ftd.gsa <- SEMgsa(kegg.pathways[j], pc1.npn, ftdDNAme$group, n_rep = 5000)

# Input graph as the union of FTD KEGG pathways

graph <- graph.union(kegg.pathways[j])
graph <- properties(graph)[[1]]

# Seed extraction

seed <- V(graph)$name[V(graph)$name %in% unique(unlist(ftd.gsa$DRN))]

# Graph weighting

W <- weightGraph(graph, pc1.npn, ftdDNAme$group, method = "r2z")

# Perturbed backbone as a Steiner tree (Kou's algorithm)

R <- activeModule(W, type = "kou", seed = seed, eweight = "pvalue")

# Entrez ID conversion

V(R)$label <- mapIds(org.Hs.eg.db, V(R)$name, 'SYMBOL', 'ENTREZID')

# Perturbation evaluation and plotting

pert <- SEMrun(graph = R, data = pc1.npn, group = ftdDNAme$group)


## Figure 5. FTD perturbed backbone. ---------------------------------##

pdf("Figure5.pdf", width = 16, height = 12)
gplot(pert$graph, fontsize = 30)
dev.off()

##--------------------------------------------------------------------##


# Backbone improvement (not shown in the manuscript).

model <- modelSearch(graph = R, data = data.npn, gnet = NULL,
                     d = 0,
                     search = "basic",
                     beta = 0.05,
                     alpha = 0.05,
                     pstop = TRUE,
                     verbose = FALSE)

V(model$graph)$label <- mapIds(org.Hs.eg.db, V(model$graph)$name, 'SYMBOL', 'ENTREZID')

pert <- SEMrun(graph = model$graph, data = pc1.npn, group = ftdDNAme$group)
gplot(pert$graph)



