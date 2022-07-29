### SEMgraph: An R Package for Causal Network Analysis of High-Throughput 
### Data with Structural Equation Models

# Replication code #


### Installation ###


#The latest stable version can be installed from CRAN:

install.packages("SEMgraph")

#The latest development version can be installed from GitHub:

install.packages("devtools")
devtools::install_github("fernandoPalluzzi/SEMgraph")

# Do not forget to install the SEMdata package too!
#It contains useful high-throughput sequencing data,
#reference networks, and pathways for SEMgraph training:

devtools::install_github("fernandoPalluzzi/SEMdata")

# Finally, install also the suggested packages

install.packages("org.Hs.eg.db")
install.packages("huge")

options(warn = -1)

# Load libraries

library(SEMgraph)
library(SEMdata)


###  Section 2: ALS  ###

## 2.1. SEM fitting functions.

# ALS model fitting (sem0: common model, no groups)

summary(alsData$graph) # ALS input graph
dim(alsData$exprs) # ALS RNA -seq expression data
table(alsData$group) # { case = 1, control = 0} vector

# Nonparanormal transform
library(huge)
data.npn <- huge.npn(alsData$exprs)

sem0 <- SEMrun(graph = alsData$graph, data = data.npn)

est0 <- parameterEstimates(sem0$fit)
head(est0)

# Other possible output (not shown in the manuscript).
summary(sem0$fit)

# ALS model fitting (sem1: common model, group influence on nodes)

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


# ALS model fitting (sem2: two models, group influence on edges)

sem2 <- SEMrun(alsData$graph, data.npn, alsData$group, fit = 2)

# Other possible output (not shown in the manuscript).
summary(sem2$fit)
parameterEstimates(sem2$fit)


# Perturbed graph elements

# Differentially Regualted Nodes (DRNs)
DRN <- sem1$gest[sem1$gest$pvalue < 0.05,]
nrow(DRN)
head(DRN)

# Differentially Regulated Edges (DREs)
DRE <- sem2$dest[sem2$dest$pvalue < 0.05,]
nrow(DRE)
head(DRE)


# RICF fitting
ricf1 <- SEMrun(alsData$graph, data.npn, alsData$group, algo = "ricf")

# Other possible output (not shown in the manuscript).
summary(ricf1$fit)
parameterEstimates(ricf1$fit)
print(ricf1$gest)

# CGGM fitting
cggm2 <- SEMrun(alsData$graph, data.npn, alsData$group, fit = 2, algo = "cggm")

# Other possible output (not shown in the manuscript).
summary(cggm2$fit$Group_0)
summary(cggm2$fit$Group_1)
parameterEstimates(cggm2$fit)
print(cggm2$dest)


## 2.2. Total effect estimation.

# Average Causal Effect (ACE) with optimal adjustement set,
# type="optimal"; source->sink ACE, effect="source2sink";
# Benjamini-Hochberg multiple test adjustment, method="BH";
# and significance level, alpha=0.05; NO bootstrap, boot=NULL

ace <- SEMace(graph = alsData$graph, data = data.npn,
			  type = "optimal", effect = "source2sink",
			  method = "BH", alpha = 0.05, boot=NULL)

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


## 2.3 Model estimation strategies

# Model Search with searc= "basic", and beta= 0.1

model <- modelSearch(graph = alsData$graph, data = data.npn,
                     gnet = NULL, d = 0,
                     search = "basic",
                     beta = 0.1,
					 method = "BH",
                     alpha = 0.05,
                     verbose = TRUE)


## Figure 2. ALS improved model (basic strategy). --------------------##

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(model$graph)$label <- mapIds(org.Hs.eg.db, V(model$graph)$name,
                               column = 'SYMBOL',
                               keytype = 'ENTREZID')

pert <- SEMrun(model$graph, model$data, alsData$group)

path <- SEMpath(model$graph, model$data, alsData$group,
                from = "6647",
                to = "4741",
                path = "directed",
                verbose = TRUE)

pdf("Figure2.pdf", width = 14, height = 19)
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
model$graph<- delete_vertex_attr(model$graph, "label")
gplot(model$graph, main="basic")
gplot(model1$graph, main="direct")
gplot(model2$graph, main="inner")
gplot(model3$graph, main="outer")


# Step 1-3 of modelSearch() with search="outer" run the following functions:

# Step 1: Bow-free Acyclic Path (BAP) search.

BAP <- SEMbap(graph = alsData$graph, data = data.npn,
              method = "BH", 
              alpha = 0.05,
              limit = 30000,
              verbose = TRUE)
			  
# Step 2: Directed Acyclic Graph (DAG) estimation.

DAG <- SEMdag(graph = alsData$graph, data = data.npn,
              LO = "TO", beta = 0.05,
			  lambdas = NA,
			  penalty = TRUE,
              verbose = TRUE)

# Step 3: Graph re-size with external interactome.

ext <- resizeGraph(g = list(alsData$graph, DAG$dag.new),
                   gnet = kegg,
				   d = 2,
				   v = TRUE,
                   verbose = TRUE)


## 2.4 Communities and factor scores

# Improved ALS model clustering and scoring, using a latent variable 
# "hidden" model (LV), edge betweeness clustering (EBC) algorithm, 
# and a minimum cluster size of 5 nodes.

LV <- clusterScore(model$graph, model$data, alsData$group,
				  type = "ebc",
				  HM = "LV",
				  size = 5)
				  
table(LV$membership)

head(parameterEstimates(LV$fit))

# Clustering only (no scores calculation)

C <- clusterGraph(model$graph,
				  type = "ebc", HM = "LV",
				  size = 5, verbose = TRUE)

# Cluster plot utility

cg <- cplot(graph = model$graph, membership = LV$membership, verbose=TRUE)
list(cg)
gplot(cg$graph)

# Cluster extraction, fitting, and perturbation evaluation

cls <- extractClusters(graph=model$graph, data=model$data, group=alsData$group, 
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


### Section 3: FTD ###

# 3.1 Gene Set Analysis (GSA)

#load libraries

library(SEMgraph)
library(SEMdata)
library(huge)

#Nonparanormal transform of DNAmePC1data

pc1.npn <- huge.npn(ftdDNAme$pc1); dim(pc1.npn)

#Set group classification

group <- ftdDNAme$group; table(group)

#Black list with n < 5 | n > 500

n <- unlist(lapply(1:length(kegg.pathways),
			 function(x) vcount(kegg.pathways[[x]])))
blacklist <- which( n < 5 | n > 500)
length(blacklist) #1

# SEM-based GSA

pathways <- kegg.pathways[-blacklist]
GSA <- SEMgsa(pathways, pc1.npn, group, method = "BH", alpha = 0.05)

#TOP10pathways

GSA$gsa[1:10,c(1:3,7)]


# 3.2 Fitting active disease modules

# FTD-related pathway selection

ftd.pathways <- c("MAPK signaling pathway",
                  "Protein processing in endoplasmic reticulum",
                  "Endocytosis",
                  "Wnt signaling pathway",
                  "Notch signaling pathway",
                  "Neurotrophin signaling pathway")
j <- which(names(kegg.pathways) %in% ftd.pathways)

# Gene set analysis (GSA) with 5000 permutations

ftd.gsa <- SEMgsa(kegg.pathways[j], pc1.npn, group, n_rep = 5000)
print(ftd.gsa$gsa)

# Seed extraction

seed <- unique(unlist(ftd.gsa$DEG))
length(seed) #60

# KEGG interactome weighting

W <- weightGraph(kegg, pc1.npn, group, method = "r2z")
summary(W) #4242 34975

# Perturbed backbone as a Steiner tree (Kou's algorithm)

ST <- SEMtree(W, data= pc1.npn, seed=seed, type="ST", eweight="pvalue")
summary(ST) # 92 92

# Perturbation evaluation with raw data

sem1 <- SEMrun(ST, pc1.npn, group)

# Perturbation evaluation with deconfounding data

adj.pc1 <- SEMbap(ST, pc1.npn, method = "bonferroni", alpha = 5E-06)$data
adj.sem1<- SEMrun(ST, adj.pc1, group)

#tree agglomerative hierarchical clustering (TAHC)

C <- clusterGraph(ST, type="tahc")
cg <- cplot(ST, membership=C)
list(cg)


## Figure 4. FTD perturbed backbone. ---------------------------------##

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(cg2)$label <- mapIds(org.Hs.eg.db, V(cg2)$name, 'SYMBOL', 'ENTREZID')

# Set node color and size
cg2 <- cg$HM2
V(cg2)$color <- ifelse(V(cg2)$name %in% seed, "green", "white")
V(cg2)$size <- 2*degree(cg2, mode="total")

pdf("Figure4.pdf", width = 16, height = 12)
gplot(cg2, l = "neato", fontsize = 30)
dev.off()

##--------------------------------------------------------------------##


# 3.3 Locating differentially connected genes

# Input graph as the union of FTD KEGG pathways

gU <- graph.union(kegg.pathways[j])
gU <- properties(gU)[[1]]
summary(gU) #586 3572

# SEM-based differential causal inference (DCI) with EBC clustering

gD <- SEMdci(gU, pc1.npn, group, type = "ebc", method = "BH", alpha = 0.05)
summary(gD) #111 103

# Perturbation evaluation of the 3th component of gD

gC <- properties(gD)

gC3 <- gC[[3]]

sem1 <- SEMrun(gC3, pc1.npn, group, fit=1, algo="cggm")

sem2 <- SEMrun(sem1$graph, pc1.npn, group, fit=2, algo="cggm")

#save(gC3, file="gC3.graph")
#load("gC3.graph")
#gplot(gC3, l = "fdp")


## Figure 5. DCI perturbed backbone. ---------------------------------##

# Set node and edge colors
gC3<- sem2$graph

# Convert Entrez identifiers to gene symbols
library(org.Hs.eg.db)
V(gC3)$label <- mapIds(org.Hs.eg.db, V(gC3)$name, "SYMBOL", "ENTREZID")

pdf("Figure5.pdf", width = 16, height = 12)
gplot(gC3, l = "fdp", fontsize = 30)
dev.off()

##--------------------------------------------------------------------##


## END !!
