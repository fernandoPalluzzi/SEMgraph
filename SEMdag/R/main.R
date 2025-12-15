rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
library(SEMgraph)
library(SEMdata)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(limma)
library(huge)
library(stringr)
library(pcalg)
library(gnlearn)
library(CAM)
library(CMA)
library(caret)
library(readxl)
library(tidyr)
library(tibble)
library(dplyr)
library(ggpubr)

options (future.globals.maxSize = 4000 * 1024^5)
options(warn = -1)

source("Help.R")
source("DLiNGAM_p.R")
Rcpp::sourceCpp('direct_lingam_funcs.cpp')


### 0) Load data ---------------------------------------------------------------

# ALS data

path <- which(names(kegg.pathways) == "Amyotrophic lateral sclerosis")
ig <- properties(kegg.pathways[[path]])[[1]]; gplot(ig)
n <- vcount(ig); n

# train: GSE124439
df_train_lingam<- alsData$exprs
df_train_all<- huge::huge.npn(df_train_lingam)
group_train<- alsData$group; table(group_train)
ig <- graph2dag(ig, df_train_all)
reactome<- upgrade_graph(reactome);reactome
#REACTOME<- weightGraph(reactome, df_train_all, group = NULL, method = "r2z")

# test: GSE153960
data <- loadRData("data/als/GSE153960_SE.RData")
df_test_all <- t(assay(data[["se"]]))
df_test_all <- huge::huge.npn(df_test_all)
group_test <- data[["group"]]; table(group_test)

idx<- colnames(df_train_all)[colnames(df_train_all) %in% colnames(df_test_all)]
df_train_all <- df_train_all[,idx]
df_train_lingam <- df_train_lingam[,idx]
df_test_all <- df_test_all[,idx]
# head(colnames(df_train_all)); head(colnames(df_test_all))

design <- cbind(Intercept=1,Group=group_train)
fit <- lmFit(t(df_train_all), design)
fit <- eBayes(fit, trend=FALSE)
summary(decideTests(fit[,-1]))
top<- topTable(fit, coef=2, number=Inf, genelist=NULL, adjust="BH")
top05<- subset(top, adj.P.Val < 0.05)
seed<- rownames(top05[1:100,]) #top100 ranked genes with p(BH)<0.05

df_train <- df_train_all[, seed]
df_train_lingam<- df_train_lingam[, seed]
df_test<- df_test_all[, seed]
namesToDrop <- caret::findCorrelation(cor(df_train), cutoff = 0.99, names=TRUE)
if (length(namesToDrop) > 0) {
  df_train<- df_train[,-which(colnames(df_train) %in% namesToDrop)]
  df_train_lingam <- df_train_lingam[,-which(colnames(df_train_lingam) %in% namesToDrop)]
  df_test<- df_test[,-which(colnames(df_test) %in% namesToDrop)]
}
summary(decideTests(fit[,-1]))

namesToDrop #0

dim(df_train_all) #160 15936
dim(df_train) #160 100
dim(df_train_lingam) #160 100

dim(df_test_all) #273 15936
dim(df_test) #273 100


# BRCA data 

path <- which(names(kegg.pathways) == "Breast cancer")
ig <- properties(kegg.pathways[[path]])[[1]]; gplot(ig)
n <- vcount(ig); n

# train: TCGA BRCA
data<- loadRData("data/brca/brcatrain_SE.RData")
df_train_lingam <- t(assay(data[["se"]]))
df_train_all <- huge::huge.npn(df_train_lingam); dim(df_train_all)
group_train <- data[["group"]]; table(group_train)
ig <- graph2dag(ig, df_train_all); summary(ig)
reactome<- upgrade_graph(reactome);reactome
#REACTOME<- weightGraph(reactome, df_train_all, group = NULL, method = "r2z")

# test: GSE81538 + GSE205725
data<- loadRData("data/brca/brcatest_SE.RData")
df_test_all <- t(assay(data[["se"]])); dim(df_test_all)
df_test_all <- huge::huge.npn(df_test_all);
group_test <- data[["group"]]; table(group_test)

idx<- colnames(df_train_all)[colnames(df_train_all) %in% colnames(df_test_all)]
df_train_all <- df_train_all[,idx]
df_train_lingam <- df_train_lingam[,idx]
df_test_all <- df_test_all[,idx]
# head(colnames(df_train_all)); head(colnames(df_test_all))

design <- cbind(Intercept=1,Group=group_train)
fit <- lmFit(t(df_train_all), design)
fit <- eBayes(fit, trend=FALSE)
summary(decideTests(fit[,-1]))
top<- topTable(fit, coef=2, number=Inf, genelist=NULL, adjust="BH")
top05<- subset(top, adj.P.Val < 0.05)
seed<- rownames(top05[1:100,]) #top100 ranked genes with p(BH)<0.05

df_train <- df_train_all[, seed]
df_train_lingam<- df_train_lingam[, seed]
df_test<- df_test_all[, seed]
namesToDrop <- caret::findCorrelation(cor(df_train), cutoff = 0.99, names=TRUE)
if (length(namesToDrop) > 0) {
  df_train<- df_train[,-which(colnames(df_train) %in% namesToDrop)]
  df_train_lingam <- df_train_lingam[,-which(colnames(df_train_lingam) %in% namesToDrop)]
  df_test<- df_test[,-which(colnames(df_test) %in% namesToDrop)]
}
summary(decideTests(fit[,-1]))

namesToDrop #0

dim(df_train_all) #224 12940
dim(df_train) # 224 100
dim(df_train_lingam) # 224 100

dim(df_test_all) #377 12940 
dim(df_test) #377 100 


# COVID-19 data 

path <- which(names(kegg.pathways) == "Coronavirus disease - COVID-19")
ig <- properties(kegg.pathways[[path]])[[1]]; gplot(ig)
n <- vcount(ig); n

# train: GSE157103
data <- loadRData("data/covid19/GSE157103_SE.RData")
df_train_lingam <- t(assay(data[["se"]]))
df_train_all <- huge::huge.npn(df_train_lingam)
group_train <- data[["group"]]; table(group_train)
ig <- graph2dag(ig, df_train_all)
reactome<- upgrade_graph(reactome);reactome
#REACTOME<- weightGraph(reactome, df_train_all, group = NULL, method = "r2z")

# test: GSE152641
data <- loadRData("data/covid19/GSE152641_SE.RData")
df_test_all <- t(assay(data[["se"]]))
df_test_all <- huge::huge.npn(df_test_all)
group_test <- data[["group"]]; table(group_test)

idx<- colnames(df_train_all)[colnames(df_train_all) %in% colnames(df_test_all)]
df_train_all <- df_train_all[,idx]
df_train_lingam <- df_train_lingam[,idx]
df_test_all <- df_test_all[,idx]
# head(colnames(df_train_all)); head(colnames(df_test_all))

design <- cbind(Intercept=1,Group=group_train)
fit <- lmFit(t(df_train_all), design)
fit <- eBayes(fit, trend=FALSE)
summary(decideTests(fit[,-1]))
top<- topTable(fit, coef=2, number=Inf, genelist=NULL, adjust="BH")
top05<- subset(top, adj.P.Val < 0.05)
seed<- rownames(top05[1:100,]) #top100 ranked genes with p(BH)<0.05

df_train <- df_train_all[, seed]
df_train_lingam<- df_train_lingam[, seed]
df_test<- df_test_all[, seed]
namesToDrop <- caret::findCorrelation(cor(df_train), cutoff = 0.99, names=TRUE)
if (length(namesToDrop) > 0) {
  df_train<- df_train[,-which(colnames(df_train) %in% namesToDrop)]
  df_train_lingam <- df_train_lingam[,-which(colnames(df_train_lingam) %in% namesToDrop)]
  df_test<- df_test[,-which(colnames(df_test) %in% namesToDrop)]
}
summary(decideTests(fit[,-1]))

namesToDrop #0

dim(df_train_all) #126 14812
dim(df_train) #126 100
dim(df_train_lingam) #126 100

dim(df_test_all) #86 14812
dim(df_test) #86 100


# STEMI data

path <- which(names(kegg.pathways) == "Lipid and atherosclerosis")
ig <- properties(kegg.pathways[[path]])[[1]]; gplot(ig)
n <- vcount(ig); n

# train: GSE59867 
# data<- loadRData("data/stemi/GSE62646_SE_new.RData")
data <- loadRData("data/stemi/GSE59867_SE_new.RData")
df_train_lingam <- t(assay(data[["se"]]))
df_train_all <- huge::huge.npn(df_train_lingam)
group_train <- data[["group"]]; table(group_train)
ig <- graph2dag(ig, df_train_all)
reactome<- upgrade_graph(reactome);reactome
#REACTOME<- weightGraph(reactome, df_train_all, group = NULL, method = "r2z")

# test: GSE62646
data<- loadRData("data/stemi/GSE62646_SE_new.RData")
# data <- loadRData("data/stemi/GSE59867_SE_new.RData")
df_test_all <- t(assay(data[["se"]]))
df_test_all <- huge::huge.npn(df_test_all)
group_test <- data[["group"]];  table(group_test)

idx<- colnames(df_train_all)[colnames(df_train_all) %in% colnames(df_test_all)]
df_train_all <- df_train_all[,idx]
df_train_lingam <- df_train_lingam[,idx]
df_test_all <- df_test_all[,idx]
# head(colnames(df_train_all)); head(colnames(df_test_all))

design <- cbind(Intercept=1,Group=group_train)
fit <- lmFit(t(df_train_all), design)
fit <- eBayes(fit, trend=FALSE)
summary(decideTests(fit[,-1]))
top<- topTable(fit, coef=2, number=Inf, genelist=NULL, adjust="BH")
top05<- subset(top, adj.P.Val < 0.05)
seed<- rownames(top05[1:100,]) #top100 ranked genes with p(BH)<0.05

df_train <- df_train_all[, seed]
df_train_lingam<- df_train_lingam[, seed]
df_test<- df_test_all[, seed]
namesToDrop <- caret::findCorrelation(cor(df_train), cutoff = 0.99, names=TRUE)
if (length(namesToDrop) > 0) {
  df_train<- df_train[,-which(colnames(df_train) %in% namesToDrop)]
  df_train_lingam <- df_train_lingam[,-which(colnames(df_train_lingam) %in% namesToDrop)]
  df_test<- df_test[,-which(colnames(df_test) %in% namesToDrop)]
}
summary(decideTests(fit[,-1]))

namesToDrop #0

dim(df_train_all) #157 19555
dim(df_train) #157  99
dim(df_train_lingam) #157  99

dim(df_test_all) #42 19555
dim(df_test) #42 99



### 1) Search DAGs and Resizing ------------------------------------------------

#load interactome & create empty graph
ig0 <- make_empty_graph(ncol(df_train))
V(ig0)$name<- colnames(df_train)
d<- floor(igraph::mean_distance(ig, weights=NA));d
# d<- round(igraph::mean_distance(ig, weights=NA));d
md<- round(mean(igraph::degree(ig,  mode = "out")));md

# EqVarDAG, KEGG topological node-based
dag <- SEMdag(graph=ig, data=df_train_all, LO="TO", beta=0.1, penalty=FALSE)
dagKEGG.1 <- dag$dag; is_dag(dagKEGG.1); summary(dagKEGG.1)

# EqVarDAG, KEGG topological layer-based 
dag <- SEMdag(graph=ig, data=df_train_all, LO="TL", beta=0.1, eta = NULL, penalty=FALSE)
dagKEGG.2 <- dag$dag; is_dag(dagKEGG.2); summary(dagKEGG.2)

# EqVarDAG, BU Residual Variance (glasso), node-based eta=0
dag <- SEMdag(graph=ig0, data=df_train, LO="BU", beta=0.1, eta=0, penalty=FALSE)
dag1.1 <- dag$dag; is_dag(dag1.1); summary(dag1.1)

# EqVarDAG, BU Residual Variance (glasso), layer-based eta=NULL
dag <- SEMdag(graph=ig0, data=df_train, LO="BU", beta=0.1, eta = NULL, penalty=FALSE)
dag1.2 <- dag$dag; is_dag(dag1.2); summary(dag1.2) 

# PC, Constraint-Based Algorithm
dag <- pcalg::pc(suffStat=list(C = cor(df_train), n = nrow(df_train)),
                 indepTest=gaussCItest, alpha=0.05, labels=colnames(df_train),
                 skel.method="stable.fast", m.max=Inf, verbose=FALSE)
dag <-  pcalg::pdag2dag(dag@graph)
dag2<- graph_from_graphnel(dag$graph); is_dag(dag2); summary(dag2) 

# GES (Greedy Equivalence Search), Score-based Algorithm
score <- new("GaussL0penObsScore", df_train) #define the score (BIC)
dag <- pcalg::ges(score, maxDegree = md)
dag <- pcalg::pdag2dag(as(dag$essgraph, "graphNEL"))
dag3<- graph_from_graphnel(dag$graph); is_dag(dag3); summary(dag3)

# ARGES (Adaptively Restricted GES), Hybrid Algorithm
score <- new("GaussL0penObsScore", df_train)
suffStat <- list(C= cor(df_train), n=nrow(df_train))
skel.fit <- pcalg::skeleton(suffStat = suffStat, indepTest = gaussCItest, 
                            alpha = 0.05, labels = colnames(df_train))
skel <- as(skel.fit@graph, "matrix")
dag <- pcalg::ges(score, fixedGaps= !skel, adaptive = "triples", maxDegree = md)
dag <- pcalg::pdag2dag(as(dag$essgraph, "graphNEL"))
dag4<- graph_from_graphnel(dag$graph); is_dag(dag4); summary(dag4)

# Max-Min Hill Climbing (MMHC), Hybrid Algorithm
dag <- bnlearn::rsmax2(as.data.frame(df_train))
dag5<- graph_from_edgelist(dag$arcs, directed = TRUE); is_dag(dag5)
summary(dag5)

# LINGAM (DirectLINGAM), Non-Gaussian SEM Algorithm
df_train_lingam<- scale(df_train_lingam); dim(df_train_lingam)
dag <- DirectLINGAM(as.matrix(df_train_lingam), beta = 0.05)$DAG #0.05
colnames(dag) <- rownames(dag) <- colnames(df_train_lingam)
dag6 <- graph_from_adjacency_matrix(dag, mode = "directed"); is_dag(dag6)
summary(dag6)

# CAM (Causal Additive Model)
start<- Sys.time()
dag <- CAM(
  X = df_train,
  parsScore = list(numBasisFcts = 10),
  maxNumParents = md,
  #output = TRUE, 
  variableSel = TRUE,
  variableSelMethod = selGamBoost
  #pruning = TRUE,
  #pruneMethod = selGam
)
end <- Sys.time()
print(end - start)
#Time difference of 1.776724 mins
dag7 <- graph_from_edgelist(dag$edgeList, directed = TRUE); is_dag(dag7)
V(dag7)$name <- colnames(df_train); summary(dag7)

# NOTEARS (gnlearn), Continuous Optimization Algorithm
start<- Sys.time()
dag <- notears(df=df_train, lambda1=0.1, loss.type="l2", max.iter=100,
               h.tol=1e-8, rho.max=1e+16, w.threshold=0.1)
dag8<- graph2dag(dag, df_train, bap=FALSE); is_dag(dag8) 
end <- Sys.time()
print(end - start)
summary(dag8)

gg <- list(dag1.1, dag1.2, dag2, dag3, dag4, dag6, dag7, dag8) #-mmhc
names(gg) <- c("EV_BU_TO","EV_BU_TL","PC","GES","ARGES","LiNGAM","CAM","NOTEARS")
gg_kegg <- list(dagKEGG.1, dagKEGG.2)
names(gg_kegg) <- c("EV_KEGG_TO", "EV_KEGG_TL")

dag <- c(gg_kegg, gg);dag

### 2) SEM-Based Out-of-Sample Predictions -------------------------------------

res_tot<- run_fda(gg=dag, df_train=df_train_all, group_train=group_train,
                  df_test=df_test_all, group_test=group_test) 
print(res_tot)

res_rf<- run_rf(df_train, group_train, df_test, group_test)
print(res_rf)

# END !

