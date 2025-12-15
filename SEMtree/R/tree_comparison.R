#load libraries
rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
#setwd("SEMtree_code_data")
library(GEOquery)
library(SEMgraph)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(clusterProfiler)
library(scales)
library(readxl)
library(dplyr)

source("./code/help_tree.R")
options(warn = -1)

getwd()


# Get data 

gse <- "GSE172114"
gse_filename <- "./data/GSE172114_series_matrix.txt.gz"
GSE <- getdata(gse_filename, gse, data.type = "rseq")
se <- GSE$se
assay(se)[1:5,1:5]
boxplot(assay(se))

group <- se@colData@listData[["GROUP"]]
table(group) #0 23; 1 46
data <- t(assay(GSE$se))
dim(data) #69 13989


# I LOAD TREES

load("./graph/BioNet.graph")
load("./graph/COSINE.graph")
load("./graph/pathfindeR.graph")
load("./graph/ST.graph")
load("./graph/STr2z.graph")
load("./graph/WalktrapGM.graph")
load("./graph/WalktrapGMl.graph")
load("./graph/KEGGCovid19.graph")
load("./graph/tree_base.graph")

# load("./graph/ST.graph")
# load("./graph/STone.graph")
# load("./graph/STr2zN.graph")
# load("./graph/STr2z.graph")
# load("./graph/STsem.graph")
# load("./graph/STcov.graph")
# load("./graph/STcfa.graph")

# setting tree list
xx<- list(BioNet, COSINE, pathfindeR, ST, STr2z, WalktrapGM, WalktrapGMl, KEGGCovid19, tree)
names(xx)<- c("BioNet", "COSINE", "pathfindeR", "ST", "STr2z", "WalktrapGM", "WalktrapGMl","KEGGCovid19", "tree")

# xx<- list(ST, STone, STr2zN, STr2z, STsem, STcov, STcfa)
# names(xx)<- c("ST", "STone", "STr2zN", "STr2z", "STsem", "STcov", "STcfa")


# 1) SEMgraph perturbation analysis

T<- length(xx); names(xx)

# SEMgsa
gsa.res<- SEMgsa(g=xx, data, group)
str(gsa.res, max.level=1)
gsa.res$gsa

# SEMace
ace.res<- lapply(xx[1:T], function(x) SEMace(graph=x, data=data, group=group, type="optimal", effect="source2sink", method="none", alpha=1, boot=NULL))
pval.ace <- scientific_format(3)(as.numeric(lapply(1:T, function(x) length(ace.res[[x]]$pvalue)*min(ace.res[[x]]$pvalue))))
names(pval.ace) <- names(xx); pval.ace

ace.res2<- lapply(xx[1:T], function(x) SEMace(graph=x, data=data, group=group, type="optimal", effect="source2sink", method="BH", alpha=0.05, boot=1000))
sig.ace <- round(as.numeric(lapply(1:T, function(x) (nrow(ace.res2[[x]])/nrow(ace.res[[x]]))*100)))
names(sig.ace) <- names(xx); sig.ace

# SEMrun
sem.res<-lapply(xx[1:T], function(x) SEMrun(graph=x, data=data, group=group, algo="ricf"))
str(sem.res, max.level=1)
lapply(sem.res[1:T], function(x) table(V(x$graph)$color))


# Tree visualization
treew<- lapply(xx[1:T], function(x) weightGraph(x, data, group=NULL, method = "r2z"))
for(i in 1:T) V(treew[[i]])$size <- igraph::degree(treew[[i]], mode = "all")*8
for(i in 1:T) E(treew[[i]])$width <- round(-log10(E(treew[[i]])$pv),3)/9

for(i in 1:T) {
  pdf(paste0(names(treew)[i],".pdf"), width = 14, height = 12)
  # V(treew[[i]])$color <- V(sem.res[[i]]$graph)$color
  # main= names(treew)[i]
  gplot(treew[[i]], psize = 50)
  dev.off()
  Sys.sleep(5)
}

# ST perturbed paths
ace <- SEMace(treew[[4]], data, group, effect = "source2sink", method="bonferroni", alpha=0.05)
ace <- ace[order(ace$pvalue),]; ace
ace$sink <- quiet(mapIds(org.Hs.eg.db, ace$sink, column = 'SYMBOL', keytype = 'ENTREZID'))
ace$source <- quiet(mapIds(org.Hs.eg.db, ace$source, column = 'SYMBOL', keytype = 'ENTREZID'))
ace$pvalue

# SEMpath
from<- "55054"; to<- "1234"
dp<- SEMpath(treew[[4]], data, group, from, to, path="directed")

# gplot(dp$map, main="X -> Y path")
for(i in 1:T) V(dp$map)$size <- igraph::degree(dp$map, mode = "all")*5
for(i in 1:T) E(dp$map)$width <- round(-log10(E(dp$map)$pv),3)/6
V(dp$map)$color<- V(treew[[4]])$color
V(dp$map)$color<- V(sem.res[[4]]$graph)$color
pdf("ST.pdf", width = 14, height = 12)
gplot(dp$map)
dev.off()


# 2) PenalizedLDA

#install.packages("penalizedLDA")
library(penalizedLDA)
library(MASS)

T <- length(xx); T
X <- data; dim(X)
y <- group + 1; table(y)
pLDA <- NULL

for (t in 1:T) {  #t=1
  # Dataset selection
  Y<-as.matrix(X[,colnames(X) %in% V(xx[[t]])$name]) # dim(Y)
  
  # Perform cross-validation for lambda and K parameters
  l1<-c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1, 10, sqrt(ncol(Y)))
  cv.out<-PenalizedLDA.cv(Y, y, lambdas=l1, K=NULL)
  out1<-PenalizedLDA(Y, y, lambda=cv.out$bestlambda, K=cv.out$bestK)
  zero<-length(which(out1$discrim==0))
  genes<-ncol(Y)
  
  d1<-out1$xproj[,1]
  #ldahist(d1,y, type="histogram", sep=TRUE)
  ldahist(d1,y, nbins= 5, type="density", sep=TRUE)
  
  dmean<-mean(d1[y==2])-mean(d1[y==1])
  dpval<-ncol(Y)*t.test(d1[y==2],d1[y==1])$p.value
  rss<-glm(d1~y,gaussian(link = "identity"))$deviance
  #summary(glm(d1~y,gaussian(link = "identity")))
  
  if (dmean > 0) { 
    class<- ifelse(d1 < (mean(d1[y==2]) + mean(d1[y==1]))/2, 1, 2)
  }else{
    class<- ifelse(d1 > (mean(d1[y==2]) + mean(d1[y==1]))/2, 1, 2)
  }
  tt<-table(y, class)
  se<-tt[2,2]/(tt[2,2]+tt[2,1])
  sp<-tt[1,1]/(tt[1,1]+tt[1,2])
  ac<-(tt[1,1]+tt[2,2])/sum(tt)
  
  pLDA<-rbind(pLDA, c(genes, zero, rss, dmean, dpval, se, sp, ac))
}

# data.frame of pLDA results 
pLDA<-as.data.frame(pLDA)
rownames(pLDA)<-names(xx)
colnames(pLDA)<-c("No. genes", "No. zero", "RSS", "dmean", "dpval", "Se", "Sp", "Acc")
pLDA


# 3) Enrichment gene/GO analysis

T<- length(xx); names(xx)
lapply(xx[1:T], function(x) table(V(x)$color))

# Import literature genes related to COVID-19 
cov1 <- read_excel("./data/ijmsv19p0402s2.xlsx") %>%
  dplyr::select("Related genes of COVID-19") %>%
  na.omit() %>%
  dplyr::pull() 
length(cov1) #245
cov1 <- quiet(mapIds(org.Hs.eg.db, cov1, column = 'ENTREZID', keytype = 'SYMBOL'))
sum(is.na(cov1))
cov1 <- as.vector(cov1[!is.na(cov1)]); length(cov1) #177

cov2 <- read.delim("./data/C0206750_disease_gda_summary.tsv") %>%
  dplyr::mutate(Gene_id = as.character(Gene_id)) %>%
  dplyr::pull(Gene_id)
length(cov2) #33

cov3 <- unique(c(cov1, cov2)); length(cov3) #196
cov <- cov3[cov3 %in% colnames(data)]
length(cov) #92

gene<- lapply(xx[1:T], function(x) gene_metrics(x, cov))
print(gene)

# Extract go terms related to COVID-19 based on literature genes
goldGO <- enrichGO(gene = cov,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable = TRUE)
goldGO <- goldGO@result

#write.table(goldGO, file=GOres.txt")
#goldGO<- read.table("./data/GOres.txt")[,1:2]; head(goldGO)

GO<- lapply(xx[1:T], function(x) GO_metrics(x, goldGO))
print(GO)


# 4) Nodes and edges Jaccard(J) similarity

nodesets<- lapply(xx[1:length(xx)], function(x) V(x)$name)
yy <- lapply(xx[1:length(xx)], function(x) as_edgelist(x))
#yy<- lapply(xx[1:length(xx)], function(x) get.data.frame(x))
edgesets <- lapply(yy[1:length(yy)], function(x) paste0(x[,1],"--",x[,2]))

J_nodes <- jaccard(nodesets); J_nodes
J_edges <- jaccard(edgesets); J_edges


# END !
