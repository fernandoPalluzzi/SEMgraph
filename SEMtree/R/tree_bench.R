#load libraries
rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
#setwd("SEMtree_code_data")
library(GEOquery)
library(SEMgraph)
library(org.Hs.eg.db)
library(COSINE)
library(BioNet)
library(pathfindR)
library(diffusr)
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

# Differential expression analysis (deAna)
se <- GSE$se; dim(se) #13989 69
group <- se@colData@listData[["GROUP"]]; table(group) #0 23; 1 46
data <- t(assay(GSE$se))

# Check DEGs
degs <- rownames(assay(se))[which(rowData(se)$ADJ.PVAL < 5E-6)]
length(degs) #1003

# Interactome 
network <- as.undirected(kegg) 
genes <-rownames(assay(se))[rownames(assay(se)) %in% V(network)$name]
network <- induced_subgraph(network,genes)
network <- properties(network)[[1]]
summary(network)
#IGRAPH 17381bb UNW- 3032 19733 -- 
#+ attr: name (v/c), weight (e/n)


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


## (1) ALL methods

# 1) Bionet
pval <- rowData(se)$PVAL 
names(pval) <- rownames(se)
fb <- fitBumModel(pval, plot=TRUE)

scores <- scoreNodes(network=network, fb=fb, fdr=5e-6)
table(ifelse(scores < 0, 0, 1)) #1=396
mod1 <- runFastHeinz(network, scores)
Vids <- V(mod1)$name[V(mod1)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod1, data=data, h=0.1) #minimax clustering
tree <- extract_tree(mod1, data, group, mmx, h=0.1, alpha=0.025, Vids)
BioNet<-tree; dev.off(); gplot(BioNet)
#save(BioNet, file="./graph/BioNet.graph")


# 2) COSINE
PPI<- igraph::as_data_frame(network)[,-3]
set<- t(assay(se)[genes,])
control_data <- set[1:23,]; case_data <- set[24:69,]

scaled_diff <- diff_gen_PPI(case_data,control_data,PPI)
scaled_node_score <- scaled_diff[[1]]
scaled_edge_score <- scaled_diff[[2]]

start<- Sys.time()
lambda<- 0.5; num_iter<- 1000
GA_result<- GA_search_PPI(lambda,scaled_node_score,scaled_edge_score, PPI,
                          num_iter, muCh=0.05, zToR=10, minsize=10)
end<- Sys.time()
print(end-start)
#Time difference of 20.14278 mins
#save(GA_result, file="./data/GA_result.Rdata")
# load("GA_result.Rdata")

seedID <- GA_result$Subnet
seed<- V(network)$name[V(network) %in% seedID] 
mod2 <- induced_subgraph(network, V(network)$name[V(network)$name %in% seed])
Vids <- V(mod2)$name[V(mod2)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod2, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod2, data, group, mmx, h=0.2, alpha=0.025, Vids)
COSINE<-tree; gplot(COSINE)
#save(COSINE, file="./graph/COSINE.graph")

# 3) pathfindeR

RA_input <- rowData(se)[, c(1,2,5)]
Symbol<- unlist(mget(RA_input[,1], org.Hs.egSYMBOL, ifnotfound=NA))
RA_input[,1] <- Symbol
RA_input <- data.frame(RA_input)
colnames(RA_input) <- c("GENE", "CHANGE", "P_VALUE")
head(RA_input)

RA_processed <- pathfindR::input_processing(input = RA_input, 
                                            p_val_threshold = 5e-6, 
                                            #pin_name_path = "KEGG_true.sif",
                                            pin_name_path = "KEGG",
                                            convert2alias = TRUE)
#Number of genes provided in input: 13989
#Number of genes in input after p-value filtering: 1003
#Could not find any interactions for 615 (61.32%) genes in the PIN
#Final number of genes in input: 387

# greedy search (GR)
n_iter <- 10 
snws_list <- list()
for (i in 1:n_iter){
  cat("Iteration", i)
  snws_file <- "active_snws_GR"
  snws_tot <- pathfindR::active_snw_search(
    input_for_search = RA_processed,
    #pin_name_path = "KEGG.sif",
    pin_name_path = "KEGG",
    snws_file = snws_file,
    score_quan_thr = 0.8, 
    sig_gene_thr = 0.01, 
    search_method = "GR", 
    silent_option = TRUE)
  snws_list[[i]] <- snws_tot
}

#Found 87 active subnetworks

seed <- lapply(1:length(snws_list), function(x) snws_list[[x]])
seed <- unique(unlist(seed))
seed <- mapIds(org.Hs.eg.db, seed, column = 'ENTREZID', keytype = 'SYMBOL')
seed <- as.character(na.omit(seed))
mod3 <- induced_subgraph(network, V(network)$name[V(network)$name %in% seed])
mod3 <- properties(mod3)[[1]]
Vids <- V(mod3)$name[V(mod3)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod3, data=data, h=0.1) #minimax clustering
tree <- extract_tree(mod3, data, group, mmx, h=0.1, alpha=0.025, Vids)
pathfindeR<-tree; gplot(pathfindeR)
#save(pathfindeR, file="./graph/pathfindeR.graph")


# 4) SEMtree (1-abs(cor))

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
mod4<- SEMtree(network, data, seed=term, type="ST", eweight = NULL)
seed <- V(mod4)$name
Vids <- V(mod4)$name[V(mod4)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod4, data=data, h=0.15) #minimax clustering
tree <- extract_tree(mod4, data, group, mmx, h=0.15, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/ST.graph")


# 5) ST with eweight(r2z)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=group, method="r2z")
mod4<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod4)$name
Vids <- V(mod4)$name[V(mod4)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod4, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod4, data, group, mmx, h=0.2, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STr2z.graph")


# 6) WalktrapGM 

WTC<- WalktrapGM(network, data, group)
ID <- WTC$scores$ID[1]
seed <- names(WTC$membership[WTC$membership == ID]) #length(seed) #166
mod6 <- induced_subgraph(network, V(network)$name[V(network)$name %in% seed])
Vids <- V(mod6)$name[V(mod6)$name %in% cov]

# tree extraction and visualization
#mmx <- mergeNodes(mod6, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod6, data, group, mmx=NULL, h=0, alpha=0.025, Vids)
WalktrapGM<-tree; gplot(WalktrapGM)
#save(WalktrapGM, file="./graph/WalktrapGM.graph")


# 7) WalktrapGMl (FC)

WTC<- WalktrapGMl(network, rowData(se))
ID <- WTC$scores$ID[1]
seed <- names(WTC$membership[WTC$membership == ID]) #length(seed) #166
mod7 <- induced_subgraph(network, V(network)$name[V(network)$name %in% seed])
Vids <- V(mod7)$name[V(mod7)$name %in% cov]

# tree extraction and visualization
#mmx <- mergeNodes(mod6, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod7, data, group, mmx=NULL, h=0, alpha=0.025, Vids)
WalktrapGMl<-tree; gplot(WalktrapGMl)
#save(WalktrapGMl, file="./graph/WalktrapGMl.graph")


# 8) Coronavirus disease - COVID-19

path <- which(names(kegg.pathways) == "Coronavirus disease - COVID-19")
graph<- as.undirected(kegg.pathways[[path]])
nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
mod8<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
seed <- V(mod8)$name
Vids <- V(mod8)$name[V(mod8)$name %in% cov]

# tree extraction and visualization
#mmx <- mergeNodes(mod6, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod8, data, group, mmx=NULL, h=0, alpha=0.025, Vids)
KEGGCovid19<-tree; gplot(KEGGCovid19)
#save(KEGGCovid19, file="./graph/KEGGCovid19.graph")


# 9) COVID-19 article - SEMtree (CAT)

load("graph/tree_base.graph")
seed <- V(tree)$name
Vids <- V(tree)$name[V(tree)$name %in% cov]

# tree extraction and visualization
#mmx <- mergeNodes(mod6, data=data, h=0.2) #minimax clustering
tree <- extract_tree(tree, data, group, mmx=NULL, h=0, alpha=0.025, Vids)
COVID19_base<-tree; gplot(COVID19_base)
#save(KEGGCovid19, file="./graph/KEGGCovid19.graph")


## (2) ST methods

# 1) SEMtree (1-abs(cor))

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
mod1<- SEMtree(network, data, seed=term, type="ST", eweight = NULL)
seed <- V(mod1)$name
Vids <- V(mod1)$name[V(mod1)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod1, data=data, h=0.16) #minimax clustering
tree <- extract_tree(mod1, data, group, mmx, h=0.16, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/ST.graph")


# 2) SEMtree (one)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
E(network)$weight <- 1
mod2<- SEMtree(network, data, seed=term, type="ST", eweight = "custom")
seed <- V(mod2)$name
Vids <- V(mod2)$name[V(mod2)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod2, data=data, h=0.16) #minimax clustering
tree <- extract_tree(mod2, data, group, mmx, h=0.16, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STone.graph")


# 3) ST with eweight(r2zN)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=NULL, method="r2z")
mod3<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod3)$name
Vids <- V(mod3)$name[V(mod3)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod3, data=data, h=0.16) #minimax clustering
tree <- extract_tree(mod3, data, group, mmx, h=0.16, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STr2zN.graph")


# 4) ST with eweight(r2z)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=group, method="r2z")
mod4<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod4)$name
Vids <- V(mod4)$name[V(mod4)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod4, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod4, data, group, mmx, h=0.2, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STr2z.graph")


# 5) ST with eweight(sem)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=group, method="sem")
mod5<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod5)$name
Vids <- V(mod5)$name[V(mod5)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod5, data=data, h=0.16) #minimax clustering
tree <- extract_tree(mod5, data, group, mmx, h=0.16, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STsem.graph")

# 6) ST with eweight(cov)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=group, method="cov")
mod6<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod6)$name
Vids <- V(mod6)$name[V(mod6)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod6, data=data, h=0.2) #minimax clustering
tree <- extract_tree(mod6, data, group, mmx, h=0.2, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STcov.graph")

# 7) ST with eweight(cfa)

term <- rowData(se)$ENTREZID[rowData(se)$ADJ.PVAL < 5e-6]
graph <- weightGraph(graph=network, data=data, group=group, method="cfa")
mod7<- SEMtree(graph, data, seed=term, type="ST", eweight = "pvalue")
seed <- V(mod7)$name
Vids <- V(mod7)$name[V(mod7)$name %in% cov]

# tree extraction and visualization
mmx <- mergeNodes(mod7, data=data, h=0.16) #minimax clustering
tree <- extract_tree(mod7, data, group, mmx, h=0.16, alpha=0.025, Vids)
ST<-tree; gplot(ST)
#save(ST, file="./graph/STcfa.graph")


## END !!
