rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
setwd("SEMgsa_code_data")
library(GEOquery)
library(SEMgraph)
library(SEMdata)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(mirIntegrator)
library(qpgraph)
library(netgsa)
library(DEGraph)
library(ROntoTools)
library(graphite)
library(SPIA)
library(PathNet)
library(clipper)
library(a4Preproc)
library(a4Base)
library(qpgraph)
library(mosaic)
library(gplots)
library(ggpubr)
library(stringr)
library(data.table)
library(plyr)
library(dplyr)
source("./R/Help.R")
options(warn = -1)

# Get data 
gse="GSE172114"
gse_filename=paste0('./Data/', gse, '_series_matrix.txt.gz')
TRUEdata_list <- getGSAdata(gse_filename, gse = gse,
                            summary_probes = "random", log = TRUE)

# Check DEGs
se_true <- TRUEdata_list[["se"]]; group <- se_true@colData@listData[["GROUP"]]
de_true<- rownames(assay(se_true))[which(rowData(se_true)$ADJ.PVAL < 0.05)]; length(de_true)
exprsHeatmap(expr=assay(se_true)[de_true,], grp=group)

# Check intersection btw nodes of pathway and genes of data
geneSE = TRUEdata_list[["geneSE"]]; ig = TRUEdata_list[["graph"]]
length(V(ig)[V(ig) %in% rownames(geneSE)])

# Enrichment analysis (SIMULATION)
ncond = 2
disease = TRUEdata_list[["disease"]]
topology = "community"  # betweenness; neighborhood
muval_vec = seq(0.5,0.7,0.1); muval= muval_vec[1]
n = 10
N = 100

path_all <- getALLPATH(GSAdata_list = TRUEdata_list, disease_check = disease,
                       maxNodes=300, minNodes = 30, maxComp = 0.6)
COVIDall <- getCOVIDpath(path_all)

start<- Sys.time()
GSA_iter = list()
for (j in seq(1:N)){
  
  set.seed(j)
  cat("Iteration", j, "/", N, "- dereg:", topology, "- muval:", muval,"\n")
  
  e2a <- getDEGS(COVIDall, n=n, j=j, topology=topology)
  TP <- getTRUEpath(path_all, COVIDall, e2a$gene)
  SIMdata <- sim_design(TRUEdata_list, muval = muval, e2a = e2a, 
                        ncond = 2, seed = j)

  out <- getGSAresults(GSAdata_list = SIMdata, path_list = path_all)
  GSA_iter[[j]] <- getSUMMARY(out_list=out, TP=TP, sim = TRUE)
  
}
end<- Sys.time() 
print(end-start) 

# Simulation metrics
p <- GSAplot(GSA_iter, y="FP",type = "errorplot"); p
p <- GSAplot(GSA_iter, y="power",type = "errorplot"); p
p <- GSAplot(GSA_iter, y="prioritization",type = "errorplot"); p

