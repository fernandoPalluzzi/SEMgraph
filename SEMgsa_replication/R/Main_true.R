rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
setwd("SEMgsa_code_data")
library(GEOquery)
library(SEMgraph)
library(SEMdata)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(edgeR)
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
library(scales)
library(qpgraph)
library(huge)
library(mosaic)
library(gplots)
library(stringr)
library(data.table)
library(plyr)
library(dplyr)
source("./R/Help.R")
options(warn = -1)

# Get data 
gse_vec = c("GSE172114", "GSE53740"); gse=gse_vec[1]
gse_filename=paste0('./Data/', gse, '_series_matrix.txt.gz')
TRUEdata_list <- getGSAdata(gse_filename, gse = gse,
                            summary_probes = "random", log = TRUE)

# Check DEGs
se_true <- TRUEdata_list[["se"]]; group <- se_true@colData@listData[["GROUP"]]
de_true<- rownames(assay(se_true))[which(rowData(se_true)$ADJ.PVAL < 0.05)]; length(de_true)
exprsHeatmap(expr=assay(se_true)[de_true,], grp=group) 

# Enrichment analysis
disease = TRUEdata_list[["disease"]]
j = 1; set.seed(j)

# Get kegg.gr
kegg.gr <- kegg.pathways
graphnel <- lapply(1:length(kegg.gr), 
                   function(x) igraph::as_graphnel(kegg.gr[[x]]))
names(graphnel) <- names(kegg.gr)

# Get kegg.gs
kegg.gs <- lapply(1:length(kegg.gr), function(x)
  as_ids(V(kegg.gr[[x]])))
names(kegg.gs) <- names(kegg.gr)

path_all = list(kegg.gr = kegg.gr, graphnel = graphnel, kegg.gs = kegg.gs)

start<- Sys.time()
out = list()
out <- getGSAresults(GSAdata_list = TRUEdata_list, path_list = path_all)
end<- Sys.time()
print(end-start) 

if(gse == "GSE53740"){
  metrics_bind = c()
  for(name in names(TRUEdata_list[["graph"]])){
    GSA <- getSUMMARY(out_list=out, disease=name, sim = FALSE) 
    results <- GSA[["out_wide"]]
    metrics <- GSA[["metrics"]]
    metrics <- dplyr::mutate(metrics, disease = name)
    metrics_bind <- rbind(metrics_bind, metrics)
  }
  
} else{
  GSA <- getSUMMARY(out_list=out, disease=disease, sim = FALSE) 
  results <- GSA[["out_wide"]]
  metrics <- GSA[["metrics"]]
}


# Plot
ggplot(metrics, aes(x=factor(method), y=p.value, label = scientific_format(1)(metrics$p.value))) +
  geom_segment( aes(x=factor(method), xend=factor(method), y=0, yend=p.value), 
                color=c("#4DBBD5FF","#4DBBD5FF","#4DBBD5FF", "#4DBBD5FF","#E64B35FF","#4DBBD5FF")) +
  geom_point(color = c("#4DBBD5FF","#4DBBD5FF","#4DBBD5FF", "#4DBBD5FF","#E64B35FF","#4DBBD5FF"),size=5) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("P-value of target pathway") +
  geom_text(nudge_y = 0.03, size = 5) +
  theme(text = element_text(size = 18))

ggplot(metrics, aes(x=factor(method), y=round(rel.ranks), label = round(rel.ranks))) +
  geom_segment( aes(x=factor(method), xend=factor(method), y=0, yend=rel.ranks), 
                color=c("#4DBBD5FF","#4DBBD5FF","#4DBBD5FF", "#4DBBD5FF","#E64B35FF","#4DBBD5FF")) +
  geom_point(color = c("#4DBBD5FF","#4DBBD5FF","#4DBBD5FF", "#4DBBD5FF","#E64B35FF","#4DBBD5FF"),size=5) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Ranking of target pathway") +
  geom_text(nudge_y = 5, size = 5) +
  theme(text = element_text(size = 18)) 

