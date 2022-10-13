#load libraries
rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
#setwd("SEMtree_code_data")
library(SEMgraph)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(COSINE)
library(BioNet)
library(pathfindR)
library(diffusr)
library(reshape)
library(ggplot2)
library(viridis)
library(dplyr)

source("./code/help_tree.R")
options(warn = -1)

# COSINE simulation 

network<- as.undirected(kegg)
seed_vec <- c(1:100)
set_vec <- c(1:4)
nset <- 4; nmet <- 7
results <- data.frame(iter = rep(1:100, each = nset*nmet), 
                      set = rep(rep(1:nset, each=nmet),100),
                      FG = rep(c(rep(50,3*nmet), rep(40, nmet)), 100),
                      method = rep(c("BioNet", "COSINE", "pathfindeR",
                                     "ST", "STr2z", "WalktrapGM", "WalktrapGMl"),400),
                      size = rep(NA,nset*nmet*100), precision = rep(NA,nset*nmet*100),
                      recall = rep(NA,nset*nmet*100), f1 = rep(NA,nset*nmet*100))
n = 1
for (j in seed_vec){
  
  set.seed(j)
  simulated_data <- generate_sim(mu = 0.75)
  
  for (set in set_vec){
    
    cat("\n Iteration", j, "- Set", set, "\n")
    
    if(set == 1) SIMdata <- t(rbind(simulated_data[[7]], simulated_data[[2]]))
    if(set == 2) SIMdata <- t(rbind(simulated_data[[7]], simulated_data[[3]]))
    if(set == 3) SIMdata <- t(rbind(simulated_data[[7]], simulated_data[[5]]))
    if(set == 4) SIMdata <- t(rbind(simulated_data[[7]], simulated_data[[6]]))
    
    colnames(SIMdata) <- c(rep("c",20), rep("d", 20))
    rownames(SIMdata) <- V(network)$name[1:nrow(SIMdata)]
    group <- c(rep(0,20), rep(1,20))
    
    if(set == 4) FG<- rownames(SIMdata)[1:40]
    if(set != 4) FG<- rownames(SIMdata)[1:50]
    
    y<- group; x<- SIMdata
    
    write.table(x, file="xdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    write.table(cbind(colnames(x),y), file="cdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    write.table(rownames(x), file="rdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    exprs.file <- file.path(getwd(), "xdat.tab")
    cdat.file <- file.path(getwd(), "cdat.tab")
    rdat.file <- file.path(getwd(), "rdat.tab")
    
    SIMdata <- EnrichmentBrowser::readSE(exprs.file, cdat.file, rdat.file)
    SIMdata <- deAna(SIMdata)
    data<- t(assay(SIMdata))
    # network <- subNetwork(rownames(SIMdata), gHumanNet)
    
    cat("Bionet...")
    
    pval <- rowData(SIMdata)$ADJ.PVAL
    names(pval) <- rownames(SIMdata)
    fb <- fitBumModel(pval, plot=FALSE)
    
    term <- rowData(SIMdata)$ENTREZID[rowData(SIMdata)$ADJ.PVAL < 0.05]
    scores<- ifelse(names(pval) %in% term, 1, -1)
    names(scores)<- names(pval)
    module <- runFastHeinz(network, scores)
    seed <- V(module)$name
    
    if (is.null(seed) | length(seed) == 1){
      
      results$size[n] <- 0
      results$precision[n] <- NA
      results$recall[n] <- NA
      results$f1[n] <- NA
      
    } else{
      
      results$size[n] <- length(seed)
      int <- length(intersect(seed, FG))
      pre <- int/length(seed) 
      rec <- int/length(FG) 
      f1 <- 2*((pre*rec)/(pre+rec))
      
      results$precision[n] <- round(pre,3)
      results$recall[n] <- round(rec,3)
      results$f1[n] <- round(f1,3)
      
    }
    
    cat("COSINE...")
    
    PPI<- igraph::as_data_frame(network)
    genes<- rownames(assay(SIMdata)) %in% V(network)$name
    SIMdatat<- t(assay(SIMdata)[genes,])
    
    control_data <- SIMdatat[1:20,]
    case_data <- SIMdatat[21:40,]
    
    scaled_diff <- quiet(diff_gen_PPI(case_data,control_data,PPI))
    scaled_node_score <- scaled_diff[[1]]
    scaled_edge_score <- scaled_diff[[2]]
    
    lambda<- 0.5
    num_iter<- 100 #300
    GA_result<- quiet(GA_search_PPI(lambda,scaled_node_score,scaled_edge_score, PPI,
                                    num_iter, muCh=0.05, zToR=10, minsize=10))
    
    results$size[n+1] <- GA_result$Subnet_size
    seedID <- GA_result$Subnet
    seed<- V(network)$name[V(network) %in% seedID]
    
    if (is.null(seed) | length(seed) ==1){
      
      results$size[n+1] <- 0
      results$precision[n+1] <- NA
      results$recall[n+1] <- NA
      results$f1[n+1] <- NA
      
    } else{
      
      results$size[n+1] <- length(seed)
      int<- length(intersect(seed, FG))
      pre  <- int/length(seed)
      rec <- int/length(FG)
      f1 <- 2*((pre*rec)/(pre+rec))
      
      results$precision[n+1] <- round(pre,3)
      results$recall[n+1] <- round(rec,3)
      results$f1[n+1] <- round(f1,3)
      
    }
    
    cat("pathfindR...")
    
    RA_input <- rowData(SIMdata)[, c(1,2,5)]
    Symbol<- unlist(mget(RA_input[,1], org.Hs.egSYMBOL, ifnotfound=NA))
    RA_input[,1] <- Symbol
    RA_input <- data.frame(RA_input)
    colnames(RA_input) <- c("GENE", "CHANGE", "P_VALUE")
    
    # greedy search (GR)
    
    snws_file <- paste0("active_snws_GR", n+2)
    snws_tot <- quiet(active_snw_search(input_for_search = RA_input,
                                        pin_name_path = "KEGG",
                                        snws_file = snws_file,
                                        score_quan_thr = 0.5, 
                                        sig_gene_thr = 0.01, 
                                        search_method = "GR", 
                                        silent_option = TRUE))
    
    if(length(snws_tot) == 0){
      
      results$size[n+2] <- 0
      results$precision[n+2] <- NA
      results$recall[n+2] <- NA
      results$f1[n+2] <- NA
      
    } else{
      
      if(length(snws_tot[[1]]) == 1){
        
        results$size[n+2] <- 0
        results$precision[n+2] <- NA
        results$recall[n+2] <- NA
        results$f1[n+2] <- NA
        
      } else{
        
        seed <- snws_tot[[1]]
        results$size[n+2] <- length(seed)
        seed <- quiet(mapIds(org.Hs.eg.db, seed, column = 'ENTREZID', keytype = 'SYMBOL'))
        seed<- as.character(na.omit(seed))
        
        int<- length(intersect(seed, FG))
        pre <- int/length(seed)
        rec <- int/length(FG)
        f1 <- 2*((pre*rec)/(pre+rec))
        
        results$precision[n+2] <- round(pre,3)
        results$recall[n+2] <- round(rec,3)
        results$f1[n+2] <- round(f1,3)
        
      }
    }
    
    cat("SEMtree...")
    
    data<- t(assay(SIMdata));group <- c(rep(0,20), rep(1,20))
    term <- rowData(SIMdata)$ENTREZID[rowData(SIMdata)$ADJ.PVAL < 0.05]
    
    if (identical(term, character(0)) | length(term) == 1){
      
      results$size[n+3] <- 0
      results$precision[n+3] <- NA
      results$recall[n+3] <- NA
      results$f1[n+3] <- NA
      
    } else{
      
      tree<- SEMtree(network, data, seed=term, type="ST")
      seed <- V(tree)$name
      
      results$size[n+3] <- length(seed)
      int <- length(intersect(seed, FG))
      pre <- int/length(seed)
      rec <- int/length(FG)
      f1 <- 2*((pre*rec)/(pre+rec))
      
      results$precision[n+3] <- round(pre,3)
      results$recall[n+3] <- round(rec,3)
      results$f1[n+3] <- round(f1,3)
      
    }
    
    
    data<- t(assay(SIMdata));group <- c(rep(0,20), rep(1,20))
    term <- rowData(SIMdata)$ENTREZID[rowData(SIMdata)$ADJ.PVAL < 0.05]
    
    if (identical(term, character(0)) | length(term) == 1){
      
      results$size[n+4] <- 0
      results$precision[n+4] <- NA
      results$recall[n+4] <- NA
      results$f1[n+4] <- NA
      
    } else{
      
      if (length(term) > 5){
        
        graph <- quiet(weightGraph(graph=network, data=data, group=group, method="r2z"))
        tree<- SEMtree(graph, data=data, seed=term, type="ST", eweight = "pvalue")
        seed <- V(tree)$name
        
        results$size[n+4] <- length(seed)
        int <- length(intersect(seed, FG))
        pre <- int/length(seed)
        rec <- int/length(FG)
        f1 <- 2*((pre*rec)/(pre+rec))
        
        results$precision[n+4] <- round(pre,3)
        results$recall[n+4] <- round(rec,3)
        results$f1[n+4] <- round(f1,3)
        
      } else{
        
        results$size[n+4] <- 0
        results$precision[n+4] <- NA
        results$recall[n+4] <- NA
        results$f1[n+4] <- NA
        
      }
    }
    
    cat("WalktrapGM...")
    
    data<- t(assay(SIMdata));group <- c(rep(0,20), rep(1,20))
    WTC<- quiet(WalktrapGM(network, data, group))
    ID <- WTC$scores$ID[1]
    seed <- names(WTC$membership[WTC$membership == ID]) #length(seed) #166
    
    if (is.null(seed) | length(seed) == 1){
      
      results$size[n+5] <- 0
      results$precision[n+5] <- NA
      results$recall[n+5] <- NA
      results$f1[n+5] <- NA
      
    } else{
      
      results$size[n+5] <- length(seed)
      int <- length(intersect(seed, FG))
      pre <- int/length(seed)
      rec <- int/length(FG)
      f1 <- 2*((pre*rec)/(pre+rec))
      
      results$precision[n+5] <- round(pre,3)
      results$recall[n+5] <- round(rec,3)
      results$f1[n+5] <- round(f1,3)
      
    }
    
    WTC<- WalktrapGMl(network, rowData(SIMdata))
    ID <- WTC$scores$ID[1]
    seed <- names(WTC$membership[WTC$membership == ID]) #length(seed) #166
    
    if (is.null(seed) | length(seed) == 1){
      
      results$size[n+6] <- 0
      results$precision[n+6] <- NA
      results$recall[n+6] <- NA
      results$f1[n+6] <- NA
      
    } else{
      
      results$size[n+6] <- length(seed)
      int <- length(intersect(seed, FG))
      pre <- int/length(seed)
      rec <- int/length(FG)
      f1 <- 2*((pre*rec)/(pre+rec))
      
      results$precision[n+6] <- round(pre,3)
      results$recall[n+6] <- round(rec,3)
      results$f1[n+6] <- round(f1,3)
      
    }
    
    n = n + 7
  }
  
}

# Summarise results 

res <- results %>% 
  group_by(method, set) %>% 
  summarise(across(c(size, precision, recall, f1), .f = list(mean = mean), na.rm = TRUE))

res <- as.data.frame(res)
res_long <- melt(res, id.vars = c("method", "set"))

res_long$method[res_long$method == "WalktrapGM"] <- "WGM_RWR"
res_long$method[res_long$method == "WalktrapGMl"] <- "WGM_FC"

cols <- c("recall_mean", "precision_mean", "f1_mean")
file=paste0('all_meas.tiff')
tiff(file=file,width=8, height=5, units="in", res=300)
setEPS()
postscript(width=8, height = 6, "fig4.eps")
res_long %>%
  dplyr::filter(variable %in% cols) %>%
  ggplot(aes(x=method, y=value, group=as.factor(variable), fill=as.factor(variable)))+
  geom_bar(position='dodge', stat='identity') +
  ggtitle("") +
  theme_minimal() +
  ylab("") +
  xlab("") +
  facet_wrap(~ set, ncol = 2, scales = "free_y") +
  theme(text = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.title = element_text(size=9)) +
  theme(legend.title=element_blank())+
  scale_fill_discrete(labels=c('Precision', 'Recall', 'F-measure')) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.3, "cm"))
dev.off()

cols <- c("size_mean")
file=paste0('size_meas.tiff')
tiff(file=file,width=8, height=5, units="in", res=300)
# setEPS()
# postscript(width=8, height = 5, "size.eps")
res_long %>%
  dplyr::filter(variable %in% cols) %>%
  ggplot(aes(x=method, y=value, group=as.factor(set), fill=as.factor(set)))+
  geom_bar(position='dodge', stat='identity') +
  ggtitle("") +
  theme_minimal() +
  ylab("") +
  xlab("") +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_discrete("Set") +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.3, "cm"))
dev.off()
