### MAIN REAL (SEMbap)

#load libraries
rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
#setwd("~/Desktop/SEMgraph/SEMbap")
#install.packages(c("RSpectra", "matrixcalc", "pcaPP", "cvTools"))
#devtools::install_github("benjaminfrot/lrpsadmm")
library(lrpsadmm)
library(mvtnorm)
library(dplyr)
library(data.table)
library(SEMgraph)
library(SEMdata)
library(org.Hs.eg.db)

source("help.R")
options(warn = -1)

getwd()

j <- 123
set.seed(j)

# raw BRCA data
load("brca.RData")
str(brca_wt)
str(brca_mt)
data <- t(cbind(brca_wt, brca_mt))
dim(data) #224 20501
data <- na.omit(data)
dim(data) #224 20501
data <- huge::huge.npn(data)
sd.res <- apply(data, 2, sd)
#hist(sd.res)
data <- data[, sd.res > 0.5]
dim(data) # 224 19247
colnames(data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(data), 'ENTREZID', 'SYMBOL'))
data <- data[, !is.na(colnames(data))]
dim(data) # 224 16277

pairwiseMatrix(x = scale(data, center = TRUE, scale = FALSE))

group <- c(rep(0, 112), rep(1, 112))

i <- which(names(kegg.pathways) == "Breast cancer");i
kegg.pathways[[i]]<- upgrade_graph(kegg.pathways[[i]])
properties(kegg.pathways[[i]])

graph <- properties(kegg.pathways[[i]])[[1]]
V(graph)$label <-
  mapIds(org.Hs.eg.db, V(graph)$name, 'SYMBOL', 'ENTREZID')
summary(graph) #IGRAPH b304e0f DNW- 133 483 --


#raw SEM data
sem0 <- SEMrun(graph, data, algo = "ricf", n_rep = 0)
g0 <- sem0$graph; g0 #IGRAPH d0ea8c9 DNW- 131 481 --
gplot(sem0$graph, l = "dot")
data0 <- sem0$data[, -1]
dim(data0) #224 131

properties(g0)
df0<- vcount(g0)*(vcount(g0)-1)/2-ecount(as.undirected(g0));df0 #8034
s0<- round(sqrt(nrow(data0))/log(ncol(data0)));s0 #3


# Transcription factors (ESR2 - FZD4; ESR2 - PIK3R2)

tfs <-
  read.table("trrust_rawdata.human.txt",
             header = TRUE,
             sep = "\t")
str(tfs)
colnames(data0) <-
  as.character(mapIds(org.Hs.eg.db, colnames(data0), 'SYMBOL', 'ENTREZID'))
tfs <-
  unique(tfs$TF)[unique(tfs$TF) %in% colnames(data0)]
length(tfs) #22
cor <- as.data.frame(as.table(cor(data0)))
cor <- subset(cor,
              abs(Freq) > 0.5 & abs(Freq) != 1 &
                cor$Var1 %in% tfs &
                !cor$Var2 %in% tfs)
dim(cor) #152 3

htf <-
  sort(table(cor$Var1), decreasing = TRUE)
# htf <- names(htf[1])
htf <- names(htf[2]); htf
genetf <- subset(cor, cor$Var1 %in% htf)
varmax <- genetf %>% slice_max(Freq, n = 1) %>%
  pull(Var2) %>% as.character()
varmin <- genetf %>% slice_min(Freq, n = 1) %>%
  pull(Var2) %>% as.character()
dataplot <-
  as.data.frame(data0[, colnames(data0) %in% c(htf, varmax, varmin)])

data1 <-
  data0[, colnames(data0)[!colnames(data0) %in% tfs]]
dim(data1) #224 108
colnames(data1) <-
  as.character(mapIds(org.Hs.eg.db, colnames(data1), 'ENTREZID', 'SYMBOL'))

sem1<- SEMrun(g0, data1, algo="ricf", n_rep=0)
g1 <- sem1$graph; g1 #IGRAPH eddc180 DNW- 109 397 -- 
#gplot(g1, l = "fdp")
data1 <- sem1$data[, -1]
dim(data1) #224 109

properties(g1)
df1<- vcount(g1)*(vcount(g1)-1)/2-ecount(as.undirected(g1));df1 #5489
s1<- round(sqrt(nrow(data1))/log(ncol(data1)));s1 #3

# plot g0+tfs

semplot <- SEMrun(graph, data, group=group, algo = "ricf")
graph <- semplot$graph
tfs <- as.character(mapIds(org.Hs.eg.db, tfs, 'ENTREZID', 'SYMBOL'))
V(graph)$color[which(V(graph)$name %in% tfs)] <- "yellow"

setEPS()
file <- "original.eps"
postscript(width=18, height = 12, file)
gplot(graph, l = "dot") #main = names(Y)[i]
dev.off()


##### SEMbap

#...with g1 & data1 
sem0 <- SEMrun(g1, data1, group=NULL, algo="ricf", n_rep=0)
sem1 <- SEMrun(g1, data1, group=group, algo="ricf", n_rep=1000)

s1<- round(sqrt(nrow(data1))/log(ncol(data1)));s1 #3
bap1<- SEMbap(g1, data1, dalgo="cggm", method="BH", alpha=5E-2, cmax=4, verbose=TRUE)

#if (dalgo == "glpc") beta <- ifelse(cls > 3, 1, 0.75)
#if (dalgo == "glpc") hcount= cls =length(cluster_leading_eigen(guu))
bap2<- SEMbap(g1, data1, dalgo="glpc", method="BH", alpha=5E-2, cmax=4, verbose=TRUE)

#bap3a<- SEMbap(g1, data1, dalgo="pc", hcount="auto")
bap3<- SEMbap(g1, data1, dalgo="pc", hcount=3, verbose=TRUE)

#bap4a<- SEMbap(g1, data1, dalgo="pcss", hcount="auto")
bap4<- SEMbap(g1, data1, dalgo="pcss", hcount=3, verbose=TRUE)

bap5<- SEMbap(g1, data1, dalgo="trim", hcount=0, verbose=TRUE)
bap6<- SEMbap(g1, data1, dalgo="lava", hcount=0, verbose=TRUE)
bap7<- SEMbap(g1, data1, dalgo="rsvp", hcount=0, verbose=TRUE)

bap8<- SEMbap(g1, data1, limit=1, verbose=TRUE) #glasso


#...with g0 & data0 
sem0 <- SEMrun(g0, data0, group=NULL, algo="ricf", n_rep=0)
sem1 <- SEMrun(g0, data0, group=group, algo="ricf", n_rep=1000)

s0<- round(sqrt(nrow(data0))/log(ncol(data0)));s0 #3
bap1<- SEMbap(g0, data0, dalgo="cggm", method="BH", alpha=5E-2, cmax=4, verbose=TRUE)

#if (dalgo == "glpc") beta <- ifelse(cls > 3, 1, 0.75)
#if (dalgo == "glpc") hcount= cls =length(cluster_leading_eigen(guu))
bap2<- SEMbap(g0, data0, dalgo="glpc", method="BH", alpha=5E-2, cmax=4, verbose=TRUE)

#bap3a<- SEMbap(g0, data0, dalgo="pc", hcount="auto")
bap3<- SEMbap(g0, data0, dalgo="pc", hcount=3, verbose=TRUE)

#bap4a<- SEMbap(g0, data0, dalgo="pcss", hcount="auto")
bap4<- SEMbap(g0, data0, dalgo="pcss", hcount=3, verbose=TRUE)

bap5<- SEMbap(g0, data0, dalgo="trim", verbose=TRUE)
bap6<- SEMbap(g0, data0, dalgo="lava", verbose=TRUE)
bap7<- SEMbap(g0, data0, dalgo="rsvp", verbose=TRUE)

bap8<- SEMbap(g0, data0, limit=1, verbose=TRUE) #glasso


#...lrpsadmml (data1)

lambda <- sqrt(log(ncol(data1)) / nrow(data1))
gamma <- 0.15 # A small value of gamma favours an L with a small rank
l1 <- lambda * gamma
l2 <- lambda * (1 - gamma)

fit <- lrpsadmm(
  Sigma = cor(data1),
  Lambda1 = l1,
  Lambda2 = l2,
  maxiter = 1000,
  print_progress = FALSE,
  zeros = NULL,
  backend = "RcppEigen"
)

S <- fit$S # The estimated Shat (sparse precision matrix)
L <- fit$L # The estimated Lhat (dense low-rank matrix)
K <- Matrix::rankMatrix(L)[1] # K

Y9 <- quiet(generate.data(
  Sest = solve(S),
  n = nrow(data1),
  p = ncol(data1)
))
colnames(Y9) <- colnames(data1)
Y9<- as.matrix(Y9)


#...lrpsadmml (data0)

lambda <- sqrt(log(ncol(data0)) / nrow(data0))
gamma <- 0.15 # A small value of gamma favours an L with a small rank
l1 <- lambda * gamma
l2 <- lambda * (1 - gamma)

fit <- lrpsadmm(
  Sigma = cor(data0),
  Lambda1 = l1,
  Lambda2 = l2,
  maxiter = 1000,
  print_progress = FALSE,
  zeros = NULL,
  backend = "RcppEigen"
)

S <- fit$S # The estimated Shat (sparse precision matrix)
L <- fit$L # The estimated Lhat (dense low-rank matrix)
K <- Matrix::rankMatrix(L)[1] # K

Y9 <- quiet(generate.data(
  Sest = solve(S),
  n = nrow(data0),
  p = ncol(data0)
))
colnames(Y9) <- colnames(data0)
Y9<- as.matrix(Y9)


# SEM fitting

Y <- list(
  list(g1, data1),
  list(bap1$dag, bap1$data),
  list(bap2$dag, bap2$data),
  list(bap3$dag, bap3$data),
  list(bap4$dag, bap4$data),
  list(bap5$dag, bap5$data),
  list(bap6$dag, bap6$data),
  list(bap7$dag, bap7$data),
  list(bap8$dag, bap8$data),
  list(g1, Y9)
)

names(Y) <-
  c(
    "unadjusted",
    "cggm_dsep",
    "glpc_dsep",
    "pc",
    "pcss",
    "trim",
    "lava",
    "rsvp",
	"glasso",
    "lrpsadmm"
  )

J <- length(Y);J

results <- data.frame(
  method = names(Y),
  srmr = rep(NA, length(Y)),
  dev_df = rep(NA, length(Y)),
  node_act = rep(NA, length(Y)),
  node_inh = rep(NA, length(Y)),
  nlog10P = rep(NA, length(Y)),
  vcountP = rep(NA, length(Y))
)

fit <- list()

for (i in 1:J) { #i=1
  print(names(Y)[[i]]);cat("\n")
  fit <-
    quiet(SEMrun(
      graph = Y[[i]][[1]],
      data = Y[[i]][[2]],
      group = group,
      algo = "ricf"
    ))
  fitIdx <- fit$fit$fitIdx
  results$srmr[i] <- fitIdx[3]
  results$dev_df[i] <- fitIdx[1] / fitIdx[2]
  results$node_act[i] <- fit$pval[1]
  results$node_inh[i] <- fit$pval[2]
  results$nlog10P[i] <- -log10(2*min(fit$pval))
  if (i == 3 | i == 4) {
   fit$graph <- delete_vertices(fit$graph,
    V(fit$graph)$name[grep("LV", V(fit$graph)$name)])
  }
  results$vcountP[i] <- sum(V(fit$graph)$color != "white")
}

class(results) <- c("lavaan.data.frame","data.frame")
results


gg <- lapply(fit[1:J], function(x) x$graph)

for (i in 1:J) {
  if (i == 3 | i == 4) {
    gg[[i]] <-
      delete_vertices(gg[[i]], V(gg[[i]])$name[grep("LV", V(gg[[i]])$name)])
  }
  gplot(gg[[i]], l = "dot", main = names(Y)[i])
  Sys.sleep(3)
}


# END !


# deconfounding scatter plot

rownames(dataplot) <-NULL
col <- colnames(dataplot[,-1])

colnames(bap1$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap1$data), 'SYMBOL', 'ENTREZID'))
colnames(bap2$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap2$data), 'SYMBOL', 'ENTREZID'))
colnames(bap3$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap3$data), 'SYMBOL', 'ENTREZID'))
# colnames(bap4$data) <-
#   as.character(mapIds(org.Hs.eg.db, colnames(bap4$data), 'SYMBOL', 'ENTREZID'))
# colnames(bap5$data) <-
#   as.character(mapIds(org.Hs.eg.db, colnames(bap5$data), 'SYMBOL', 'ENTREZID'))
colnames(bap6$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap6$data), 'SYMBOL', 'ENTREZID'))
colnames(bap7$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap7$data), 'SYMBOL', 'ENTREZID'))
colnames(bap8$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap8$data), 'SYMBOL', 'ENTREZID'))
colnames(bap9$data) <-
  as.character(mapIds(org.Hs.eg.db, colnames(bap9$data), 'SYMBOL', 'ENTREZID'))
colnames(Y10) <-
  as.character(mapIds(org.Hs.eg.db, colnames(Y10), 'SYMBOL', 'ENTREZID'))

Y1 <- bap1$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y2 <- bap2$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y3 <- bap3$data %>%
  as.data.frame() %>%
  dplyr::select(col)
# Y4 <- Y4 %>%
#   as.data.frame() %>%
#   dplyr::select(col)
# Y5 <- Y5 %>%
#   as.data.frame() %>%
#   dplyr::select(col)
Y6 <- bap6$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y7 <- bap7$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y8 <- bap8$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y9 <- bap9$data %>%
  as.data.frame() %>%
  dplyr::select(col)
Y10 <- Y10 %>%
  as.data.frame() %>%
  dplyr::select(col)

dataplot$method <- "Unadjusted"
Y1$method <- "CGGM"
Y2$method <- "pcor"
Y3$method <- "gLASSO"
# Y4$method <- "PCA"
# Y5$method <- "gLPCA"
Y6$method <- "PCSS"
Y7$method <- "Trim"
Y8$method <- "Lava"
Y9$method <- "rvsp"
Y10$method <- "LRpS"

Y1$TCF7L2 <- dataplot$TCF7L2
Y2$TCF7L2 <- dataplot$TCF7L2
Y3$TCF7L2 <- dataplot$TCF7L2
# Y4$TCF7L2 <- dataplot$TCF7L2
# Y5$TCF7L2 <- dataplot$TCF7L2
Y6$TCF7L2 <- dataplot$TCF7L2
Y7$TCF7L2 <- dataplot$TCF7L2
Y8$TCF7L2 <- dataplot$TCF7L2
Y9$TCF7L2 <- dataplot$TCF7L2
Y10$TCF7L2 <- dataplot$TCF7L2

Y1$ESR2 <- dataplot$ESR2
Y2$ESR2 <- dataplot$ESR2
Y3$ESR2 <- dataplot$ESR2
# Y4$ESR2 <- dataplot$ESR2
# Y5$ESR2 <- dataplot$ESR2
Y6$ESR2 <- dataplot$ESR2
Y7$ESR2 <- dataplot$ESR2
Y8$ESR2 <- dataplot$ESR2
Y9$ESR2 <- dataplot$ESR2
Y10$ESR2 <- dataplot$ESR2

# 4-5

LV4 <- bap4$data[, grep("LV", colnames(bap4$data))]; rownames(LV4) <- NULL
LV5 <- bap5$data[, grep("LV", colnames(bap5$data))]; rownames(LV5) <- NULL

fit4 <- SEMrun(graph = bap4$dag, data = bap4$data, group = group, algo = "ricf")
fit5 <- SEMrun(graph = bap5$dag, data = bap5$data, group = group, algo = "ricf")

B4 <- fit4$fit$parameterEstimates
B5 <- fit5$fit$parameterEstimates

col_entrez <-  as.character(mapIds(org.Hs.eg.db, col, 'ENTREZID', 'SYMBOL'))
B4 <- B4 %>%
  dplyr::filter(lhs %in% col_entrez & grepl('LV', rhs))
B5 <- B5 %>%
  dplyr::filter(lhs %in% col_entrez & grepl('LV', rhs))

B4_1 <- B4 %>% dplyr::filter(lhs == col_entrez[1]) %>% pull()
B4_2 <- B4 %>% dplyr::filter(lhs == col_entrez[2]) %>% pull()

B5_1 <- B5 %>% dplyr::filter(lhs == col_entrez[1]) %>% pull()
B5_2 <- B5 %>% dplyr::filter(lhs == col_entrez[2]) %>% pull()

B4_LV_1 <- matrix(NA, nrow(LV4), ncol(LV4))
B4_LV_2 <- matrix(NA, nrow(LV4), ncol(LV4))
B5_LV_1 <- matrix(NA, nrow(LV5), ncol(LV5))
B5_LV_2 <- matrix(NA, nrow(LV5), ncol(LV5))

for (i in 1:ncol(LV4)){
  B4_LV_1[,i] <- B4_1[i] %*% LV4[,i]
  B4_LV_2[,i] <- B4_2[i] %*% LV4[,i]
  
  B5_LV_1[,i] <- B5_1[i] %*% LV5[,i]
  B5_LV_2[,i] <- B5_2[i] %*% LV5[,i]
}

orig_1 <- dataplot %>% dplyr::select(col[1])
orig_2 <- dataplot %>% dplyr::select(col[2])

B4_LV_1 <-  B4_LV_1 %>% rowSums(na.rm = TRUE)
B4_LV_2 <-  B4_LV_2 %>% rowSums(na.rm = TRUE)
Y4 <- data.frame(orig_1 - B4_LV_1, orig_2 - B4_LV_2); colnames(Y4) <- col

B5_LV_1 <-  B5_LV_1 %>% rowSums(na.rm = TRUE)
B5_LV_2 <-  B5_LV_2 %>% rowSums(na.rm = TRUE)
Y5 <- data.frame(orig_1 - B5_LV_1, orig_2 - B5_LV_2); colnames(Y5) <- col

Y4$method <- "PCA"
Y5$method <- "gLPCA"

Y4$TCF7L2 <- dataplot$TCF7L2
Y5$TCF7L2 <- dataplot$TCF7L2

Y4$ESR2 <- dataplot$ESR2
Y5$ESR2 <- dataplot$ESR2

data <- rbind(dataplot, Y1, Y3, Y4, Y5, Y6, Y7, Y10)

# data_long <- tidyr::gather(Y1, EGFR, GRB2, factor_key=TRUE)
# data_long
# 
# Y1 <- melt(setDT(Y1), id.vars = c("EGFR","GRB2"), variable.name = "value")
pal <- wes_palette(name = "GrandBudapest2", type = "continuous")

pal <- RColorBrewer::brewer.pal(n = 20,name="Set3")[4:11]
# rm <- c("gLPCA", "PCA")
data$method <- factor(data$method, levels=c('Unadjusted','CGGM','gLPCA',
                                            'PCA',"Trim", "PCSS", 'LRpS', "gLASSO"))
library(grid)

setEPS()
file <- paste0("tfs1.eps")
postscript(width=13, height = 3, file) #14,12

data %>%
  # dplyr::filter(!method %in% rm) %>%
  ggplot(aes(
    x = ESR2, #ESR2 #TCF7L2
    y = FZD4, #(EGFR;GRB2) #(FZD4; PIK3R2)
    color = method
  )) +
  # geom_smooth(method='lm', formula= y~x ) +
  geom_point(size = 2) +
  #stat_summary(fun.y=mean, geom="point", position="identity") +
  # facet_wrap( ~ graph, ncol = 1, scales="free_y") +
  facet_grid(. ~ method, scales="free_y") +
  theme_light() +
  xlab("Latent ESR2")+
  ylab("FZD4") +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 0)) +
  theme(legend.title = element_text(size = 9)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "") +
  theme(legend.key.size = unit(0.3, "cm")) +
  # scale_fill_manual(values=wes_palette(n=9, name="GrandBudapest2")) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(text = element_text(size = 17)) +
  scale_color_brewer(palette="Set2") +
  ylim(-2.5, 2.5) 
dev.off()

setEPS()
file <- paste0("tfs2.eps")
postscript(width=13, height = 3, file) #14,12
data %>%
  # dplyr::filter(!method %in% rm) %>%
  ggplot(aes(
    x = ESR2, #ESR2 #TCF7L2
    y = PIK3R2, #(EGFR;GRB2) #(FZD4; PIK3R2)
    color = method
  )) +
  # geom_smooth(method='lm', formula= y~x ) +
  geom_point(size = 2) +
  #stat_summary(fun.y=mean, geom="point", position="identity") +
  # facet_wrap( ~ graph, ncol = 1, scales="free_y") +
  facet_grid(. ~ method, scales="free_y") +
  theme_light() +
  xlab("Latent ESR2")+
  ylab("PIK3R2") +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 0)) +
  theme(legend.title = element_text(size = 9)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "") +
  theme(legend.key.size = unit(0.3, "cm")) +
  # scale_fill_manual(values=wes_palette(n=9, name="GrandBudapest2")) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(text = element_text(size = 17)) +
  scale_color_brewer(palette="Set2") +
  ylim(-2.5, 2.5) 

dev.off()
