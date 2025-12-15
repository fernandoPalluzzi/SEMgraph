
generate_sim <- function (mu = mu){
  
  # cat("Generating covariance matrix......")
  set1_matrix <- matrix(0, 500, 500)
  for (i in 1:500) {
    for (j in 1:500) {
      if (i == j) 
        set1_matrix[i, j] = 1
    }
  }
  set2_matrix <- matrix(0, 500, 500)
  for (i in 1:500) {
    for (j in 1:500) {
      if (i == j) 
        set2_matrix[i, j] = 1
      if (i <= 50 & j <= 50 & i != j) 
        set2_matrix[i, j] = 0.6
    }
  }
  set3_matrix <- set1_matrix
  set4_matrix <- set2_matrix
  set5_matrix <- matrix(0, 500, 500)
  for (i in 1:500) {
    for (j in 1:500) {
      if (i == j) 
        set5_matrix[i, j] = 1
      if (i <= 50 & j <= 50 & i != j) {
        if ((i <= 25 & j <= 25) | (i > 25 & j > 25)) 
          set5_matrix[i, j] = 0.6
        else set5_matrix[i, j] = -0.6
      }
    }
  }
  # cat("Simulating data......")
  set1_data <- mvrnorm(n = 20, rep(0, 500), Sigma = set1_matrix) 
  set2_data <- mvrnorm(n = 20, c(rep(mu, 50), rep(0, 450)), 
                       Sigma = set2_matrix)
  set3_data <- mvrnorm(n = 20, c(rep(mu, 50), rep(0, 450)), 
                       Sigma = set3_matrix)
  set4_data <- mvrnorm(n = 20, rep(0, 500), Sigma = set4_matrix) 
  set5_data <- mvrnorm(n = 20, c(rep(mu, 25), rep(-mu, 
                                                  25), rep(0, 450)), Sigma = set5_matrix) 
  set6_data <- cbind(set2_data[, 1:10], set3_data[, 1:10], 
                     set4_data[, 1:10], set5_data[, 1:10], set1_data[, 1:460]) 
  control_data <- mvrnorm(n = 20, rep(0, 500), Sigma = set1_matrix) 
  # cat("Done!")
  return(list(set1_data, set2_data, set3_data, set4_data, 
              set5_data, set6_data, control_data))
}

#gse_filename="GSE172114_series_matrix.txt.gz"; gse=NULL; data.type="rseq"

getdata <- function(gse_filename, gse = NULL, data.type = "rseq"){
	
	GSE<- quiet(getGEO(file=gse_filename))
	
	phenodata <- Biobase::pData(GSE@phenoData)
	
	phenodata  <- phenodata  %>%
		dplyr::select(title,'disease state:ch1') %>%
		dplyr::rename(Sample = title,
					Group = "disease state:ch1") %>%
		dplyr::mutate(Group = case_when(Group == "Non-critical" ~ 0,
										Group == "Critical" ~ 1),
					Sample = gsub("\\-.*","",Sample))
	
	group <- phenodata$Group
	
	X<- read.csv("./data/GSE172114_rsem_gene_count_matrix_TMM_69samples.csv")
	
	ENSEMBL<- strsplit(X[,1], "[.]")
	genes<- NULL
	for(i in 1:nrow(X)) genes[i]<- ENSEMBL[[i]][[1]]
	ENTREZID<- quiet(mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'ENSEMBL'))
	
	X<- na.omit(cbind(ENTREZID, X[,-1]))
	del<- which(duplicated(X$ENTREZID) == TRUE)
	x<- X[-del,-1]
	rownames(x)<- X$ENTREZID[-del] 
	geneSE <- x

	disease = "Coronavirus disease - COVID-19"
	
	path <- which(names(kegg.pathways) == disease)
	ig <- as.undirected(kegg.pathways[[path]])
	
	# Get SummarizedExperiment format 
	y<- group; x<- geneSE
	
	write.table(x, file="xdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	write.table(cbind(colnames(x),y), file="cdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	write.table(rownames(x), file="rdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	
	exprs.file <- file.path(getwd(), "xdat.tab")
	cdat.file <- file.path(getwd(), "cdat.tab")
	rdat.file <- file.path(getwd(), "rdat.tab")
	
	se <- readSE(exprs.file, cdat.file, rdat.file, data.type)

	# Log2 transformation
	if (data.type == "ma") assay(se) <- l2t(se)
	 
	# DEA
	se <- deAna(se)
	#if the specified assay contains the *raw* read counts: 
	#se <- deAna(se, de.method = "limma", padj.method = "BH", assay = "raw")
	#se <- deAna(se, de.method = "edgeR", assay = "raw")
	#se <- deAna(se, de.method = "DESeq2", assay = "raw")

	return(data = list(graph = ig, se = se, group = group, 
					   phenodata = phenodata, disease = disease))
}


l2t <- function(se){
	
	# Log2 transformation
	ex <- assay(se)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
	if (LogC) {
	 ex[which(ex <= 0)] <- NaN
	 ex <- log2(ex)
	}
	 
	return(ex)
}

#graph=network; data=data; group=group; alpha=0.05; maxsize=200

WalktrapGM <- function(graph, data, group, alpha=0.05, maxsize=200, ...)
{
	# nodes and edge weighting
	graph <- weightGraph(graph, data, group, method="r2z")
	W <- as_adjacency_matrix(graph, attr=NULL, sparse=FALSE)
	p0 <- ifelse(p.adjust(V(graph)$pv, method="BH") < alpha, 1, 0)
	V(graph)$score <- diffusr::random.walk(p0, W, r=0.5)$p.inf
	vweights <- qnorm(1-V(graph)$pv/2)*V(graph)$score
	eweights <- qnorm(1-E(graph)$pv/2) #qnorm(1-0.05/2)

	# Bootstrap node samples
	dist <- matrix (ncol = 5000, nrow = maxsize)  
	for (i in 3:maxsize) {  
	 drawsDS1 <- matrix(sample(vweights, size = i * 5000, replace = TRUE), i)  
	 drawsSUM <- apply(drawsDS1, 2, sum)  
	 drawsW <- sqrt(sum((drawsSUM)^2))
	 drawsSCORE <- drawsSUM/drawsW
	 dist[i,] <- drawsSCORE
	}

	#Community structure via short random walks
	wtc <- igraph::cluster_walktrap(graph, weights=eweights, steps=5)
	memb <- igraph::membership(wtc) # table(memb)
	csize <- igraph::sizes(wtc) # length(csize)

	#Module scores
	resMat <- NULL  
	cvector <- which(csize >= 5 & csize <= maxsize) # length(cvector)
	for (i in cvector) { #i=cvector[1]
	 comm.total.weight <- sqrt(sum((vweights[which(memb == i)])^2))
	 comm.total.mean <- sum(vweights[which(memb == i)])/comm.total.weight 
	 boot.total.mean <- mean(dist[csize[i],])
	 boot.total.sd <- sd(dist[csize[i],])
	 comm.total.tscore <- abs(comm.total.mean-boot.total.mean)/boot.total.sd
	 resMat <- rbind(resMat, c(i, comm.total.tscore, csize[i]))  
	}
	colnames(resMat) <- c("ID", "Score", "Size")
	resMat <- as.data.frame(na.omit(resMat))
	resMat <- resMat[order(resMat$Score, decreasing=TRUE),]
	memb <- memb[memb %in% names(csize)[cvector]] # length(table(memb))

	return(list(scores=resMat, membership=memb))
}

#graph=network; limma=rowData(se); alpha=0.05; maxsize=200

WalktrapGMl <- function(graph, limma, alpha= 0.05, maxsize=200, ...)
{
	# Nodes weighting:
	V <- rownames(limma)[rownames(limma) %in% V(graph)$name]
	FC <- limma[V,2] 
	names(FC) <- limma[V,1]
	# head(FC); length(FC)
	vweights <- c()   
	for (i in 1:vcount(graph)) {
	 vertex <- V(graph)$name[i]
	 vw <- abs(FC[which(names(FC) == vertex)])
	 if (length(vw) != 0){
	   vweights[i] <- vw
	 } else {
	   vweights[i] <- NA
	 }
	 #vweights[i] <- abs(FC[which(names(FC) == vertex)])
	}
	vweights[is.na(vweights)] <- .01
	# Edge weighting:
	ftm <- as_edgelist(graph)
	eweights <- c()
	for (i in 1:nrow(ftm)) {
	 from <- ftm[i,1]  
	 to <-  ftm[i,2]  
	 ew <- (abs(FC[which(names(FC) == from)]) + abs(FC[which(names(FC) == to)]))/2
	 if (length(ew) != 0){
	   eweights[i] <- ew
	 } else {
	   eweights[i] <- NA
	 }
	 #eweights[i] <- (abs(FC[which(names(FC) == from)]) + abs(FC[which(names(FC) == to)]))/2
	}
	eweights[is.na(eweights)] <- .01

	# Bootstrap node samples:
	dist <- matrix (ncol = 5000, nrow = maxsize)  
	for (i in 3:maxsize) {  
	 drawsDS1 <- matrix(sample(vweights, size = i * 5000, replace = TRUE), i)  
	 drawsSUM <- apply(drawsDS1, 2, sum)  
	 drawsW <- sqrt(sum((drawsSUM)^2))
	 drawsSCORE <- drawsSUM/drawsW
	 dist[i,] <- drawsSCORE
	}

	# Community structure via short random walks:
	wtc <- igraph::cluster_walktrap(graph, weights=eweights, steps=5)
	memb <- igraph::membership(wtc) # table(memb)
	csize <- igraph::sizes(wtc) # length(csize)

	# Module scores:
	resMat <- NULL  
	cvector <- which(csize >= 5 & csize <= maxsize) # length(cvector)
	for (i in cvector) { #i=cvector[1]
	 comm.total.weight <- sqrt(sum((vweights[which(memb == i)])^2))
	 comm.total.mean <- sum(vweights[which(memb == i)])/comm.total.weight 
	 boot.total.mean <- mean(dist[csize[i],])
	 boot.total.sd <- sd(dist[csize[i],])
	 comm.total.tscore <- abs(comm.total.mean-boot.total.mean)/boot.total.sd
	 resMat <- rbind(resMat, c(i, comm.total.tscore, csize[i]))  
	}
	colnames(resMat) <- c("ID", "Score", "Size")
	resMat <- as.data.frame(na.omit(resMat))
	resMat <- resMat[order(resMat$Score, decreasing=TRUE),]
	memb <- memb[memb %in% names(csize)[cvector]] # length(table(memb))

	return(list(scores=resMat, membership=memb))
}

#mod=mod6; data=data; group=group; mmx=mmx; h=0; alpha=0.025; Vids=Vids

extract_tree <- function(mod, data, group, mmx, h, alpha, Vids, ...)
{
	if (h > 0){
	 T <- SEMtree(graph=NULL, data=data, seed=V(mmx$gLM)$name, type="CAT")
	 #T <- SEMtree(graph=NULL, data=data, seed=V(mmx$gLM)$name, type="CPDAG")
	}else{
     T <- SEMtree(graph=NULL, data=data, seed=V(mod)$name, type="CAT")
	 #T <- SEMtree(graph=NULL, data=data, seed=V(mod)$name, type="CPDAG")
	}
	
	E(T)$color <- "gray60"
    V(T)$color <- "white"
	if (h > 0){
	 clm <- table(mmx$membership)
	 clm <- gsub("p","", names(clm))
	 clv <- table(mmx$membership[names(mmx$membership) %in% Vids])
     clv <- gsub("p","", names(clv))
	 V(T)$color[which(V(T)$name %in% clm)] <- "green"
	 V(T)$color[which(V(T)$name %in% clv)] <- "orange"
     V(T)$color[which(V(T)$name %in% Vids)] <- "yellow"
	}else{
     V(T)$color[which(V(T)$name %in% Vids)] <- "yellow"
	}
	
	#library(org.Hs.eg.db)
	V(T)$label <- mapIds(org.Hs.eg.db, V(T)$name, 'SYMBOL', 'ENTREZID')

	return(T)
}

jaccard <- function(x, ...)
{
	# initialize similarity matrix
	n<- length(x)
	jaccard<- matrix(NA, nrow=n, ncol=n, dimnames=list(names(x),names(x)))

	# compute J index for pairwise vectors of a list  
	for(i in 1:(n-1)) { 
	 for(j in (i+1):n) {
		jaccard[i,j]= length(intersect(x[[i]],x[[j]])) / length(union(x[[i]],x[[j]]))
		jaccard[j,i]=jaccard[i,j]
	 }
	}
	diag(jaccard)<-1

	return(jaccard)
}

gene_metrics<- function(tree, goldgene,...) #tree=xx[[1]]; goldgene=cov
{
  seed <- V(tree)$name
  sel <- V(tree)$color == "yellow" | V(tree)$color == "orange"
  seedcov <- length(V(tree)$name[sel])
  
  size_genes <- length(seed)
  gene_pre <- seedcov/length(seed)
  gene_rec <- seedcov/length(goldgene)
  gene_f1 <- 2*(gene_pre*gene_rec)/(gene_pre+gene_rec)
  
  return(data.frame(size_genes, seedcov, gene_pre, gene_rec, gene_f1))
}

GO_metrics<- function(tree, goldGO, ...) 
{
  seed <- V(tree)$name
  ego <- enrichGO(gene = seed,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable = TRUE)
  
  all.GO <- ego@result
  nsel <- sum(rownames(all.GO) %in% rownames(goldGO))
  
  GO_pre <- nsel/nrow(all.GO)
  GO_rec <- nsel/nrow(goldGO)
  GO_f1 <- 2*(GO_pre*GO_rec)/(GO_pre+GO_rec)
  
  return(data.frame(nrow(all.GO), nsel, GO_pre, GO_rec, GO_f1))
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}
