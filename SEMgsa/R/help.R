# Get GEO data -----------------------------------------------------------------

getGEOdata <- function(GSEfilename, GSEformat="matrix"){
  
  # Get system enviroment
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
  
  # Choose download format and transform GSE in a dataset 
  if(GSEformat == "soft"){
    
    GSE <- getGEO(filename = GSEfilename, GSEMatrix = FALSE, 
                  AnnotGPL = TRUE, parseCharacteristics = FALSE)
    
    # Check the number of platforms
    gsmplatforms <- lapply(GSMList(GSE), function(x) {Meta(x)$platform})
    n_platform <- length(table(unlist(gsmplatforms)))
    
    # Get probeset data from platforms (if not unique)
    if (n_platform != 1){
      gList = gData = list()
      p = list()
      f = list()
      for(i in (1:n_platform)){ 
        
        p_platform = data.frame()
        id_platform <- unique(gsmplatforms)[[i]]
        id_samples <- names(gsmplatforms[gsmplatforms == id_platform])
        gList[[i]] <- lapply(GSMList(GSE)[id_samples], function(x) Table(x)[,2])
        gData[[i]] <- do.call(cbind, gList[[i]])
        rownames(gData[[i]]) <- Table(GSMList(GSE)[id_samples][[1]])$ID_REF
        names(gData)[[i]] <- id_platform
        
        for(j in colnames(gData[[i]])){
          
          # Get phenoData
          psample <- data.frame(GSE@gsms[[j]]@header)[1,]
          p_platform <- rbind.fill(p_platform,psample)
          
        }
        rownames(p_platform) <- colnames(gData[[i]])
        p[[i]] <- p_platform
        names(p)[[i]] <- id_platform
        
        # Get geneAnnot
        f[[i]] <- GSE@gpls[[id_platform]]@dataTable@table
        names(f)[[i]] <- id_platform
        
      }} else {
        p = data.frame()
        gList <- lapply(GSMList(GSE), function(x) Table(x)[,2])
        gData <- do.call(cbind, gList)
        rownames(gData) <- Table(GSMList(GSE)[[1]])$ID_REF
        
        for(j in colnames(gData)){
          
          # Get phenoData
          psample <- data.frame(GSE@gsms[[j]]@header)[1,]
          p <- rbind.fill(p,psample)
          
        }
        rownames(p) <- colnames(gData)
        
        # Get geneAnnot
        f <- GSE@gpls[[1]]@dataTable@table
      }
    
  }else if(GSEformat == "matrix"){
    
    # Get GEO object
    GSE <- getGEO(filename = GSEfilename, GSEMatrix = TRUE, 
                  AnnotGPL= TRUE, parseCharacteristics = FALSE)
    # Get phenoData
    p<- pData(GSE) 
    # Get geneAnnot
    f<- fData(GSE) 
    # Get expression data
    gData<- exprs(GSE) 
  } 
  
  return(list(GSE = GSE, phenoData = p, geneAnnot = f, exprsData = gData))
}

# Probe2gene conversion --------------------------------------------------------

probe2gene <- function(GPLdata, GPLannot, geneID = "entrez",
                       summary_probes = "random"){
  
  #1) Assigning new identifiers to rownames (entrez or symbol)
  if ( is.character(GPLannot) ){
    GPLannot<- getGEO(file=GPLannot, parseCharacteristics = FALSE)
    table_annot<- Table(GPLannot) #colnames(table_annot)
  } else {
    table_annot<- GPLannot
  }
  probeset<- table_annot$ID
  
  if(length(probeset) == length(rownames(GPLdata))){
    if (geneID == "entrez"){
      if (("Gene ID" %in% colnames(table_annot)) == TRUE){
        row_annotations<- as.vector(table_annot$"Gene ID")
      } else if (("ENTREZ_GENE_ID" %in% colnames(table_annot)) == TRUE){
        row_annotations<- as.vector(table_annot$"ENTREZ_GENE_ID")
      }} else if (geneID == "symbol"){
        row_annotations<- as.vector(table_annot$"Gene Symbol")
      }
    sel<- which(probeset %in% rownames(GPLdata))
    rownames(GPLdata)<-row_annotations[sel]
    # rownames(GPLdata)<- make.names(row_annotations[sel], unique=TRUE)
    # rownames(GPLdata) = substring(rownames(GPLdata), 2)
  } else {
    message("n.probes sample NE n.probes plastform")
    return( NULL)
  }
  
  #2) Delete the number of probes that present ID= " " and "///"
  del<- c(which(rownames(GPLdata)==""), grep("///",rownames(GPLdata)))
  GPLdata1<- GPLdata[-del,] 
  
  #3) Conversion from probeset ID to gene IDs
  if (summary_probes == "random"){
    #a) Remove duplicate probes
    del<- which(duplicated(rownames(GPLdata1))) 
    GPLdata2<- GPLdata1[-del,]
  } else {
    #b) Take mean over duplicated probes
    tmp <- unique(rownames(GPLdata1))
    forder <- tmp[order(tmp)]
    forder <- forder[!is.na(forder)]
    GPLdata2 <- data.frame(matrix(nrow=length(forder),ncol=ncol(GPLdata1)))
    for (i in 1:length(forder)){
      ind <- which(rownames(GPLdata1)==forder[i])
      if (length(ind) == 1){
        GPLdata2[i,] <- GPLdata1[ind,]
      } else{
        GPLdata2[i,] <- apply(GPLdata1[ind,], 2, summary_probes)
      }
    }
    rownames(GPLdata2) <- forder; colnames(GPLdata2) <- colnames(GPLdata1)
  }
  
  return(GPLdata2)
}

# Log2 transformation ----------------------------------------------------------

l2t <- function(se){
  
  # Log2 transformation
  ex <- assay(se)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  assay(se) <- log2(ex) }
  
  return(se)
}

# Import (and convert) data from GEO -------------------------------------------

getGSAdata <- function(gse_filename, gse = NULL, summary_probes = "random", log = TRUE){
  
  if(gse == "GSE35570"){
    
    probeSE <- quiet(getGEOdata(GSEfilename = gse_filename,
                                GSEformat="matrix"))
    phenodata <- probeSE[["phenoData"]]
    annot <- probeSE[["geneAnnot"]]
    ex <- probeSE[["exprsData"]]
    
    phenodata  <- phenodata  %>%
      dplyr::select("characteristics_ch1") %>%
      dplyr::rename(Group = "characteristics_ch1") %>%
      dplyr::mutate(Group = case_when(Group %like% "normal" ~ 0,
                                      Group %like% "PTC" ~ 1),
                    Sample = rownames(phenodata)) %>%
      arrange(desc(Group))
    
    group <- phenodata$Group
    ex <- ex[, phenodata$Sample]
    disease <- "Thyroid cancer"
    
    geneSE <- probe2gene(GPLdata=ex, GPLannot =annot, summary_probes= summary_probes)
    path <- which(names(kegg.pathways) == disease)
    # ig  <- kegg.pathways[[path]]
    ig <- quiet(properties(kegg.pathways[[path]])[[1]])
    
  } else if (gse == "GSE172114"){
    
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
    
    X<- read.csv("data/GSE172114_rsem_gene_count_matrix_TMM_69samples.csv")
    
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
    ig <- quiet(properties(kegg.pathways[[path]])[[1]])
    
  } else if (gse == "GSE6956"){
    probeSE <- quiet(getGEOdata(GSEfilename = gse_filename,
                                GSEformat="matrix"))
    phenodata <- probeSE[["phenoData"]]
    annot <- probeSE[["geneAnnot"]]
    ex <- probeSE[["exprsData"]]
    
    phenodata <- phenodata %>%
      dplyr::select("source_name_ch1" ) %>%
      dplyr::rename(Group = "source_name_ch1" ) %>%
      dplyr::mutate(Sample = rownames(phenodata),
                    Group = case_when(Group %like% "Normal" ~ 0,
                                      Group %like% "Adenocarcinoma" ~ 1))
    # ex <- exprs(probeSE)
    group <- phenodata$Group
    disease <- "Prostate cancer"
    
    geneSE <- probe2gene(GPLdata=ex, GPLannot =annot, summary_probes= summary_probes)
    path <- which(names(kegg.pathways) == disease)
    # ig  <- kegg.pathways[[path]]
    ig <- quiet(properties(kegg.pathways[[path]])[[1]])
    
  } else if (gse == "GSE68719"){
    
    Sys.setenv("VROOM_CONNECTION_SIZE" = 262144*2) 
    GSE<- quiet(getGEO(file=gse_filename))
    
    phenodata <- Biobase::pData(GSE@phenoData)
    phenodata  <- phenodata  %>%
      dplyr::select(title) %>%
      dplyr::rename(Group = title) %>%
      dplyr::mutate(Group = case_when(Group %like% "C" ~ 0,
                                      Group %like% "P" ~ 1),
                    Sample = rownames(phenodata)) %>%
      arrange(Group)
    
    group <- phenodata$Group
    
    expr<- read.table("data/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz", 
                      header = TRUE, row.names = 1, sep = "\t")
    
    counts <- expr[,-1] 
    # Create a DGEList object for storing read counts and associated information
    # d <- DGEList(counts=counts, group=group)
    
    # Filter really lowly expressed genes
    # keep <- filterByExpr(d)
    # d <- d[keep,,keep.lib.sizes=FALSE]
    # str(d) #15709 73
    
    # Calculate normalization factors to scale the raw library sizes.
    # Normalization may heavily affect calculation by eliminating artificial expression
    # levels due to library size or sequencing depth. Default method = "TMM"
    # d <- calcNormFactors(d, method="TMM") 
    # head(d$samples)
    
    # Voom transforms count data to log2-counts per million (logCPM),... 
    voom <- limma::voom(counts)$E
    # boxplot(voom)
    
    genes <- expr[,1]
    ENTREZID<- quiet(mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL'))
    del<- which(duplicated(ENTREZID) == TRUE | is.na(ENTREZID))
    voom<- voom[-del,]
    rownames(voom)<- ENTREZID[-del]
    colnames(voom) <- phenodata$Sample
    geneSE <- voom
    
    disease <- "Parkinson disease"
    path <- which(names(kegg.pathways) == disease)
    # ig  <- kegg.pathways[[path]]
    ig <- quiet(properties(kegg.pathways[[path]])[[1]])
  } else if (gse == "GSE53740"){
    
    group <- ftdDNAme$group
    disease <- "Frontotemporal Dementia (FTD)"
    
    ftd.pathways <- c("MAPK signaling pathway",
                      "Protein processing in endoplasmic reticulum",
                      "Endocytosis",
                      "Wnt signaling pathway",
                      "Notch signaling pathway",
                      "Neurotrophin signaling pathway")
    path <- which(names(kegg.pathways) %in% ftd.pathways)
    ig <- kegg.pathways[path] 
    
    geneSE <- t(quiet(huge.npn(ftdDNAme$pc1)))
    phenodata <- data.frame(Sample = colnames(geneSE),
                            Group = group)
    
  }
  
  # Get SummarizedExperiment format 
  y<- group; x<- geneSE
  
  write.table(x, file="xdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(cbind(colnames(x),y), file="cdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(rownames(x), file="rdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  
  exprs.file <- file.path(getwd(), "xdat.tab")
  cdat.file <- file.path(getwd(), "cdat.tab")
  rdat.file <- file.path(getwd(), "rdat.tab")
  
  se <- EnrichmentBrowser::readSE(exprs.file, cdat.file, rdat.file)
  if (gse == "GSE172114" & gse == "GSE68719") se <- readSE(exprs.file, cdat.file, rdat.file, data.type = "rseq")
  
  if (log == TRUE & gse != "GSE3678" & gse != "GSE68719" & gse != "GSE53740"){
    # Log2 transformation
    se <- l2t(se)
    geneSE <- assay(se)
  }
  
  # DEA
  se <- deAna(se)
  
  return(data = list(graph = ig, geneSE = geneSE, se = se, group = group, 
                     phenodata = phenodata, disease = disease))
}


# Pathway deregulation ---------------------------------------------------------
## For a given pathway, compute centrality measures and choose a subset of genes
# to be affected with the desired mean changes (choose DC).
# DC: detection call, a percentage between 0% and 100% indicating the proportion
# of genes to be affected

## Community pathway deregulation
## For a given pathway, find communities and choose a subset of communities such that
## genes in these communities are affected with the desired mean changes. 

community.dereg <- function(pathway = NULL, n = 10, seed=1){
  
  set.seed(seed)
  adj <- as.matrix(igraph::get.adjacency(pathway, attr = "weight"))
  adj <- colSums(adj)
  genes <- ifelse(adj>=1, 1, ifelse(adj==0, 0, -1))
  genes_weight <- genes[genes != 0]
  genes <- names(genes_weight)
  size <- length(genes)
  
  # genes <- V(pathway)
  # genes <- as_ids(genes)
  # size <- length(genes)
  el <- data.frame(as_edgelist(pathway))
  colnames(el) <- c("src","dest")
  
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  lec <- igraph::cluster_edge_betweenness(gr)
  tmp <- as.numeric(table(lec$membership))
  proportion <- tmp/size
  names(proportion) <- names(table(lec$membership))
  proportion <- sort(proportion, decreasing = T)
  mem <- names(proportion)[1]
  subset <- lec$names[which(lec$membership==mem)]
  subset <- sample(subset, n)
  weight <- as.vector(genes_weight[names(genes_weight) %in% subset])
  return(list(geneID=subset, geneW=weight))
}

betweenness.dereg <- function(pathway = NULL, mode = "vertex", n = 10, 
                              seed=1){
  set.seed(seed)
  adj <- as.matrix(igraph::get.adjacency(pathway, attr = "weight"))
  adj <- colSums(adj)
  genes <- ifelse(adj>=1, 1, ifelse(adj==0, 0, -1))
  genes_weight <- genes[genes != 0]
  genes <- names(genes_weight)
  size <- length(genes)
  
  # genes <- V(pathway)
  # genes <- as_ids(genes)
  # size <- length(genes)
  el <- data.frame(as_edgelist(pathway))
  colnames(el) <- c("src","dest")
  
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  ## calculate edge betweenness based on gr
  if (mode == "vertex"){
    vb <- igraph::betweenness(gr, directed = F)
  }else if (mode == "edge"){
    vb <- igraph::edge_betweenness(gr, directed = F)
  }
  
  vb <- sort(vb, decreasing = TRUE)
  subset <- names(vb[1:n])
  weight <- as.vector(genes_weight[names(genes_weight) %in% subset])
  return(list(geneID=subset, geneW=weight))
}

# nei: neighborhood distance
## Choose the node that has the largest number of neighbors, together with its neighborhood

neighborhood.dereg <- function(pathway = NULL, n=10, nei=2, seed=1){
  set.seed(seed)
  adj <- as.matrix(igraph::get.adjacency(pathway, attr = "weight"))
  adj <- colSums(adj)
  genes <- ifelse(adj>=1, 1, ifelse(adj==0, 0, -1))
  genes_weight <- genes[genes != 0]
  genes <- names(genes_weight)
  size <- length(genes)
  
  # genes <- V(pathway)
  # genes <- as_ids(genes)
  # size <- length(genes)
  el <- data.frame(as_edgelist(pathway))
  colnames(el) <- c("src","dest")
  
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  ## calculate edge betweenness based on gr
  ngbh <- igraph::ego(gr, nei)
  ngbh_len <- sapply(ngbh, length)
  names(ngbh_len) <- V(gr)$name
  id <- which.max(ngbh_len)
  subset <- ngbh[[id]]$name
  
  ngbh_len <- sort(ngbh_len, decreasing = TRUE)
  subset <- names(ngbh_len[1:n])
  weight <- as.vector(genes_weight[names(genes_weight) %in% subset])
  return(list(geneID=subset, geneW=weight))
}

getDEGS <- function(COVIDall, n, j,topology = "betweenness"){
  
  if (topology == "community"){
    e2a <- lapply(seq(1:length(COVIDall)),
                  function(x) community.dereg(COVIDall[[x]], n = n, seed=j))
  } else if(topology == "betweenness"){
    e2a <- lapply(seq(1:length(COVIDall)),
                  function(x) betweenness.dereg(COVIDall[[x]], n = n, seed=j))
  } else if(topology == "neighborhood"){
    e2a <- lapply(seq(1:length(COVIDall)),
                  function(x) neighborhood.dereg(COVIDall[[x]], n = n, seed=j))
  }
  # e2a <- unlist(unique(e2a))
  df <- quiet(lapply(1:length(e2a), function(x) bind_cols(e2a[[x]])))
  df <-ldply(df); colnames(df) <- c("gene", "weight")
  dupl <- subset(df, duplicated(gene))
  dupl <- df %>% dplyr::filter(gene %in% dupl$gene)
  not_dupl <- subset(df,!duplicated(gene)) 
  dist <- dupl %>%
    distinct()
  e2a <- rbind(not_dupl, dist)
  return(e2a)
}

# Create simulation design (choose muval value & entrez2affect) ----------------

sim_design <- function(GSAdata_list, muval = 0.1, e2a = NULL, ncond = 2,
                       seed = 1){
  
  x <- GSAdata_list[["geneSE"]]
  y <- GSAdata_list[["group"]]
  phenodata <- GSAdata_list[["phenodata"]]
  
  # scaling data within group after t()
  x1<- t(scale(t(x[,(phenodata$Sample[phenodata$Group == 1])])))
  x0<- t(scale(t(x[,(phenodata$Sample[phenodata$Group == 0])])))
  z<- cbind(x1,x0)
  
  # set seed and add N(0,1) noise
  set.seed(seed)
  Z <- matrix(0, ncol = ncol(z), nrow = nrow(z))
  for (i in 1:nrow(z)) {
    Z[i,] <- z[i,] + rnorm(n=ncol(z), mean=0, sd=1) #0.6 or 0.1
  }
  colnames(Z)<- colnames(z)
  rownames(Z)<- rownames(z)
  
  # select DEGs and mean signal:
  e2aID <- e2a$gene
  e2aID <- unique(e2aID[e2aID %in% rownames(Z)])
  varZ<- apply(Z, 1, var)
  v<- mean(varZ[e2aID])
  
  # ind.1 <- sample(e2a, length(e2a)*0.5)
  # ind.2 <- e2a[!e2a %in% ind.1]
  # e2aW <- e2a$weight
  ind.1 <- e2a[e2a$weight == 1,]$gene; ind.1 <- ind.1[ind.1 %in% rownames(Z)]
  ind.2 <- e2a[e2a$weight == -1,]$gene; ind.2 <- ind.2[ind.2 %in% rownames(Z)]
  
  Z[ind.1, (phenodata$Sample[phenodata$Group == 1])] <- 
    Z[ind.1,(phenodata$Sample[phenodata$Group == 1])] + muval*v
  
  Z[ind.2, (phenodata$Sample[phenodata$Group == 1])] <- 
    Z[ind.2,(phenodata$Sample[phenodata$Group == 1])] - muval*v
  
  group <- phenodata$Group[match(colnames(Z), phenodata$Sample)]
  sample <- phenodata$Sample[match(colnames(Z), phenodata$Sample)]
  phenodata = data.frame(Group = group, Sample= sample)
  
  y<- group; x<- Z
  
  write.table(x, file="xdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(cbind(colnames(x),y), file="cdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(rownames(x), file="rdat.tab", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  
  exprs.file <- file.path(getwd(), "xdat.tab")
  cdat.file <- file.path(getwd(), "cdat.tab")
  rdat.file <- file.path(getwd(), "rdat.tab")
  
  se <- EnrichmentBrowser::readSE(exprs.file, cdat.file, rdat.file)
  se <- deAna(se)
  
  return(data = list(graph = graph, geneSE = x, se = se, group = y,
                     phenodata = phenodata))
  
}

# Get KEGG (hsa) pathways for all selected GSA methods -------------------------

getALLPATH <- function(GSAdata_list, disease_check = NULL, maxNodes=NULL,
                       minNodes = NULL, maxEdges = NULL,
                       maxComp = NULL){
  
  # mutualEdges = FALSE
  geneSE = GSAdata_list[["geneSE"]]
  
  # Get kegg.gr
  pathways <- kegg.pathways
  N1<-length(pathways)
  pathway_disease <- pathways[which(names(pathways) %in% disease_check)]
  
  if (!is.null(minNodes)) pathways<-SmallPaths(pathways, minNodes)
  if (!is.null(maxNodes)) pathways<-BigPaths(pathways, maxNodes)
  if (!is.null(maxComp)) pathways<-GraphComp(pathways, maxComp)
  # if (isTRUE(mutualEdges)) pathways<-MutualEdges(pathways)
  if (!is.null(maxEdges)) pathways<-DensePaths(pathways, maxEdges)
  
  if(isFALSE(disease_check %in% names(pathways))){
    pathways_igraph <- append(pathways, pathway_disease)
  } else{
    pathways_igraph <- pathways
  }
  
  kegg.gr <- lapply(1:length(pathways_igraph), function(x)
    quiet(properties(pathways_igraph[[x]])[[1]]))
  names(kegg.gr) <- names(pathways_igraph)
  
  graphnel <- lapply(1:length(kegg.gr), 
                          function(x) igraph::as_graphnel(kegg.gr[[x]]))
  names(graphnel) <- names(kegg.gr)
  
  # Get kegg.gs
  # kegg.gs <- lapply(1:length(pathways_igraph), function(x)
  #   as_ids(V(quiet(properties(pathways_igraph[[x]])[[1]]))))
  kegg.gs <- lapply(1:length(kegg.gr), function(x)
    as_ids(V(kegg.gr[[x]])))
  names(kegg.gs) <- names(kegg.gr)
  
  message(paste(length(kegg.gr), "pathways were filtered out"))
  
  return(list(kegg.gr = kegg.gr, graphnel = graphnel, kegg.gs = kegg.gs))
}

MutualEdges<-function (pathways) {
  nmu <- unlist(lapply(1:length(pathways), function(x) sum(which_mutual(pathways[[x]]))))
  blacklist <- which(nmu != 0);
  if (length(blacklist) != 0) pathways <- pathways[-blacklist]
  return(pathways)
}

BigPaths<-function (pathways, maxNodes) {
  pathways <- Filter(function(p) length(V(p)) <= maxNodes, pathways)
  return(pathways)
}

SmallPaths<-function (pathways, minNodes) {
  pathways <- Filter(function(p) length(V(p)) >= minNodes, pathways)
  return(pathways)
}

DensePaths<-function (pathways, maxEdges) {
  pathways <- Filter(function(p) length(E(p)) <= maxEdges, pathways)
  return(pathways)
}

GraphComp<-function (pathways, maxComp) {
  ig <- quiet(lapply(names(pathways), function(x) properties(pathways[[x]])[[1]]))
  nodes <- lapply(names(pathways), function(x) vcount(pathways[[x]]))
  prop <- lapply(seq(1:length(ig)), function (x) vcount(ig[[x]])/nodes[[x]] >= maxComp)
  prop_true <- sapply(prop, any)
  # filter <- lapply(which(prop_true), function (x) ig[[x]])
  filter <- lapply(which(prop_true), function (x) pathways[[x]])
  names <- names(pathways)[which(prop_true)]; names(filter) <- names
  return(filter)
}

# Get TP pathways (COVID related) for simulation -------------------------------

getCOVIDpath <- function(path_all){
  
  covid_related <- c("mRNA surveillance pathway", "Endocytosis",
                     "Vascular smooth muscle contraction", "Complement and coagulation cascades",
                     "Platelet activation", "Neutrophil extracellular trap formation",
                     "Renin-angiotensin system", "Toll-like receptor signaling pathway",
                     "NOD-like receptor signaling pathway", "RIG-I-like receptor signaling pathway",
                     "Cytosolic DNA-sensing pathway", "JAK-STAT signaling pathway",
                     "Natural killer cell mediated cytotoxicity", "Fc gamma R-mediated phagocytosis",
                     "Leukocyte transendothelial migration") #TNF signaling pathway
  
  covid_related <- covid_related[covid_related %in% names(path_all$kegg.gr)]
  COVIDrelated <- path_all$kegg.gr[covid_related]
  # COVIDrelated <- path_all$kegg.gr[which(covid_related %in% names(path_all$kegg.gr))]; 
  names <- names(COVIDrelated)
  COVIDrelated <- lapply(1:length(COVIDrelated), function(x) quiet(properties(COVIDrelated[[x]])[[1]]))
  COVID <- quiet(properties(path_all$kegg.gr[["Coronavirus disease - COVID-19"]])[[1]])
  COVIDall <- append(list(COVID), COVIDrelated)
  names(COVIDall) <- c("Coronavirus disease - COVID-19", 
                       names)
  message(paste(length(COVIDall), "pathways were filtered out"))
  
  return(COVIDall)
} 

# Run SEMgsa -------------------------------------------------------------------

SEMgsa<- function(g=list(), data, group, method="BH",
                  alpha=0.05, n_rep = 1000)
{
  # Set SEM objects:
  gs<- names(g)
  K<- length(g)
  res.tbl<- NULL
  DRN<- list()

  for (k in 1:K){ #k=1
    cat( "k=", k, gs[k], "\n" )
    #ig <- g[[k]]
    ig <- simplify(g[[k]], remove.loops = TRUE)
    err <- paste(" ValueError: none of pathway=", k,
                 " variables are present in the dataset.\n", sep = "")

    # SEM fitting
    fit <- NULL
    tryCatch(SEMgraph:::quiet(fit <- SEMgraph:::SEMricf(
      ig, data, group, random.x = FALSE, n_rep)),
      error = function(c) cat(err))
    if( length(fit[[1]]) == 0 ) {
      res.tbl<- rbind(res.tbl, rep(NA,3))
      next
    }
    p<- ncol(fit$dataXY)
    B<- (diag(p)-fit$fit$ricf$Bhat)[-1,-1] #B
    if( sum(B) == 0 ) {
      res.tbl<- rbind(res.tbl, rep(NA,3))
      next
    }
    pval<- fit$gest$pvalue[-c(1:3)]
    theta <- fit$gest$Stat
    genes<- gsub("X", "", rownames(fit$gest))[-c(1:3)]
    genes<- genes[p.adjust(pval, method=method) < alpha]
    DRN<- c(DRN, list(genes))

    pNa_fit = fit$pval[1]; pNi_fit = fit$pval[2]
    PVAL <- 2*min(pNa_fit,pNi_fit)

    No.DEGs <- sum(p.adjust(pval, method= method) < alpha)
    res.tbl <- rbind(res.tbl, cbind(vcount(ig), No.DEGs,
                                    PVAL))
  }

  res.tbl<- data.frame(res.tbl)
  colnames(res.tbl)<- c("No.nodes","No.DEGs",
                        "PVAL")
  rownames(res.tbl) <- names(g)

  # ADJP <- p.adjust(res.tbl$PVAL, method="BH")
  ADJP <- p.adjust(res.tbl$PVAL, method="bonferroni")
  gsa <- cbind(res.tbl, ADJP)

  # Sorting data frames (ex. by pD)
  # gsa <- gsa[order(gsa$ADJP),]

  return( list(gsa=gsa, DRN=DRN) )
}


run_SEMgsa <- function(GSAdata_list, path_list, n_rep = 2000){
  
  se = GSAdata_list[["se"]]
  data<- t(assay(se))
  group<- se$GROUP
  
  res <- SEMgsa(g=path_list, data, group, method = "BH", 
                alpha=0.05, n_rep=n_rep)
  out <- data.frame(Name = rownames(res[["gsa"]]), 
                    PVAL = res$gsa$ADJP)
  
  return(out)
}


# Run NetGSA -------------------------------------------------------------------

run_NetGSA <- function(GSAdata_list, path_list, alpha=0.05, perm=2000) 
{ 
  se = GSAdata_list[["se"]]
  xx <- assay(se)
  group <- se$GROUP + 1
  K <- length(path_list)
  res.tbl <- NULL
  
  for (k in 1:K){ #k=1
    cat( "k=", k, names(path_list)[k],"\n" )
    ig <- path_list[[k]]
    # quiet(ig<- properties(path_list[[k]])[[1]])
    common <- intersect(names(V(ig)), rownames(xx))
    pp <- length(common) #gene number
    x <- xx[rownames(xx) %in% common,] #gene data for two conditions
    
    A <- as_adjacency_matrix(ig, sparse = FALSE); A <- A[common,common]
    A <- A[order(as.character(rownames(A))),]
    
    # Run NetGSA assuming the network is shared between the two conditions
    B <- matrix(rep(1,pp), nrow=1)
    colnames(B) <- rownames(A)
    rownames(B) <- names(path_list[k])
    x <- x[order(as.character(rownames(x))), ]
    
    lambda <- sqrt(log(pp)/ncol(x))
    rho <- 0.1 * lambda
    eta <- 0.01 * (1+pp/2)
    wAdj <- tryCatch(netgsa::netEst.undir(x=x, one=A, lambda=lambda, 
                                          rho=rho, eta=eta)$Adj, error=function(e) NA)
    fit <- tryCatch(netgsa::NetGSA(A=list(wAdj,wAdj), x=x, group=group,
                                   pathways=B, lklMethod="REHE")$results,
                    error=function(e) NA)		
    if(is.na(wAdj) | is.na(fit)) {
      res.tbl<- rbind(res.tbl, rep(NA,3))
      next
    }
    res.tbl <- rbind(res.tbl, c(fit$pSize, fit$teststat, fit$pval))
  }
  
  res.tbl <- data.frame(res.tbl)
  colnames(res.tbl) <- c("No.genes", "WaldTest", "PVAL")
  rownames(res.tbl) <- names(path_list)
  
  out <- res.tbl %>%
    mutate(Name = rownames(res.tbl),
           PVAL = p.adjust(PVAL, method="BH")) %>%
    dplyr::select(Name, PVAL)
  rownames(out) <- NULL
  
  return(out)
}

# Run DEGraph ------------------------------------------------------------------

#Permutation test of mean difference using ANOVA-type test
#with N(0,1) moment-base approximation (Larson & Owen, 2015)

D.test <- function(X, group, S, n_rep, ...)
{ 
  EE<- eigen(S) # Eigenvalues-Eigenvectors of Psi
  R<- EE$vectors%*%diag(1/sqrt(EE$values))%*%t(EE$vectors)
  D<- as.matrix(X)%*%R%*%rep(1, nrow(S))
  B<- list(B=n_rep+1, seed=123)
  gest<- flip::flip(D, ~group, perms=B, statTest="coeff")
  D<- gest@res$Stat
  # pD<- gest@res$"p-value"
  pp<- sum(gest@permT[-1,] >= D)/n_rep
  pm<- sum(gest@permT[-1,] <= D)/n_rep
  pD<- 2*min(pp, pm)
  # N(0,1) approximate pvalues
  aveT<- mean(gest@permT[-1,])
  sdT<- sd(gest@permT[-1,])
  z<- abs(gest@permT[1,]- aveT)/sdT
  pz<- 2*(1-pnorm(z))
  
  return(c(D, pD, pz))
}

### DEGraph analysis 

run_DEGraph <- function(GSAdata_list, path_list, alpha=0.05, perm=2000) 
{ 
  se = GSAdata_list[["se"]]
  xx <- t(assay(se))
  yy <- se$GROUP
  K <- length(path_list)
  res.tbl <- NULL
  
  for (k in 1:K){ #k=1
    cat( "k=", k, names(path_list)[k],"\n" )
    # Graph of the first (max) connected component
    g1 <- path_list[[k]]
    # quiet(g1<- properties(path_list[[k]])[[1]])
    selx <- which(colnames(xx) %in% V(g1)$name)
    selg <- which(V(g1)$name %in% colnames(xx)[selx])
    g1 <- induced_subgraph(g1, selg)
    #gplot(g1); E(g1)$weight
    if (ecount(g1) < 2) {
      res.tbl<- rbind(res.tbl, c(vcount(g1), 0, NA, NA, NA))
      next
    }
    # Laplacian associated to an weighted (-1,0,1) adjacency matrix
    #ltype=c("meanInfluence", "normalized", "unnormalized", "totalInfluence")
    Adj <- as.matrix(as_adjacency_matrix(g1, attr="weight"))
    lfA <- DEGraph::laplacianFromA(Adj, k=1, ltype="meanInfluence") #str(lfA)
    U <- lfA$U # dim(U)
    #l <- lfA$l; length(l)
    rk <- max(which(lfA$kIdx)) #rk
    # DEGraph covariance (precision) matrix S(C)
    Xd <- xx[,selx]%*%U[,1:rk, drop=FALSE] #dim(Xd)
    if (ncol(Xd) < 2) {
      res.tbl<- rbind(res.tbl, c(vcount(g1), rk, NA, NA, NA))
      next
    }
    if (corpcor::is.positive.definite(cor(Xd), method = "chol")){
      Sd <- cor(Xd)
    }else{
      Sd <- corpcor::cor.shrink(Xd,verbose=TRUE)[1:rk,1:rk]
    }
    # DEgraph D2 statistics
    #D2 <- D2.test(xx = Xd, yy = yy, S = Sd, perm = perm, X2 = TRUE)
    D <- tryCatch(D.test(X = Xd, group = yy, S = Sd, n_rep = perm), 
                  error=function(e) c(NA,NA,NA))
    res.tbl <- rbind(res.tbl, c(vcount(g1), rk, D))
  }
  
  res.tbl<- data.frame(res.tbl) #sum(is.na(res.tbl$PVAL)) #38
  colnames(res.tbl) <- c("No.genes", "No.dimension", "D.statistic", "pD", "PVAL")
  rownames(res.tbl) <- names(path_list)
  
  out <- res.tbl %>%
    mutate(Name = rownames(res.tbl),
           PVAL = p.adjust(PVAL, method="BH")) %>%
    dplyr::select(Name, PVAL)
  rownames(out) <- NULL
  
  return(out)
}


# Run PathwayExpress -----------------------------------------------------------
# run_SPIA <- function(GSAdata_list, path_list, alpha=0.05, perm=2000)
# {
#   # set pathways, DEs with their fcs and entrez names:
#   se = GSAdata_list[["se"]]
#   top <- as.data.frame(rowData(se)) #head(top)
#   entrez <- top$ENTREZID
#   selDE <- which(top$ADJ.PVAL < alpha)
#   if (length(selDE) == 0){
#     out <- data.frame(Name = names(path_list) ,
#                       PVAL = rep(NA, length(path_list)))
#   } else {
#     DE <- top$ENTREZID[selDE]
#     fc <- top$FC[selDE]
#     names(fc) <- entrez[selDE]
#     
#     # compute the number of genes and DEs, and perform SPIA
#     nGS<- unlist(lapply(1:length(path_list), function(x) length(nodes(path_list[[x]]))))
#     nDE<- unlist(lapply(1:length(path_list), function(x) length(intersect(nodes(path_list[[x]]),DE))))
#     # graphnel_list <- lapply(1:length(gr), function(x) igraph::as_graphnel(gr[[x]]))
#     # names(graphnel_list) <- names(gr)
#     spia<- ROntoTools::pe(x=fc, graphs=path_list, ref=top$ENTREZID, nboot=perm, verbose=TRUE)
#     #Summary(spia); colnames(Summary(spia))
#     
#     res.tbl <- data.frame(Name = rownames(Summary(spia)), Summary(spia)[,c(1,6:8)])
#     colnames(res.tbl)[5] <- "PVAL"
#     rownames(res.tbl)<-NULL
#     
#     res.tbl <- dplyr::add_row(res.tbl,Name = names(path_list)[which((names(path_list) %in% res.tbl$Name) == FALSE)], PVAL = NA)
#     #res.tbl<- data.frame(No.genes=nGS, No.DE=nDE, Summary(spia)[,c(1,6:8)])
#     #rownames(res.tbl) <- NULL
#     
#     out <- res.tbl %>%
#       mutate(PVAL = p.adjust(PVAL, method="BH")) %>%
#       dplyr::select(Name, PVAL)
#   }
#   
#   return(out)
# }

run_PathwayExpress <- function(GSAdata_list, path_list, alpha=0.05, perm=2000) 
{ 
  # set pathways, DEs with their fcs and entrez names:
  se = GSAdata_list[["se"]]
  top <- as.data.frame(rowData(se)) #head(top)
  fc <- top$FC
  names(fc) <- top$ENTREZID
  
  # perform PathwayExpress
  pe<- ROntoTools::pe(x=fc, graphs=path_list, ref=top$ENTREZID, nboot=perm, verbose=TRUE)
  
  res.tbl <- data.frame(Name = rownames(Summary(pe)), Summary(pe)[,c(1,6:8)])
  colnames(res.tbl)[5] <- "PVAL"
  rownames(res.tbl)<-NULL
  
  res.tbl <- dplyr::add_row(res.tbl,Name = names(path_list)[which((names(path_list) %in% res.tbl$Name) == FALSE)], PVAL = NA)
  out <- res.tbl %>%
    mutate(PVAL = p.adjust(PVAL, method="BH")) %>%
    dplyr::select(Name, PVAL)
  
  return(out)
}


# Run ORA ----------------------------------------------------------------------

run_ORA <- function(GSAdata_list, path_list, n_rep = 2000){
  
  se = GSAdata_list[["se"]]
  sbea.res <- sbea(method="ora", se=se, gs=path_list,
                   alpha = 0.05, perm = n_rep, padj.method = "none", beta = 1) 
  out <- data.frame(Name = sbea.res[["res.tbl"]]@listData[["GENE.SET"]], 
                    PVAL = sbea.res[["res.tbl"]]@listData[["PVAL"]])
  
  out <- out %>%
    dplyr::add_row(Name = names(path_list)[which((names(path_list) 
                                                  %in% out$Name) == FALSE)],
                   PVAL = NA) %>%
    mutate(PVAL = p.adjust(PVAL, method="BH"))
  
  return(out)
}

# Run topologyGSA --------------------------------------------------------------

run_topologyGSA <- function(GSAdata_list, path_list, perm.num=2000) 
{ 
  se = GSAdata_list[["se"]]
  xx <- t(assay(se))
  yy <- se$GROUP
  K <- length(path_list)
  res.tbl <- NULL
  
  for (k in 1:K){ #k=22
    cat( "k=", k, names(path_list)[k],"\n" )
    # g <- SEMgraph:::quiet(SEMgraph::properties(path_list[[k]]))[[1]]
    g <- path_list[[k]]
    selx <- which(colnames(xx) %in% V(g)$name)
    selg <- which(V(g)$name %in% colnames(xx)[selx])
    sg <- induced_subgraph(g, selg)
    p <- vcount(sg)
    q <- ecount(sg)
    # gplot(g); gplot(sg)
    
    disease = c("Thyroid cancer","Coronavirus disease - COVID-19",
                "Prostate cancer", "Parkinson disease",
                "MAPK signaling pathway", "Protein processing in endoplasmic reticulum",
                "Endocytosis", "Wnt signaling pathway", "Notch signaling pathway",
                "Neurotrophin signaling pathway")
    
    #Select, moralize and triangulate sg:
    if (q > 300 & !names(path_list)[k] %in% disease) {
      res.tbl<- rbind(res.tbl, c(p, q, NA, NA, NA))
      next
    }
    sg <- SEMgraph::graph2dag(sg, xx)
    sg <- igraph::as_graphnel(sg)
    ug <- gRbase::moralize(sg) #plot(ug)
    tug <- gRbase::triangulate(ug) #plot(tug)
    
    #Get cliques of an undirected triangulate graph
    
    clqs <- withTimeout({qpgraph::qpGetCliques(tug, verbose=FALSE)}, 
                        timeout = 1.08, onTimeout = "silent")
    
    if (is.null(clqs)){
      res.tbl<- rbind(res.tbl, c(p, q, NA, NA, NA))
      next
    }
    # clqs <- qpgraph::qpGetCliques(tug, verbose=FALSE)
    # MLE of the regolarized sample covariance matrix using IPF 
    X <- xx[,which(colnames(xx) %in% graph::nodes(sg))]
    if (corpcor::is.positive.definite(cor(X),method = "chol")){
      S <- cor(X)
    }else{
      S <- corpcor::cor.shrink(X,verbose=TRUE) #[1:p,1:p]
    }
    S_ipf <- qpgraph::qpIPF(S, clqs)
    
    # get the adjacency matrix and put the diagonal to one
    #A <- as(ug, "matrix");	diag(A) <- 1
    # entries in S and S_ipf for present edges in g should coincide
    #sum(abs(S_ipf[A==1] - S[A==1])) # 0
    # entries in the inverse of S_ipf for missing edges in g should be zero
    #sum(C[A==0]) #0
    
    #D2 <- D2.test(xx = X, yy = yy, S = S_ipf, perm = perm, X2 = TRUE)
    D <- D.test(X = X, group = yy, S = S_ipf, n_rep = perm.num)
    res.tbl <- rbind(res.tbl, c(p, q, D))
  }
  
  res.tbl <- data.frame(res.tbl)
  colnames(res.tbl) <- c("No.genes", "No.edges", "D.statistic", "pD", "PVAL")
  rownames(res.tbl) <- names(path_list)
  
  out <- res.tbl %>%
    mutate(Name = rownames(res.tbl),
           PVAL = p.adjust(PVAL, method="BH")) %>%
    dplyr::select(Name, PVAL)
  rownames(out) <- NULL
  
  return(out)
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

# RUN all GSA methods ----------------------------------------------------------

getGSAresults <- function(GSAdata_list, path_list){
  
  out = list()
  
  cat('... Running topologyGSA ...', "\n")
  out[["topologyGSA"]] <- quiet(run_topologyGSA(GSAdata_list,
                                                path_list = path_list$kegg.gr))
  
  cat('... Running DEGraph ...', "\n")
  out[["DEGraph"]] <- quiet(run_DEGraph(GSAdata_list,
                                        path_list = path_list$kegg.gr))
  
  
  cat('... Running NetGSA ...', "\n")
  out[["NetGSA"]] <- quiet(run_NetGSA(GSAdata_list,
                                      path_list = path_list$kegg.gr))
  
  cat('... Running ORA ...', "\n")
  out[["ORA"]] <- quiet(run_ORA(GSAdata_list,
                                path_list = path_list$kegg.gs))
  
  cat('... Running SEMgsa ...', "\n")
  out[["SEMgsa"]] <- quiet(run_SEMgsa(GSAdata_list,
                          path_list = path_list$kegg.gr))
  
  cat('... Running PathwayExpress ...', "\n")
  out[["PathwayExpress"]] <- quiet(run_PathwayExpress(GSAdata_list,
                                  path_list = path_list$graphnel))
  
  return(out)
}

# Get p-value of target pathway ------------------------------------------------

getPVALdisease <- function(results, disease=NULL) {
  
  method <- sub(".*\\.", "", colnames(results))
  method <- method[method != "Name"]
  size = nrow(results)
  pval_df = c()
  
  for (i in method){
    
    pval <- results %>%
      dplyr::select(Name, contains(i)) 
    colnames(pval) <- c("Name", "p.value")
    
    pval <- pval %>%
      dplyr::filter(Name == disease) %>% 
      dplyr::mutate(method = i) %>%
      dplyr::select(method, p.value)
    
    pval_df <- rbind(pval_df, pval)
    
  }
  
  pval_df = pval_df %>% arrange(method)
  return(pval_df)
  
}

# Get rank of target pathway ---------------------------------------------------

getRANKdisease <- function(results, disease=NULL) {
  
  method <- sub(".*\\.", "", colnames(results))
  method <- method[method != "Name"]
  size = nrow(results)
  rank_tot = c()
  
  for (i in method){
    
    rank_df <- results %>%
      dplyr::select(Name, contains(i)) %>%
      dplyr::arrange(across(where(is.numeric))) 
    colnames(rank_df) <- c("Name", "PVAL")
    
    # rank_na <- rank_df %>%
    #   dplyr::filter(is.na(PVAL))
    rank_df <- rank_df %>%
      dplyr::filter(!is.na(PVAL))
    if(nrow(rank_df)== 0){
      rank_tot <- rbind(rank_tot, c(i, rep(NA,3)))
    }
    
    PVAL <- rank_df$PVAL
    
    comp.ranks <- vapply(PVAL, function(p) mean(PVAL <= p, na.rm=TRUE) * 100, numeric(1))
    ucats <- unique(PVAL)
    abs.ranks <- match(PVAL, ucats)
    rel.ranks <- abs.ranks / length(ucats) * 100
    
    # rank = base::rank(PVAL, na.last = TRUE,
    #                   ties.method = "min")
    rank <- rank_df %>%
      dplyr::mutate(abs.ranks = abs.ranks,
                    rel.ranks = rel.ranks,
                    comp.ranks = comp.ranks) %>%
      dplyr::filter(Name == disease) %>% 
      dplyr::mutate(method = i) %>%
      dplyr::select(method,  abs.ranks,
                    rel.ranks, comp.ranks)
    
    rank_tot <- rbind(rank_tot, rank)
    
  }
  
  rank_tot = rank_tot %>% arrange(method)
  return(rank_tot)
  
}

# Get DEGs path (also for cross-tak) -------------------------------------------

getTRUEpath <- function(path_all, COVIDall, e2a = NULL){
  
  e2a_path <- lapply(1:length(path_all$kegg.gs), 
                     function (x) sum(e2a %in% (path_all$kegg.gs[[x]])))
  names(e2a_path) <- names(path_all$kegg.gs)
  df_e2a_path <- data.frame(e2a_path=unlist(e2a_path))
  df_e2a_path$Name <- rownames(df_e2a_path); rownames(df_e2a_path) <- NULL
  TP <- df_e2a_path %>% 
    dplyr::mutate(TP = ifelse(e2a_path <= 1, 0,1)) %>%
    dplyr::select(Name, TP) %>% dplyr::arrange(Name) 
  TPsubset <- TP[TP$Name %in% names(COVIDall),]
  TP <- TP %>% 
    dplyr::filter(TP == 0) %>%
    bind_rows(TPsubset)
  
  return(TP = TP)
}

# Get TP,TN,FP,FN --------------------------------------------------------------

getERROR <- function(results, TP = NULL, alpha = 0.05, eps = 1e-06) {
  
  method <- sub(".*\\.", "", colnames(results))
  method <- method[method != "Name"]
  size = nrow(results)
  power_bind = c()
  for (i in method){
    
    res <- results %>%
      dplyr::select(Name, contains(i)) %>%
      dplyr::filter(Name %in% TP$Name) %>%
      dplyr::arrange(across(where(is.numeric))) 
    colnames(res) <- c("Name", "PVAL")
    
    TP_gsa <- res %>%
      dplyr::mutate(TP = ifelse(PVAL <= alpha, 1,0)) %>%
      dplyr::arrange(Name) 
    #pull(TP)
    
    if(sum(is.na(TP_gsa$PVAL)) > 0){
      na<-TP_gsa[is.na(TP_gsa$PVAL),]$Name
      Amat = TP %>% dplyr::filter(!Name %in% na)%>%
        dplyr::arrange(Name) %>% pull(TP)
      Ahat = TP_gsa %>% dplyr::filter(!Name %in% na)%>%
        dplyr::arrange(Name) %>% pull(TP)
    } else{
      Ahat <- TP_gsa %>% pull(TP)
      Amat = TP %>% dplyr::arrange(Name) %>% pull(TP)
    }
    
    
    power_res = data.frame(method = i) %>%
      dplyr::mutate(FP = round(sum((abs(Ahat) > eps) * (abs(Amat) <= eps), na.rm = T)/sum(Amat==0),4),
                    FN = round(sum((abs(Ahat) <= eps) * (abs(Amat) > eps), na.rm = T)/sum(Amat),4), 
                    power = round((1-FN),4))
    
    # TP = sum((abs(Ahat) > eps) * (abs(Amat) > eps), na.rm = T)/sum(Amat),
    # TN = sum((abs(Ahat) <= eps) * (abs(Amat) <= eps), na.rm = T)/sum(Amat==0)
    power_bind <- rbind(power_bind, power_res)
    
  }
  
  return(power_bind)
}        

getPRIORITIZATION <- function(results, TP = NULL) {
  
  method <- sub(".*\\.", "", colnames(results))
  method <- method[method != "Name"]
  size = nrow(results)
  prior_tot = c()
  
  for (i in method){
    
    rank_df <- results %>%
      dplyr::select(Name, contains(i)) %>%
      dplyr::filter(Name %in% TP$Name) %>%
      dplyr::arrange(across(where(is.numeric))) 
    colnames(rank_df) <- c("Name", "PVAL")
    
    rank_df <- rank_df %>%
      dplyr::filter(!is.na(PVAL)) %>%
      arrange(PVAL)
    TP <- TP %>%
      dplyr::filter(Name %in% rank_df$Name) %>%
      arrange(Name)
    
    rank_df <- rank_df %>%
      dplyr::mutate(rank = seq(1:nrow(rank_df))) %>% 
      dplyr::arrange(Name) %>%
      dplyr::mutate(TP = TP$TP)
    
    # rank_df = rank_df %>%
    #   mutate(rank_new = ifelse(rank <= sum(TP) & TP == 1
    #                            & PVAL <= 0.05, 1,0))
    
    rank_df = rank_df %>%
      mutate(rank_new = ifelse(rank <= sum(TP) & TP == 1, 1,0))
    
    prioritization = sum(rank_df$rank_new == rank_df$TP)/nrow(rank_df)
    prior = data.frame(method = i, prioritization = prioritization)
    prior_tot <- rbind(prior_tot, prior)
    
  }
  
  prior_tot = prior_tot %>% arrange(method)
  return(prior_tot)
}


# Get summary of GSA results (target) ------------------------------------------

getSUMMARY <- function(out_list, disease=NULL, sim = FALSE, TP = NULL, alpha = 0.05){
  
  out_long <- as.data.frame(do.call(rbind, lapply(out_list, as.vector)))
  out_long <- cbind(method=sub("\\..*", "", rownames(out_long)), out_long)
  rownames(out_long) <- NULL
  
  out_wide <- reshape(out_long, idvar = "Name", timevar = "method", direction = "wide")
  out_wide <- out_wide %>% 
    dplyr::select(where(~mean(is.na(.))< 0.5))
  
  if (sim == FALSE){
    
    rank <- getRANKdisease(out_wide, disease)
    pval <- getPVALdisease(out_wide, disease)
    
    results <- pval %>%
      dplyr::mutate(abs.ranks = ifelse(is.na(p.value), NA, rank$abs.ranks),
                    rel.ranks = ifelse(is.na(p.value), NA, rank$rel.ranks),
                    comp.ranks = ifelse(is.na(p.value), NA, rank$comp.ranks))
    
    out <- list(out_wide = out_wide, metrics = results)
    
  } else{
    
    results1 <- getERROR(out_wide, TP, alpha)
    results2 <- getPRIORITIZATION(out_wide, TP)
    
    out <- list(out_wide = out_wide, metrics = results1, prioritization = results2)
    
  }
  
  return(out)
}

# Simulation metrics -----------------------------------------------------------

getMEDIANERROR <- function(GSA_sim_out) {
  
  method <- sub(".*\\.", "", colnames(GSA_sim_out[[1]][["out_wide"]]))
  method <- method[method != "Name"]
  res_bind = c()
  for (i in method){
    for (j in (1:length(GSA_sim_out))){
      
      res <- GSA_sim_out[[j]][["metrics"]]
      res <- res %>%
        dplyr::filter(method == i) %>%
        dplyr::mutate(N.iter = j) %>%
        dplyr::select(N.iter, method, TP, TN, FP, FN) 
      
      res_bind <- rbind(res_bind, res)
      
    }
  } 
  
  error <- res_bind %>%
    dplyr::group_by(method) %>%
    dplyr::summarise_at(c("TP","TN","FP","FN"), median, na.rm = TRUE) %>%
    arrange(method)
  
  
  return(error)
}

# Boxplot ----------------------------------------------------------------------

GSAplot <- function(GSA_sim_list, y="power", type = "errorplot"){
  # x => method
  # y => power/FP
  
  method <- sub(".*\\.", "", colnames(GSA_sim_list[[1]][["out_wide"]]))
  method <- method[method != "Name"]
  res_bind = c()
  for (i in method){
    for (j in (1:length(GSA_sim_list))){
      
      if(y == "prioritization"){
        res <- GSA_sim_list[[j]][["prioritization"]]
        res <- res %>%
          dplyr::filter(method == i) %>%
          dplyr::mutate(N.iter = j) %>%
          dplyr::select(N.iter, method, prioritization) 
      } else{
        res <- GSA_sim_list[[j]][["metrics"]]
        res <- res %>%
          dplyr::filter(method == i) %>%
          dplyr::mutate(N.iter = j) %>%
          dplyr::select(N.iter, method, FP, power) 
      }
      
      res_bind <- rbind(res_bind, res)
      
    }
  } 
  
  # res_bind = res_bind %>%
  #   dplyr::mutate(type=ifelse(method=="SEMgsa_bonf" | method== "SEMgsa_brown",
  #                             "Highlighted","Normal")) %>%
  #   arrange(method)
  
  res_bind = res_bind %>%
    dplyr::mutate(type=ifelse(method=="SEMgsa",
                              "Highlighted","Normal")) %>%
    arrange(method)
  
  del <- c("ORA")
  res_bind = res_bind %>%
    dplyr::filter(!method %in% del) 
    # dplyr::mutate(method= recode(method, "SEMgsa_brown"="SEMgsa"))
  
  if(type == "boxplot"){
    if(y == "FP"){
      ggplot(res_bind, aes(x=method, y=FP, fill=type, alpha=type)) + 
        geom_boxplot() +
        scale_fill_manual(values=c("#69b3a2", "grey")) +
        scale_alpha_manual(values=c(1,0.1)) +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(text = element_text(size = 19))  +
        xlab("") 
    } else if (y == "power"){
      ggplot(res_bind, aes(x=method, y=power, fill=type, alpha=type)) + 
        geom_boxplot() +
        scale_fill_manual(values=c("#69b3a2", "grey")) +
        scale_alpha_manual(values=c(1,0.1)) +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(text = element_text(size = 19))  +
        xlab("") + ylab("Power")
    } else if (y == "prioritization"){
      ggplot(res_bind, aes(x=method, y=prioritization, fill=type, alpha=type)) + 
        geom_boxplot() +
        scale_fill_manual(values=c("#69b3a2", "grey")) +
        scale_alpha_manual(values=c(1,0.1)) +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(text = element_text(size = 19))  +
        xlab("") + ylab("Prioritization")
    }
  } else if(type == "errorplot"){
    if(y == "FP"){
      ggerrorplot(res_bind, x = "method", y = "FP", color =  "type",
                  palette = "npg",
                  desc_stat = "mean_sd",
                  error.plot = "errorbar",            
                  add = "mean", ggtheme = theme_minimal()) + 
        theme(legend.position = "none")  +
        theme(text = element_text(size = 19))  +
        xlab("") 
    } else if (y == "power"){
      ggerrorplot(res_bind, x = "method", y = "power", color =  "type",
                  palette = "npg", 
                  desc_stat = "mean_sd",
                  error.plot = "errorbar",            
                  add = "mean", ggtheme = theme_minimal()) + 
        theme(legend.position = "none")  +
        theme(text = element_text(size = 19)) +
        xlab("") + ylab("Power")
    } else if (y == "prioritization"){
      ggerrorplot(res_bind, x = "method", y = "prioritization", color =  "type",
                  palette = "npg", 
                  desc_stat = "mean_sd",
                  error.plot = "errorbar",            
                  add = "mean", ggtheme = theme_minimal()) + 
        theme(legend.position = "none")  +
        theme(text = element_text(size = 19)) +
        xlab("") + ylab("Prioritization")
    }
  }
  
}

