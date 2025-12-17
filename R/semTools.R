#  SEMgraph library
#  Copyright (C) 2019-2021 Mario Grassi; Fernando Palluzzi; Barbara Tarantino 
#  e-mail: <mario.grassi@unipv.it>
#  University of Pavia, Department of Brain and Behavioral Sciences
#  Via Bassi 21, 27100 Pavia, Italy

#  SEMgraph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  SEMgraph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #
#' @title SEM-based gene set analysis
#'
#' @description Gene Set Analysis (GSA) via self-contained test for group
#' effect on signaling (directed) pathways based on SEM. The core of the
#' methodology is implemented in the RICF algorithm of \code{SEMrun()},
#' recovering from RICF output node-specific group effect p-values, and
#' Brown’s combined permutation p-values of node activation and inhibition.
#'
#' @param g A list of pathways to be tested.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, method is set
#' to "BH" (i.e., Benjamini-Hochberg correction).
#' @param alpha Gene set test significance level (default = 0.05).
#' @param n_rep Number of randomization replicates (default = 1000).
#' @param ... Currently ignored.
#'
#' @details For gaining more biological insights into the functional roles
#' of pre-defined subsets of genes, node perturbation obtained from RICF
#' fitting has been combined with up- or down-regulation of genes from a
#' reference interactome to obtain overall pathway perturbation as follows: 
#' \itemize{
#' \item The node perturbation is defined as activated when the minimum among
#' the p-values is positive; if negative, the status is inhibited. 
#' \item Up- or down- regulation of genes is computed from the weighted adjacency
#' matrix of each pathway as column sum of weights(-1,0,1) over each source node.
#' If the overall sum of node weights is below 1, the pathway is flagged as
#' down-regulated, otherwise as up-regulated. 
#' \item The combination between these two quantities allows to define the
#' direction (up or down) of gene perturbation. Up- or down regulated gene status,
#' associated with node inhibition, indicates a decrease in activation (or
#' increase in inhibition) in cases with respect to control group. Conversely,
#' up- or down regulated gene status, associated with node activation, indicates
#' an increase in activation (or decrease in inhibition) in cases with
#' respect to control group.
#' }
#' 
#' @return A list of 2 objects:
#' \enumerate{
#' \item "gsa", A data.frame reporting the following information for each
#' pathway in the input list:
#' \itemize{
#' \item "No.nodes", pathway size (number of nodes);
#' \item "No.DEGs", number of differential espression genes (DEGs) within
#' the pathway, after multiple test correction with one of the methods
#' available in \code{\link[stats]{p.adjust}};
#' \item "pert", pathway perturbation status (see details);
#' \item "pNA", Brown's combined P-value of pathway node activation;
#' \item "pNI", Brown's combined P-value of pathway node inhibition;
#' \item "PVAL", Bonferroni combined P-value of pNA, and pNI; i.e.,
#' 2* min(pNA, PNI);
#' \item "ADJP", Adjusted Bonferroni P-value of pathway perturbation; i.e.,
#' min(No.pathways * PVAL; 1).
#' }
#' \item "DEG", a list with DEGs names per pathways.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Grassi, M., Tarantino, B. (2022). SEMgsa: topology-based pathway enrichment analysis
#' with structural equation models. BMC Bioinformatics, 17 Aug, 23, 344.
#' <https://doi.org/10.1186/s12859-022-04884-8>
#' 
#' @examples
#'
#' \dontrun{
#'
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
#'
#' # Selection of FTD-ALS pathways from KEGG pathways
#'
#' paths.name <- c("MAPK signaling pathway",
#'                 "Protein processing in endoplasmic reticulum",
#'                 "Endocytosis",
#'                 "Wnt signaling pathway",
#'                 "Neurotrophin signaling pathway",
#'                 "Amyotrophic lateral sclerosis")
#' 
#' j <- which(names(kegg.pathways) %in% paths.name)
#'
#' GSA <- SEMgsa(kegg.pathways[j], als.npn, alsData$group,
#'               method = "bonferroni", alpha = 0.05,
#'               n_rep = 1000)
#' GSA$gsa
#' GSA$DEG
#' 
#' }
#'
SEMgsa<- function(g=list(), data, group, method = "BH", alpha = 0.05, n_rep = 1000, ...)
{
 	# set loop objects:
	gs <- names(g)
    K <- length(g)
	res.tbl <- NULL
	DEG <- list()
		
	for (k in 1:K){
	 cat( "k =", k, gs[k], "\n" )
	 ig <- simplify(g[[k]], remove.loops = TRUE)
	 if (length(E(ig)$weight) == 0) E(ig)$weight <- 1
	 adj <- as.matrix(get.adjacency(ig, attr = "weight"))
	 adj <- colSums(adj)
	 nodes <- ifelse(adj >= 1, 1, ifelse(adj == 0, 0, -1))
	 status <- ifelse(sum(nodes) >= 1, 1, ifelse(sum(nodes) == 0, 0, -1))
	 
	 # RICF fitting:
	 fit <- NULL
	 err <- paste(" ValueError: Model converged = FALSE for k =",k,"\n")
	 tryCatch(fit <- quiet(SEMricf(ig, data, group, n_rep)),
	 					    error = function(c) cat(err))
	 if (length(fit[[1]]) == 0) {
	  res.tbl<- rbind(res.tbl, rep(NA,6))
	  colnames(res.tbl) <- c("No.nodes","No.DEGs","pert","pNa","pNi","PVAL")
	  DEG <- c(DEG, list(NULL))
	  next
	 }
	 pval<- fit$gest$pvalue
	 genes<- gsub("X", "", rownames(fit$gest))
	 genes<- genes[p.adjust(pval, method=method) < alpha]
	 DEG<- c(DEG, list(genes))
	 
	 # data.frame of combined SEM results :
	 cpval <- fit$fit$pval
	 names(cpval) <- c("pNa", "pNi")
	 pvmin <- names(cpval)[which.min(cpval)]
	 sign <- ifelse(pvmin == "pNa", "+", "-")
	 if ( status != 0 ) {
	  if (status == -1 & sign == "-") pert <- "up inh"
	  else if (status == -1 & sign == "+") pert <- "down inh"
	  else if (status == 1 & sign == "+") pert <- "up act"
	  else if (status == 1 & sign == "-") pert <- "down act"
	 }else{pert <- NA}
	 df <- data.frame(
	   No.nodes = vcount(ig),
	   No.DEGs = sum(p.adjust(pval, method=method) < alpha),
	   pert = pert,
	   pNa = cpval[1],
	   pNi = cpval[2],
	   PVAL = min(2 * min(cpval[1], cpval[2]), 1)
	 )
	 res.tbl <- rbind(res.tbl,df)
	}
 
 	ADJP <- p.adjust(res.tbl$PVAL, method="bonferroni")
	res.tbl<- cbind(res.tbl, ADJP)
	rownames(res.tbl) <- names(DEG) <- names(g)
	gsa <- res.tbl[order(res.tbl$ADJP),]
	DEG <- DEG[rownames(gsa)]
	
	return( list(gsa=gsa, DEG=DEG) )
}

#' @title SEM-based differential network analysis
#'
#' @description Differential Connected Inference (DCI) via a sub-network with
#' perturbed edges obtained from the output of \code{\link[SEMgraph]{SEMace}},
#' comparable to the procedure in Jablonski et al (2022), or \code{\link[SEMgraph]{SEMrun}}
#' with two-group and CGGM solver, comparable to the algorithm 2 in Belyaeva et al (2021). 
#' To increase the efficiency of computations for large graphs, users can
#' select to break the network structure into clusters, and select the
#' topological clustering method (see \code{\link[SEMgraph]{clusterGraph}}).
#' The function \code{\link[SEMgraph]{SEMrun}} is applied iteratively on
#' each cluster (with size min > 10 and max < 500) to obtain the graph
#' with the full list of perturbed edges. 
#' 
#' @param graph Input network as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#' @param type  Average Causal Effect (ACE) with two-group, "parents"
#' (back-door) adjustement set, and "direct" effects (\code{type = "ace"},
#' default), or CGGM solver with two-group using a clustering method.
#' If \code{type = "tahc"}, network modules are generated using the tree
#' agglomerative hierarchical clustering method, or non-tree clustering
#' methods from igraph package, i.e., \code{type = "wtc"} (walktrap community
#' structure with short random walks), \code{type ="ebc"} (edge betweeness
#' clustering), \code{type = "fgc"} (fast greedy method), \code{type = "lbc"}
#' (label propagation method), \code{type = "lec"} (leading eigenvector method),
#' \code{type = "loc"} (multi-level optimization), \code{type = "opc"} (optimal
#' community structure), \code{type = "sgc"} (spinglass statistical mechanics),
#' \code{type = "none"} (no breaking network structure into clusters).
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, method is set
#' to "BH" (i.e., FDR multiple test correction).
#' @param alpha Significance level (default = 0.05) for edge set selection.
#' @param ... Currently ignored.
#'
#' @return An igraph object.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Belyaeva A, Squires C, Uhler C (2021). DCI: learning causal differences
#' between gene regulatory networks. Bioinformatics, 37(18): 3067–3069.
#' <https://doi: 10.1093/bioinformatics/btab167>
#'
#' Jablonski K, Pirkl M, Ćevid D, Bühlmann P, Beerenwinkel N (2022).
#' Identifying cancer pathway dysregulations using differential
#' causal effects. Bioinformatics, 38(6):1550–1559.
#' <https://doi.org/10.1093/bioinformatics/btab847>
#'
#' @examples
#'
#' \dontrun{
#'
#' #load SEMdata package for ALS data with 17K genes:
#' #devtools::install_github("fernandoPalluzzi/SEMdata")
#' #library(SEMdata)
#' 
#' # Nonparanormal(npn) transformation
#' data.npn<- transformData(alsData$exprs)$data
#' dim(data.npn)
#' 
#' # Selection of FTD-ALS pathways from KEGG pathways
#' 
#' paths.name <- c("MAPK signaling pathway",
#'                 "Protein processing in endoplasmic reticulum",
#'                 "Endocytosis",
#'                 "Wnt signaling pathway",
#'                 "Neurotrophin signaling pathway",
#'                 "Amyotrophic lateral sclerosis")
#' 
#' j <- which(names(kegg.pathways) %in% paths.name)
#' 
#' # Neuro interactome (max component)
#' gU <- Reduce(igraph::union, kegg.pathways[j])
#' gU <- properties(gU)[[1]]
#' #summary(gU)
#' 
#' # Modules with ALS perturbed edges using fast gready clustering
#' gD<- SEMdci(gU, data.npn, alsData$group, type="fgc")
#' summary(gD)
#' gcD<- properties(gD)
#' 
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(2,4))
#' gplot(gcD[[1]], l="fdp", main="max component")
#' gplot(gcD[[2]], l="fdp", main="2nd component")
#' gplot(gcD[[3]], l="fdp", main="3rd component")
#' gplot(gcD[[4]], l="fdp", main="4th component")
#' par(old.par)
#' 
#' }
#'
SEMdci<- function (graph, data, group, type = "ace", method = "BH", alpha = 0.05, ...) 
{
    if (type == "ace") {
        dest <- SEMace(graph, data, group, type = "parents", 
            effect = "direct", method = method, alpha = alpha, 
            boot = NULL)
        ftm <- data.frame(from = dest$source, to = dest$sink)
        return(gD = graph_from_data_frame(ftm))
    }
    if (type != "none") {
        C <- clusterGraph(graph, type = type, size = 10)
        K <- as.numeric(names(table(C)))
        gL <- NULL
        for (k in 1:length(K)) {
            cat("fit cluster =", K[k], "\n")
            V <- names(C)[C == K[k]][names(C)[C == K[k]] %in% colnames(data)]
            g <- simplify(induced_subgraph(graph, vids = V))
            if (vcount(g) > 500 | ecount(g) < 9) next
			dest <- tryCatch(quiet(SEMggm2(g, data, group)$dest),
                             error = function(err) NULL)
			if (is.null(dest)) next
			dsub <- subset(dest, p.adjust(dest$pvalue, method = method) < alpha)
            if (nrow(dsub) == 0) next
            ftm <- data.frame(from = dsub$rhs, to = dsub$lhs)
            gC <- graph_from_data_frame(ftm)
            if (ecount(gC) > 0) gL <- c(gL, list(gC))
        }
        cat("Done.\n")
        if (is.null(gL)) 
            return(gD = make_empty_graph(n = length(K)))
        #gD <- graph.union(gL)
		gD <- Reduce(union, gL)
    }
    else if (type == "none") {
        dest <- quiet(SEMggm2(g, data, group)$dest)
        dsub <- subset(dest, p.adjust(dest$pvalue, method = method) < alpha)
        ftm <- data.frame(from = dsub$rhs, to = dsub$lhs)
        gD <- graph_from_data_frame(ftm)
    }
    return(gD)
}

#' @title Graph weighting methods
#'
#' @description Add data-driven edge and node weights to the input graph.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects. By default, \code{group = NULL}. If group is not NULL, also
#' node weighting is actived, and node weights correspond to the attribute:
#' V(graph)$pv (P-value of the z-test = b/SE(b) from simple linear regression
#' y ~ x, i.e., lm(node ~ group)) and V(graph)$sign (-1 if z<-2, +1 if z>2,
#' 0 otherwise).
#' @param method Edge weighting method. It can be one of the following:
#' \enumerate{
#' \item "r2z", weight edges are defined using Fisher's r-to-z transform
#' (Fisher, 1915) to test the correlation coefficient of pairs of interacting
#' nodes, if \code{group=NULL}. Otherwise, the difference between group of
#' the r-to-z trasform will be tested. Edge weights correspond to the attribute:
#' E(graph)$pv (P-value of the z-test) and E(graph)$sign (-1 if z<-2, +1 if z>2,
#' 0 otherwise).
#' \item "sem", edge weights are defined by a SEM model that implies 
#' testing the group effect simultaneously on source and sink nodes.
#' A new parameter w is defined as the weighted sum of the total effect 
#' of the group on source and sink nodes, adjusted by node degree centrality. 
#' Edge weights correspond to the attribute: E(graph)$pv (P-value of the
#' z-test = w/SE(w)) and E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise). 
#' Not available if \code{group=NULL}.
#' \item "cov", edge weights are defined by a new parameter w combining 
#' the group effect on the source node (mean group difference, adjusted 
#' by source degree centrality), the sink node (mean group difference, 
#' adjusted by sink degree centrality), and the source--sink interaction 
#' (correlation difference). Edge weights correspond to the attribute:
#' E(graph)$pv (P-value of the z-test = w/SE(w) of the combined difference
#' of the group over source node, sink node, and their connection) and
#' E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise).
#' Not available if \code{group=NULL}.
#' \item "cfa", edge weights are defined by a CFA1 model that implies 
#' testing the group effect, w on a latent variable (LV) with observed
#' indicators two interacting nodes, fixing loading coefficients and residual
#' variances for model identification. Edge weights correspond to the
#' attribute: E(graph)$pv (P-value of the z-test = w/SE(w) of the group
#' effect on the LV) and E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise).
#' Not available if \code{group=NULL}.
#' }
#' @param limit An integer value corresponding to the number of graph 
#' edges. Beyond this limit, multicore computation is enabled to reduce 
#' the computational burden. By default, \code{limit = 10000}.
#' @param ... Currently ignored.
#'
#' @return A weighted graph, as an igraph object.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Grassi M, Tarantino B (2023). [Supplementary material of] SEMtree: tree-based structure
#' learning methods with structural equation models. 
#' Bioinformatics, 39 (6), 4829–4830 <https://doi.org/10.1093/bioinformatics/btad377>
#'
#' Fisher RA (1915). Frequency Distribution of the Values of the Correlation
#' Coefficient in Samples from an Indefinitely Large Population. Biometrika,
#' 10(4), 507–521. <doi:10.2307/2331838>
#' 
#' @examples
#' 
#' # Graph weighting
#' G <- weightGraph(graph = sachs$graph,
#'                  data = log(sachs$pkc),
#'                  group = sachs$group,
#'                  method = "r2z")
#' 
#' # New edge attributes
#' head(E(G)$pv); summary(E(G)$pv)
#' head(E(G)$zsign); table(E(G)$zsign)
#' 
#' # New node attributes
#' head(V(G)$pv); summary(V(G)$pv)
#' head(V(G)$zsign); table(V(G)$zsign)
#' 
weightGraph<- function(graph, data, group = NULL, method = "r2z", limit = 10000, ...) 
{
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	ig <- quiet(properties(ig)[[1]])
	degree <- igraph::degree(ig, v = V(ig), mode = "all")
	ftm <- igraph::as_data_frame(ig)
	Y <- scale(data[, nodes])
	if (is.null(group) | method == "r2z")
		ew <- ew.r2z(ftm, Y, group)
	if (method == "sem")
		ew <- ew.sem(ftm, Y, group, degree, limit = limit)
	if (method == "cov") 
		ew <- ew.cov(ftm, Y, group, degree, limit = limit)
	if (method == "cfa")
		ew <- ew.cfa(ftm, Y, group, limit = limit)
	if (method == "lmi") 
		ew <- ew.lmi(ftm, Y, group, limit = limit)
	zsign <- ew[[1]]
	pv <- ew[[2]]
	pv[is.na(pv)] <- runif(n = sum(is.na(pv)), min = 0.05, max = 0.95)
	pv[pv <= 0] <- 1e-15
	pv[pv >= 1] <- 1 - 1e-15
	ftm <- cbind(ftm, zsign, pv)
	gdf <- graph_from_data_frame(ftm, directed = is_directed(graph))
	Vattr <- vertex_attr(graph)
	if (length(Vattr) > 1) {
		idx <- match(V(gdf)$name, Vattr$name)
		for (i in 2:length(Vattr)) gdf <- set_vertex_attr(gdf, 
			names(Vattr)[i], value = Vattr[[i]][idx])
	}
	if (!is.null(group)) {
		graph <- vw.lm(gdf, data[, nodes], group)
	}
	else {
		graph <- gdf
	}
	return(graph)
}

vw.lm<- function(graph, data, group, ...)
{
	est<- lapply(1:ncol(data), function(x) lm(data[,x] ~ group))
	B<- sapply(1:length(est), function(x) summary(est[[x]])$coefficients[2,3])
	zsign<- ifelse(abs(B) < 2, 0, sign(B))
	pv<- sapply(1:length(est), function(x) summary(est[[x]])$coefficients[2,4])
	names(zsign)<- names(pv)<- colnames(data)
	V(graph)$pv<- pv[V(graph)$name]
	V(graph)$zsign<- zsign[V(graph)$name]
	return(graph)
}

ew.sem <- function(ftm, Y, group, degree, limit, ...)
{	
	local <- function(x) {
		df <- data.frame(cbind(Y[, c(x[[1]], x[[2]])], group))
		dx <- degree[x[[1]]]
		dy <- degree[x[[2]]]
		colnames(df)[1:2] <- c("x", "y")
		model <- paste0(
		 'y~ b0*1+b1*group
		 x~ a0*1+a1*group
		 w:=a1/',dx,' + b1/',dy)
		#cat(model)
		try(fit <- sem(model, data = df, fixed.x = TRUE))
		try(res <- parameterEstimates(fit))
		try(res[res$label == "w", -c(1:4)])
	}
	
	x <- split(ftm, f = seq(nrow(ftm)))
	message("Edge weighting via SEM of ", length(x), " edges ...")
	op <- pbapply::pboptions(type = "timer", style = 2)
	if (length(x) > limit) {
		n_cores <- parallel::detectCores(logical = FALSE)
		cl <- parallel::makeCluster(n_cores)
		parallel::clusterExport(cl, c("local", "Y", "degree", "group"),
		                        envir = environment())
		est <- pbapply::pblapply(x, local, cl = cl)
		parallel::stopCluster(cl)
	} else {
		est <- pbapply::pblapply(x, local, cl = NULL)
	}
	
	B <- sapply(1:length(est), function(x) est[[x]]$z)
	zsign <- ifelse(abs(B) < 1.96, 0, sign(B))
	pv <- sapply(1:length(est), function(x) est[[x]]$pvalue)
	cat("\n")
	return (list(zsign, pv))
}

ew.cov <- function(ftm, Y, group, degree, limit, ...)
{
	local <- function(x) {
		df <- data.frame(cbind(Y[, c(x[[1]], x[[2]])], group))
		dx <- degree[x[[1]]]
		dy <- degree[x[[2]]]
		colnames(df)[1:2] <- c("x", "y")
		model <- paste0(
		 'x ~ c(a1,a2)*1
		 y ~ c(b1,b2)*1
		 x ~~ c(c1,c2)*y
		 w:= (a2-a1)/', dx, '+(b2-b1)/', dy, '+(c2-c1)')
		#cat(model)
		try(fit <- sem(model, data = df, group = "group",
		    fixed.x = TRUE))
		try(res<- parameterEstimates(fit))
		try(res[res$label == "w", -c(1:5)])
	}
	
	x <- split(ftm, f = seq(nrow(ftm)))
	message("Edge weighting via COV of ", length(x), " edges ...")
	op <- pbapply::pboptions(type = "timer", style = 2)
	
	if (length(x) > limit) {
		n_cores <- parallel::detectCores(logical = FALSE)
		cl <- parallel::makeCluster(n_cores)
		parallel::clusterExport(cl, c("local", "Y", "degree", "group"),
		                        envir = environment())
		est<- pbapply::pblapply(x, local, cl = cl)
		parallel::stopCluster(cl)
	} else {
		est <- pbapply::pblapply(x, local, cl = NULL)
	}
	
	B <- sapply(1:length(est), function(x) est[[x]]$z)
	zsign <- ifelse(abs(B) < 1.96, 0, sign(B))
	pv <- sapply(1:length(est), function(x) est[[x]]$pvalue)
	cat("\n")
	return (list(zsign, pv))
}

ew.cfa <- function(ftm, Y, group, limit, ...)
{
	local <- function(x) {
		df <- data.frame(cbind(Y[, c(x[[1]], x[[2]])], group))
		colnames(df)[1:2] <- c("y1", "y2")
		if(cov(df$y1, df$y2) <0) df$y1 <- -1*df$y1
		a <- sqrt(cov(df$y1, df$y2))
		model <- paste0(
		 'f =~ ',a,'*y1+',a,'*y2
		 f ~ group
		 y1~~',1-a^2,'*y1
		 y2~~',1-a^2,'*y2')
		#cat(model)
		suppressWarnings(try(fit <- cfa(model, data = df,
		                     fixed.x = TRUE)))
		try(res<- parameterEstimates(fit))
		try(res[c(3, 6),])
	}
	
	x <- split(ftm, f = seq(nrow(ftm)))
	message("Edge weighting via 1CFA of ", length(x), " edges ...")
	op <- pbapply::pboptions(type = "timer", style = 2)
	
	if (length(x) > limit) {
		n_cores <- parallel::detectCores(logical = FALSE)
		cl <- parallel::makeCluster(n_cores)
		parallel::clusterExport(cl, c("local", "Y", "degree", "group"),
		                        envir = environment())
		est <- pbapply::pblapply(x, local, cl = cl)
		parallel::stopCluster(cl)
	} else {
		est <- pbapply::pblapply(x, local, cl = NULL)
	}
	
	var <- sapply(1:length(est), function(x) est[[x]]$est[2])
	B <- sapply(1:length(est), function(x) est[[x]]$z[1])
	zsign <- ifelse(abs(B) < 1.96, 0, sign(B))
	pv <- sapply(1:length(est), function(x) est[[x]]$pvalue[1])
	if (sum(var < 0) > 0) {
		pv[var < 0] <- NA
		message("WARNING ", sum(var < 0), " of ", nrow(ftm),
		" estimated residual var(LV) are negatives")
	}
	return (list(zsign, pv))
}

ew.lmi<- function(ftm, Y, group, limit, ...)
{
	local<- function(x) {
	 df<- data.frame(cbind(Y[,c(x[[1]],x[[2]])],group))
	 colnames(df)[1:2]<- c("x", "y")
	 try(fit<- lm(y ~ x*group, data= df))
	 try(res<- summary(fit)$coefficients)
	 try(res[4,])
	}
	
	x<- split(ftm, f = seq(nrow(ftm)))
	message("Edge weigthing via lm() of ", length(x), " edges...")
	op<- pbapply::pboptions(type = "timer", style = 2)
	if (length(x) > limit){
	 n_cores <- parallel::detectCores(logical = FALSE)
	 cl<- parallel::makeCluster(n_cores)
	 parallel::clusterExport(cl,
 	  c("local", "Y", "degree", "group"), envir = environment())
	 est<- pbapply::pblapply(x, local, cl=cl)
	 parallel::stopCluster(cl)
	}else{
	 est<- pbapply::pblapply(x, local, cl=NULL)
	}		
	
	B<- sapply(1:length(est), function(x) est[[x]][3])
	zsign<- ifelse(abs(B) < 1.96, 0, sign(B))
	pv<- sapply(1:length(est), function(x) est[[x]][4])
	cat("\n")
	return ( list(zsign, pv) )
}

ew.r2z <- function(ftm, Y, group, ...)
{
	zsign <- vector()
	pv <- vector()
	if(!is.null(group)) {
		n1 <- length(group[group == 1])
		n0 <- length(group[group == 0])
		for(i in 1:nrow(ftm)) {
			x <- Y[, ftm[i, 1]]
			y <- Y[, ftm[i, 2]]
			x1 <- Y[group == 1, ftm[i, 1]]
			y1 <- Y[group == 1, ftm[i, 2]]
			x0 <- Y[group == 0, ftm[i, 1]]
			y0 <- Y[group == 0, ftm[i, 2]]
			z1 <- 0.5*log((1 + cor(x1, y1))/(1 - cor(x1, y1)))
			z0 <- 0.5*log((1 + cor(x0, y0))/(1 - cor(x0, y0)))
			u <- (z1 - z0)/sqrt(1/(n1 - 3) + 1/(n0 - 3))
			zsign[i] <- ifelse(abs(u) < 1.96, 0, sign(u))
			pv[i] <- 2*pnorm(-abs(u))
		}
	}
	if(is.null(group)) {
		for(i in 1:nrow(ftm)) {
			x <- Y[, ftm[i, 1]]
			y <- Y[, ftm[i, 2]]
			n <- nrow(Y)
			z <- 0.5*log((1 + cor(x, y))/(1 - cor(x, y)))
			u <- z/sqrt(1/(n - 3))
			zsign[i] <- ifelse(abs(u) < 1.96, 0, sign(u))
			pv[i] <- 2*pnorm(-abs(u))
		}
    }
	return (list(zsign, pv))
}

#' @title Transform data methods
#'
#' @description Implements various data trasformation methods with
#' optimal scaling for ordinal or nominal data, and to help relax
#' the assumption of normality (gaussianity) for continuous data.
#' 
#' @param x A matrix or data.frame (n x p). Rows correspond to subjects, and
#' columns to graph nodes.
#' @param method Trasform data method. It can be one of the following:
#' \enumerate{
#' \item "npn" (default), performs nonparanormal(npn) or semiparametric
#' Gaussian copula model (Liu et al, 2009), estimating the Gaussian copula
#' by marginally transforming the variables using smooth ECDF functions.
#' The npn distribution corresponds to the latent underlying multivariate
#' normal distribution, preserving the conditional independence structure
#' of the original variables.
#' \item "spearman", computes a trigonometric trasformation of Spearman
#' rho correlation for estimation of latent Gaussian correlations
#' parameter of a nonparanormal distribution (Harris & Dorton (2013),
#' and generates the data matrix with the exact same sample covariance
#' matrix as the estimated one.
#' \item "kendall", computes a trigonometric trasformation of Kendall
#' tau correlation for estimation of latent Gaussian correlations
#' parameter of a nonparanormal distribution (Harris & Dorton (2013),
#' and generates the data matrix with the exact same sample covariance
#' matrix as the estimated one.
#' \item "polychoric", computes the polychoric correlation matrix and
#' generates the data matrix with the exact same sample covariance matrix
#' as the estimated one. The polychoric correlation (Olsson, 1974) is a
#' measure of association between two ordinal variables. It is based on the
#' assumption that two latent bivariate normally distributed random variables
#' generate couples of ordinal scores. Tetrachoric (two binary variables) and
#' biserial (an ordinal and a numeric variables) correlations are special cases.
#' \item "lineals", performs optimal scaling in order to achieve linearizing
#' transformations for each bivariate regression between pairwise variables for
#' subsequent structural equation models using the resulting correlation
#' matrix computed on the transformed data (de Leeuw, 1988).
#' \item "mca", performs optimal scaling of categorical data by Multiple
#' Correspondence Analysis (MCA, a.k.a homogeneity analysis) maximizing
#' the first eigenvalues of the trasformed correlation matrix. The estimates
#' of the corresponding structural parameters are consistent if the underlying
#' latent space of the observed variables is unidimensional.
#' }
#' @param ... Currently ignored.
#'
#' @details Nonparanormal trasformation is computationally very efficient
#' and only requires one ECDF pass of the data matrix. Polychoric correlation
#' matrix is computed with the \code{lavCor()} function of the \code{lavaan}
#' package. Optimal scaling (lineals and mca) is performed with the
#' \code{lineals()} and \code{corAspect()} functions of the \code{aspect}
#' package (Mair and De Leeuw, 2008). To note, SEM fitting of the generate data
#' (fake data) must be done with a covariance-based method and bootstrap SE,
#' i.e., with \code{SEMrun(..., algo="ricf", n_rep=1000)}.
#' 
#' @return A list of 2 objects is returned:
#' \enumerate{
#' \item "data", the matrix (n x p) of n observations and p transformed
#' variables or the matrix (n x p) of simulate observations based on the
#' selected correlation matrix.  
#' \item "catscores", the category weights for "lineals" or "mca"
#' methods or NULL otherwise.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' 
#' Liu H, Lafferty J, and Wasserman L (2009). The Nonparanormal: Semiparametric Estimation of
#' High Dimensional Undirected Graphs. Journal of Machine Learning Research 10(80): 2295-2328
#' 
#' Harris N, and Drton M (2013). PC Algorithm for Nonparanormal Graphical Models.
#' Journal of Machine Learning Research 14 (69): 3365-3383
#' 
#' Olsson U (1979). Maximum likelihood estimation of the polychoric correlation coefficient.
#' Psychometrika, 44(4), 443-460.
#' 
#' Mair P, and De Leeuw J (2008). Scaling variables by optimizing correlational and
#' non-correlational aspects in R. Journal of Statistical Software, 32(9), 1-23.
#' 
#' de Leeuw J (1988). Multivariate analysis with linearizable regressions. Psychometrika,
#' 53, 437-454.
#'
#' @examples
#' 
#' #... with continuous ALS data
#' graph<- alsData$graph
#' data<- alsData$exprs; dim(data)
#' X<- data[, colnames(data) %in% V(graph)$name]; dim(X)
#'
#' npn.data<- transformData(X, method="npn")
#' sem0.npn<- SEMrun(graph, npn.data$data, algo="cggm")
#'
#' mvnS.data<- transformData(X, method="spearman")
#' sem0.mvnS<- SEMrun(graph, mvnS.data$data, algo="cggm")
#'
#' mvnK.data<- transformData(X, method="kendall")
#' sem0.mvnK<- SEMrun(graph, mvnK.data$data, algo="cggm")
#' 
#' #...with ordinal (K=4 categories) ALS data
#' Xord <- data.frame(X)
#' Xord <- as.data.frame(lapply(Xord, cut, 4, labels = FALSE))
#' colnames(Xord) <- sub("X", "", colnames(Xord))
#' 
#' mvnP.data<- transformData(Xord, method="polychoric")
#' sem0.mvnP<- SEMrun(graph, mvnP.data$data, algo="cggm")
#'
#' #...with nominal (K=4 categories) ALS data
#' mca.data<- transformData(Xord, method="mca")
#' sem0.mca<- SEMrun(graph, mca.data$data, algo="cggm")
#' mca.data$catscores
#' gplot(sem0.mca$graph, l="fdp", main="ALS mca")
#' 
#' # plot colored graphs
#' #par(mfrow=c(2,2), mar=rep(1,4))
#' #gplot(sem0.npn$graph, l="fdp", main="ALS npm")
#' #gplot(sem0.mvnS$graph, l="fdp", main="ALS mvnS")
#' #gplot(sem0.mvnK$graph, l="fdp", main="ALS mvnK")
#' #gplot(sem0.mvnP$graph, l="fdp", main="ALS mvnP")
#'
transformData <- function (x, method = "npn", ...)
{
	n <- nrow(x)
	p <- ncol(x)
	x.col <- colnames(x)
	x.row <- rownames(x)
	catscores <- NULL
	
	if (method == "npn") {
		cat("Conducting the nonparanormal transformation via shrunkun ECDF...")
		z <- qnorm(apply(x, 2, rank)/(n + 1))
		z <- z/sd(z[, 1])	
	}
	if (method == "spearman") { 
		cat("Simulating gaussian data via Spearman correlations...")
		x <- 2 * sin(pi/6 * cor(x, method = "spearman"))
		z <- generateData(Sest = x, n = n, p = p)
	}
	if (method == "kendall") { 
		cat("Simulating gaussian data via Kendall correlations...")
		x <- sin(pi/2 * cor(x, method = "kendall"))
		z <- generateData(Sest = x, n = n, p = p)
    }
	if (method == "polychoric") { 
		cat("Simulating gaussian data via polychoric correlations...")
		#x <- suppressWarnings(lavCor(x, ordered = names(x))[1:p,1:p])
		x <- lavCor(x, ordered=names(x), cor.smooth=TRUE)[1:p,1:p]
		z <- generateData(Sest = x, n = n, p = p)
	}
	if (method == "lineals") {
		cat("Conducting the optimal (ordinal) linearizing transformation...")
		z <- aspect::lineals(x, level = "ordinal")
		catscores<- z$catscores
		z <- sqrt(n-1)*z$scoremat
    }
	if (method == "mca") {
		cat("Conducting the first solution of Multiple Correspondence Analysis...")
		z <- aspect::corAspect(x, aspect = "aspectEigen")
		catscores<- z$catscores
		z <- sqrt((n-1)/n)*z$scoremat
	}

	cat("done.\n")
	colnames(z) <- x.col
	rownames(z) <- x.row
	v0 <- which(apply(z, 2, var) == 0)
	if (length(v0) > 0) z <- cbind(z[,-v0], x[,v0])

	return(list(data = z, catscores = catscores))
}

generateData <- function(Sest, n, p, ...)
{
	if (!corpcor::is.positive.definite(Sest)){
	 Sest <- corpcor::cor.shrink(Sest, verbose = FALSE)[1:p,1:p]
	}
	e <- eigen(Sest)
	sqrt.true.cov.mat <- e$vectors%*%sqrt(diag(e$values))

	fake.data <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sest)
	samp.cov.mat <- cov(fake.data)
	e <- eigen(samp.cov.mat)
	sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))

	fake.data <- t(sqrt.true.cov.mat%*%solve(sqrt.samp.cov.mat,t(fake.data)))
	fake.data <- as.data.frame(fake.data)

	return(fake.data)
}

#' @title Import pathways and generate a reference network
#'
#' @description Utility to create pathway lists as igraph objects and 
#' interaction networks from Reactome, KEGG, and other pathway databases.
#'
#' @param db String indicating the database name. Please, check the 
#' \code{\link[graphite]{pathways}} function from \pkg{graphite} to list
#' the available datasets.
#' @param organism A string indicating the source organism. Please, check 
#' the \code{\link[graphite]{pathways}} function from \pkg{graphite} to list
#' the available datasets (default = "hsapiens") 
#' @param id_type Gene ID type. The default is set to "ENTREZID" (standard 
#' SEM fitting nomenclature). A common choice could be "SYMBOL", for 
#' official gene symbols.
#' @param lcc A logical value. If TRUE (default), the reference network 
#' will only include the largest connected component. It will include all 
#' disconnected components otherwise.
#' @param ... Currently ignored.
#'
#' @details This function uses \code{graphite} to download and preprocess 
#' network data from pathway databases. The output is then created using 
#' igraph and SEMgraph utilities.
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item a list of pathways ad igraph objects;
#' \item the union of graphs in the pathway list.
#' }
#'
#' @export
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @references
#'
#' Sales G, Calura E, Cavalieri D, Romualdi C (2012). graphite - a Bioconductor 
#' package to convert pathway topology to gene network. 
#' BMC Bioinformatics. 
#' <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-20>
#'
#' @examples
#'
#' \dontrun{
#' 
#' # Create KEGG reference pathway list and reference network for Homo sapiens
#' kegg.hs <- loadPathways("kegg", "hsapiens", "ENTREZID")
#' 
#' # Inspect results
#' names(kegg.hs$pathways)
#' kegg.hs$network
#'
#' }
#'
loadPathways <- function(db, organism = "hsapiens", id_type = "ENTREZID", lcc = TRUE, ...)
{
	# Import pathway from database
	pathways <- graphite::pathways(organism, db)
	pathway.list <- list()
	n <- length(pathways)
	del <- NULL

	for (i in 1:n) { #i=51
	  pw.name <- names(pathways)[[i]]
	  #if (i %% 10 == 0) message("  Processed ", i, "/", n, " pathways")
	  #message(paste0("Pathway ", i, " of ", n, ": ", pw.name))

	  # Node ID conversion
	  pw <- graphite::convertIdentifiers(pathways[[i]], id_type)
	  
	  # Conversion pathway -> GraphNEL -> igraph
	  G <- igraph::graph_from_graphnel(graphite::pathwayGraph(pw))
	  G <- igraph::simplify(G, remove.loops = TRUE)
	  G <- G - E(G)[igraph::which_mutual(G)]

	  # Skip small graph
	  if (igraph::vcount(G) <= 5 || igraph::ecount(G) == 0) {
	  message(paste0("Delete pathway ", i, " of ", n, ": ", pw.name))
	  del <- c(del, i); next
	  }

	  #Removing node name prefix and edge weight
	  V(G)$name <- sub(paste0(id_type, ":"), "", V(G)$name)
	  G <- igraph::delete_edge_attr(G, "weight")
	  
	  # Adding graph to the list
	  pathway.list[[i]] <- G
	  names(pathway.list)[[i]] <- pw.name
	}

	# Generating reference network
	reference <-  Reduce(igraph::union, pathway.list)
	if (lcc) reference <- properties(reference)[[1]] #SEMgraph

	# Add node labels to all graphs
	db <- org.Hs.eg.db::org.Hs.eg.db
	add_labels <- function(g) {
	 V(g)$label <- suppressMessages(AnnotationDbi::mapIds(db, V(g)$name, "SYMBOL", id_type))
	 g
	}
	pathway.list <- lapply(pathway.list[-del], add_labels)
	reference <- add_labels(reference)
	
	return(list(pathways = pathway.list, network = reference))
}
