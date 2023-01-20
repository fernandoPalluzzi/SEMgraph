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

#' @title Bow-free covariance search and data de-correlation
#'
#' @description \code{SEMbap()} function implements different deconfounding
#' methods to adjust the data matrix by removing latent sources of confounding
#' encoded in them. The selected methods are either based on: (i) Bow-free
#' Acyclic Paths (BAP) search, (ii) LVs proxies as additional source nodes of
#' the data matrix, X or (iii) spectral transformation of X. 
#' 
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If \code{NULL} (default), confouding within group
#' will not be considered.
#' @param dalgo Deconfounding method. Five algorithms are available: 
#' \itemize{
#' \item "cggm" (default). The "cggm" algorithm make un exhaustive search of
#' missing edges with significant covariance (see details). The inverse of
#' the selected covariance matrix (i.e. the precision matrix, W) is fitted by
#' a constrained gaussian graphical model (cggm), and a de-correlated data
#' matrix, Z is obtained multiplying the data matrix, X rightward by the
#' square root of the estimated precision matrix, Z=XW^(1/2) as suggested by
#' Grassi, Palluzzi and Tarantino (2022).
#' \item "glpc". Similarly to "cggm", "glpc" algorithm first makes an 
#' exhaustive search of missing with significant covariance (see details).
#' Once obtained the adjacency matrix of the covariances, Graph-Laplacian PCA
#' (gLPCA) algorithm learns a low dimensional representation of the observed
#' data matrix that incorporates graph structures (Jiang et al., 2013).
#' Then, the DAG is extended by including the confounding  proxies, i.e. LVs,
#' as additional source nodes defined by last q principal component scores
#' of gLPCA and these LV scores are added to the data matrix, Z=cbind(LV,X). 
#' \item "ml". The procedure of "ml" algorithm is analogous to add additional 
#' source nodes to DAG as in "glpc" algorithm, but confounding proxies are
#' the q factor scores extracted via Factor Analysis (FA) with quasi-maximum
#' likelihood estimates, iteratively solved by Expectation-Maximization (EM)
#' algorithm using the PCA solution as the initial value (Bai and Li, 2012).
#' \item "trim". Ćevid et al. (2020) suggest multiplying the data matrix, X
#' leftward by a well selected spectrum transformation matrix, F which modifies
#' the singular values of X, while keeping its singular vectors intact, Z=FX.
#' Trim transform limits all singular values to be at most some costant (t),
#' usually their median is a good choice of t.
#' \item "pcss". This procedure is analogous to applying a spectral transformation,
#' Z=FX but the new matrix of singular values is obtained by mapping the first
#' q singular values to 0 (Ćevid et al., 2020).
#' }
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, \code{method}
#' is set to "BH" (i.e., Benjamini-Hochberg multiple test correction).
#' @param alpha Significance level for false discovery rate (FDR) used
#' for Shipley's local d-separation tests. This argument is used to
#' control data de-correlation. A higher \code{alpha} level includes more
#' hidden covariances, thus considering more sources of confounding.
#' If \code{alpha = 0}, data de-correlation is disabled.
#' By default, \code{alpha = 0.05}.
#' @param limit An integer value corresponding to the number of missing
#' edges of the extracted acyclic graph. Beyond this limit, multicore
#' computation is enabled to reduce the computational burden.
#' By default, \code{limit = 30000}.
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen.
#' @param ... Currently ignored.
#'
#' @details Missing edges in causal network inference using a directed acyclic
#' graph (DAG) are frequently hidden by unmeasured confounding variables.
#' A Bow-free Acyclic Paths (BAP) search is performed with d-separation tests
#' between all pairs of variables with missing connection in the input DAG,
#' adding a bidirected edge (i.e., bow-free covariance) to the DAG when there
#' is an association between them. The d-separation test evaluates if two
#' variables (Y1, Y2) in a DAG are conditionally independent for a given
#' conditioning set, C represented in a DAG by the union of the parent sets
#' of Y1 and Y2 (Shipley, 2000). A new bow-free covariance is added if there
#' is a significant (Y1, Y2) association at a significance level \code{alpha},
#' after multiple testing correction. The selected covariance between pairs of
#' nodes is interpreted as the effect of a latent variable (LV) acting on both
#' nodes; i.e., the LV is an unobserved confounder. BAP-based algorithms
#' adjust (or de-correlate) the observed data matrix by conditioning out the
#' latent triggers responsible for the nuisance edges. For LV algorithms the
#' number of hidden proxies, q is determined through scree plot rules: look
#' the eigenvalues, i.e., singular values^2, at the knee point (if dalgo =
#' "glpc" or "ml") or look at the first eingenvalue cluster (if dalgo = "pcss").
#' If the input graph is not acyclic, a warning message will be raised, and a
#' cycle-breaking algorithm will be applied (see \code{\link[SEMgraph]{graph2dag}}
#' for details).
#'
#' @return A list of four objects:
#' \itemize{
#' \item "dag", the directed acyclic graph (DAG) extracted from input graph.
#' If (dalgo = "glpc" or "ml"), the DAG also includes LVs as source nodes.
#' \item "dsep", the data.frame of all d-separation tests over missing edges
#' in the DAG. If (dalgo != "cggm" or "glpc"), dsep dataframe is equal to NULL.
#' \item "adj", the adjacency matrix of selected covariances; i.e, the 
#' missing edges selected after multiple testing correction. If (dalgo != "cggm"
#' or "glpc"), adj matrix is equal to NULL.
#' \item "data", the adjusted (de-correlated) data matrix or if (dalgo = "glpc",
#' or "ml"), the combined data matrix, where the first columns represent LVs
#' scores and the other columns are the raw data.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs.
#' Structural Equation Modeling, 7(2), 206-218.
#' <https://doi.org/10.1207/S15328007SEM0702_4>
#'
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38(20), 4829–4830.
#' <https://doi.org/10.1093/bioinformatics/btac567>
#' 
#' Jiang B, Ding C, Bin L, Tang J (2013). Graph-Laplacian PCA: 
#' Closed-Form Solution and Robustness. IEEE Conference on Computer
#' Vision and Pattern Recognition, 3492-3498. 
#' <https://doi.org/10.1109/CVPR.2013.448>
#' 
#' Jushan Bai and Kunpeng Li (2012). Statistical Analysis of Factor Models of High
#' Dimension. The Annals of Statistics, 40 (1), 436-465
#' <https://doi.org/10.1214/11-AOS966> 
#' 
#' Ćevid D,  Bühlmann P, Meinshausen N (2020). Spectral deconfounding via
#' perturbed sparse linear models. J. Mach. Learn. Res, 21 (232), 1-41.
#' <http://jmlr.org/papers/v21/19-545.html>
#'
#' @examples
#'
#' # Model fitting
#' sem0 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc))
#' sem1 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc), group = sachs$group)
#'
#' # BAP search with default method (dalgo="bap")
#' bap <- SEMbap(graph = sachs$graph, data = log(sachs$pkc), verbose = TRUE)
#' #gplot(bap$dag)
#'
#' # Model fitting with de-correlated data
#' sem0 <- SEMrun(graph = bap$dag, data = bap$data)
#'
#' # BAP search within group with gLPCA scores
#' glpc <- SEMbap(graph = sachs$graph, data = log(sachs$pkc), group = sachs$group, 
#'                dalgo= "glpc", verbose = FALSE)
#' gplot(glpc$dag)
#'
#' # Model fitting (node perturbation) with gLPCA scores
#' sem1 <- SEMrun(graph = glpc$dag, data = glpc$data, group = sachs$group)
#'
SEMbap <- function(graph, data, group=NULL, dalgo="cggm", method="BH",
				   alpha=0.05, limit=30000, verbose=FALSE, ...)
{
	# Set graph and data objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	dataY<- as.matrix(data[,nodes])
	dag<- graph2dag(graph, dataY)

	if (dalgo == "cggm" | dalgo == "glpc") {
	  # d-separation local tests (B_U):
	  dsep<- Shipley.test(dag, dataY, limit=limit, verbose=FALSE)$dsep
	  d_sep<- subset(dsep, p.adjust(dsep$p.value, method = method) < alpha)
	  guu<- graph_from_data_frame(d_sep[,1:2], directed=FALSE)
	  if (ecount(guu) > 0) {
	   guu<- difference(guu, as.undirected(dag))
	   cat("Number of significant local tests:", nrow(d_sep), "/", nrow(dsep),"\n\n")
	  }else{
	   return(message("NULL covariance graph: ALL adjusted pvalues > ", alpha, "!"))
	  }

	  # Set the complete adjacency matrix of guu
	  Adj <- as_adj(as.undirected(dag), type = "both", sparse = FALSE)
	  adj <- as_adj(guu, type = "both", sparse = FALSE)
	  idx <- which(rownames(Adj) %in% rownames(adj) == FALSE)
	  if (length(idx) > 0) {
	    R <- matrix(0, length(idx), ncol(adj))
	    C <- matrix(0, nrow(adj), length(idx))
	    D <- matrix(0, length(idx), length(idx))
	    adj <- rbind(cbind(D,R), cbind(C,adj))
	    rownames(adj)[1:length(idx)] <- rownames(Adj)[idx]
	    colnames(adj)[1:length(idx)] <- rownames(Adj)[idx]
	  } #colnames(adj)

	  if (verbose) {
	   # Covariance and latent variables graphs (guu, gLV)
	   ftm<- as_edgelist(as.undirected(guu))
	   ftmLV<- NULL
	   V(guu)$color<- "white"
	   for (i in 1:nrow(ftm)) ftmLV<- rbind(ftmLV, cbind(rep(paste0("L",i), 2),ftm[i,]))
	   gLV<- graph_from_data_frame(ftmLV, directed=TRUE)
	   V(gLV)$color<- ifelse(substr(V(gLV)$name,1,1)=="L","yellow","white")

	   plot(guu, main="extended covariance graph (guu)")
	   Sys.sleep(3)
	   plot(gLV, main="extended latent variables graph (gLV)")
	   Sys.sleep(3)
	  }

	  if (dalgo == "cggm")
 	   Z <- estimatePsi(adj, dataY, group=group, dalgo="cggm")
	  if (dalgo == "glpc") {
	   Z<- estimateLV(adj, dataY, group=group, dalgo="glpc")
	   dag <- hiddenGraph(dag, Z)
	  }
	}

	else if (dalgo == "pc" | dalgo == "ml"){
	 Z <- estimateLV(adj=as_adj(dag), dataY, group=group, dalgo=dalgo)
	 dag <- hiddenGraph(dag, Z)
	 dsep <- adj <- NULL
	}else{
	 Z <- estimateFX(adj=as_adj(dag), dataY, group=group, dalgo=dalgo)
	 dsep <- adj <- NULL
	}

	# SEM fitting with adjusted bow-free covariances:
	if (verbose) fit<- SEMrun(dag, Z, algo = "ricf", n_rep = 0)

	return( list(dag=dag, dsep=dsep, adj=adj, data=Z) )
}

estimatePsi <- function(adj, data, group, dalgo, ...)
{
	# Set data objects
	if (is.null(group)){
	 Y <- data[, which(colnames(data) %in% colnames(adj))]
	}else{
	 Y_0 <- data[group == 0, which(colnames(data) %in% colnames(adj))]
	 Y_1 <- data[group == 1, which(colnames(data) %in% colnames(adj))]
	}

	# Estimate deconfounding data via covariance/precision matrix

	estimate_Psi <- function(X)
	{
	 X <- scale(X)
	 p <- ncol(X)
	 
	 if (dalgo == "lavaan") {
	   diag(adj)[rowSums(adj) == 0] <- 1
	   ug<- graph_from_adjacency_matrix(adj, mode="undirected")#plot(ug)
	   sem0<- quiet(SEMrun(graph=ug, data=X, algo="lavaan", SE="none"))
	   wi<- inspect(sem0$fit, "est")$theta[1:p,1:p]
	   colnames(wi)<- rownames(wi)<- gsub("z", "", colnames(wi))
	   wi<- wi[rownames(adj),colnames(adj)]
	 }
	 if (dalgo == "cggm") { #HTF
	   S <- cov(X[, colnames(adj)])
	   wi <- ggm::fitConGraph(adj, S, n=nrow(X))$Shat
	   colnames(wi) <- rownames(wi) <- colnames(adj)
	 }
	 #if (verbose) image(solve(wi) != 0, main=dalgo)

	 if(!corpcor::is.positive.definite(wi)){
	   wi<- corpcor::cov.shrink(wi, verbose = FALSE)
	   #wi<- corpcor::cor.shrink(wi, verbose = TRUE)
	   #wi<- corpcor::make.positive.definite(wi)
	 }
	 E<- eigen(wi) # Eigenvalues-eigenvectors of w
	 #R<- E$vectors%*%diag(sqrt(E$values))%*%t(E$vectors) #dim(R)
	 #sum(wi - R %*% R)
	 R<- E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors) #dim(R)
	 #sum(solve(wi) - R %*% R)
	 colnames(R) <- rownames(R) <- colnames(adj)
	 X <- X[, colnames(adj)]
	 return( as.matrix(X)%*%R )
	}
	
	if (is.null(group)){
     Z <- estimate_Psi(Y)
	}else{
	 Z <- rbind(
	  estimate_Psi(Y_0),
	  estimate_Psi(Y_1))
	}

	return(data = Z)
}

estimateFX <- function(adj, data, group, dalgo, ...)
{
	# Set data objects
	K <- 0
	if (is.null(group)){
	 Y <- data[, which(colnames(data) %in% colnames(adj))]
	 if (dalgo == "pcss") K<- screePlot(Y, method="cluster")
	}else{
	 Y_0 <- data[group == 0, which(colnames(data) %in% colnames(adj))]
	 Y_1 <- data[group == 1, which(colnames(data) %in% colnames(adj))]
	 if (dalgo == "pcss") K <- screePlot(rbind(Y_0, Y_1), method="cluster")
	}

	# Estimate deconfounding data via svd(X, nu, nv)

	estimate_FX <- function(X, q)
	{
	 X <- scale(X)
	 n <- nrow(X)  

	 if (dalgo == "trim") {
	  x <- svd(X, nv = 0)
	  U <- sqrt(n)*x$u #round(t(U)%*%U,6)
	  M <- median(x$d)
	  d <- vector()
	  for (i in 1:length(x$d)) d[i]<- min(x$d[i],M)/x$d[i]
	  return( U%*%diag(d/n)%*%t(U)%*%X )
	 }
	
	 if (dalgo == "lava") {
	  x <- svd(X, nv = 0)
	  U <- sqrt(n)*x$u #round(t(U)%*%U,6)
	  M <- median(x$d)
	  d <- sqrt(n*M*x$d^2/(n*M+x$d^2))/x$d
	  return( U%*%diag(d/n)%*%t(U)%*%X )
	 }

	 if (dalgo == "rsvp") {
	  x <- svd(X, nv = 0)
	  U <- sqrt(n)*x$u #round(t(U)%*%U,6)
	  d <- 1/x$d
	  return( U%*%diag(d/n)%*%t(U)%*%X )
	 }
    
     if (dalgo == "pcss") {
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  U <- sqrt(n)*x$u #round(t(U)%*%U,6)
	  d <- c(rep(0,q),x$d[-c(1:q)])/x$d
	  return( U%*%diag(d/n)%*%t(U)%*%X )
	 }
	}
	
	if (is.null(group)){
	 Z <- estimate_FX(Y, q=K)
	}else{
	Z <- rbind(
	  estimate_FX(Y_0, q=K),
	  estimate_FX(Y_1, q=K))
	}

	return(data = Z)
}

estimateLV <- function(adj, data, group, dalgo, ...)
{
	# Set data objects
	if (is.null(group)){
	 Y <- data[, which(colnames(data) %in% colnames(adj))]
	 K <- screePlot(Y, method="knee")
	}else{
	 Y_0 <- data[group == 0, which(colnames(data) %in% colnames(adj))]
	 Y_1 <- data[group == 1, which(colnames(data) %in% colnames(adj))]
	 K <- screePlot(rbind(Y_0, Y_1), method="knee")
	}

	# Estimate latent variables via svd(X, nu, nv)

	estimate_LV <- function(X, q)
	{
	 X <- scale(X)
	 n <- nrow(X)
	 p <- ncol(X)

	 if (dalgo == "pc") {
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  LV <- sqrt(n-1)*as.matrix(x$u[,1:q]) #round(t(LV)%*%LV,6)
	  colnames(LV) <- paste0("LV", seq_len(q))
	  return(cbind(LV,X))
	 }

	 if (dalgo == "ml") {
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  LV <- factor.analysis(X, r=q, method="ml")$Z
	  colnames(LV) <- paste0("LV", seq_len(q))
	  return(cbind(LV,X))
	 }

	 if (dalgo == "glpc") {
	  S <- cov(X)*(n-1)/n
	  E <- eigen(S)
	  e1 <- E$values[1] # e1
	  
	  dj <- rowSums(adj, na.rm = FALSE)
	  D <- diag(dj)
	  L <- D - adj
	  E <- eigen(L)
	  f1 <- E$values[1] # f1

	  one <-rep(1,p)
	  I <- diag(one)
	  if (sum(D) == 0) {
	   b <- 0; f1 <- 1
	  }else{b <- 0.5}

	  G <- (1-b)*(I-S/e1) + b*(L/f1 + one%*%t(one)/p)
	  E <- eigen(G)
	  e <- E$values
	  W <- E$vectors #sum(round(t(W)%*%W, 6))
	  U <- X%*%W

	  #ev<- sort(1-e, decreasing = TRUE)
	  #q <- screePlot(ev, method="knee")
	  if (q == 0) return(X)
	  LV <- scale(U[,(p-q+1):p])
	  colnames(LV) <- paste0("LV", seq_len(q))
	  return(cbind(LV,X))
     }
	}
	
	if (is.null(group)){
	 Z <- estimate_LV(Y, q=K)
	}else{
	 Z <- rbind(
	  estimate_LV(Y_0, q=K),
	  estimate_LV(Y_1, q=K))
	}

	return(data = Z)
}

screePlot <- function(X, method, ...)
{
	ev <- svd(X, nu = 0, nv = 0)$d^2
	
	if (method == "knee") {
	 # look at the knee point of a scree plot
	 scree <- ev
	 scree <- scree[seq_len(round(length(scree) / 2))]
	 values <- seq(length(scree))

	 d1 <- diff(scree) / diff(values) # first derivative
	 d2 <- diff(d1) / diff(values[-1]) # second derivative
	 idx <- which.max(abs(d2))
	}

	if (method == "cluster") {
	 # look at the first cluster of k-means clustering
	 scree <- ev
	 clust <- kmeans(
		scree,
		centers = c(
		 scree[1],
		 scree[2],
		 scree[round(length(scree) / 2 + 1)],
		 scree[length(scree)]
		)
	 )$cluster
	 idx <- sum(clust == clust[1])
	}

	# visualize the scree plot
	r <- length(ev)
	plot(c(1:r), ev[1:r])
	abline(h=ev[idx]-0.006, lty=2)
	pve <- round(cumsum(ev[1:idx])[idx]/sum(ev),2)
	#logger::log_info("Estimated {idx} latent confounders")
	message(paste0("Estimated ",idx," (", pve,") latent confounders"))

	return(idx)
}

hiddenGraph <- function(graph, data, cg=NULL, EntrezID=FALSE, verbose=FALSE, ...)
{
	VH <- colnames(data)[grepl("LV",colnames(data))]
	gH <- graph + igraph::vertices(VH)
	if (EntrezID) {
	 V(gH)$label <- mapIds(org.Hs.eg.db, V(gH)$name, 'SYMBOL', 'ENTREZID')
	 V(gH)$label[is.na(V(gH)$label)] <- VH
	}
	E <- NULL
	 for(v in 1:length(VH)){
	  if (!is.null(cg)) graph<- cg[[v]]
	  for(i in 1:vcount(graph)){ #i=1
		E <- c(E, VH[v], V(graph)$name[i])
	  }
	 }
	gH <- gH + igraph::edges(E)
	if (verbose) gplot(gH, l="fdp")
	
	return(gH)	
}

#' @title Missing edge testing implied by a graph with Shipley's basis-set
#'
#' @description Compute all the P-values of the d-separation tests
#' implied by the missing edges of a given acyclic graph (DAG).
#' The conditioning set Z is represented, in a DAG, by the union of the
#' parent sets of X and Y (Shipley, 2000). 
#' The results of every test, in a DAG, is then combined using the
#' Fisher’s statistic in an overall test of the fitted model
#' C = -2*sum(log(P-value(k))), where C is distributed as a chi-squared
#' variate with df = 2k, as suggested by Shipley (2000).
#'
#' @param graph A directed graph as an igraph object.
#' @param data A data matrix with subjects as rows and variables as
#' columns.
#' @param MCX2 If TRUE, a Monte Carlo P-value of the combined C test is enabled.
#' @param verbose If TRUE, Shipley's test results will be showed to
#' screen (default = TRUE).
#' @param limit An integer value corresponding to the number of missing
#' edges of the extracted acyclic graph. Beyond this limit, multicore
#' computation is enabled to reduce the computational burden.
#' By default, \code{limit = 30000}.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of three objects: (i) "dag":  the DAG used to perform the Shipley
#' test (ii) "dsep": the data.frame of all d-separation tests over missing edges in
#' the DAG and (iii) "ctest": the overall Shipley's' P-value.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs.
#' Structural Equation Modeling, 7(2): 206-218.
#' <https://doi.org/10.1207/S15328007SEM0702_4>
#'
#' @examples
#'
#' #\donttest{
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' sem <- SEMrun(alsData$graph, als.npn)
#' C_test <- Shipley.test(sem$graph, als.npn, MCX2 = FALSE)
#' #MC_test <- Shipley.test(sem$graph, als.npn, MCX2 = TRUE)
#'
#' #}
#'
Shipley.test<- function(graph, data, MCX2=FALSE, limit=30000, verbose=TRUE,...)
{
	# graph to DAG conversion :
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	df1<- vcount(graph)*(vcount(graph)-1)/2-ecount(as.undirected(graph))
	dataY<- as.matrix(data[, nodes])
	
	if (!is_dag(graph)){
	 cat("WARNING: input graph is not acyclic !\n")
	 cat(" Applying graph -> DAG conversion...\n")
	 dag<- graph2dag(graph, dataY, bap=FALSE) #del cycles & all <->
	 df2<- vcount(dag)*(vcount(dag)-1)/2-ecount(as.undirected(dag))
	 cat(" \nDegrees of freedom:\n Input graph  =", 
            df1, "\n Output graph =", df2, "\n\n")
	}else{
	 dag <- graph
	 df2 <- df1
	}

	# d-separation local tests (B_U) & Shipley's overall pvalue
	dsep<- dsep.test(dag=dag, S=cov(dataY), n=nrow(dataY), limit=limit)
	#Combining p-values with Fisher's procedure:
	X2<- -2 * sum(log(dsep$p.value + 1E-16))
	df<- 2 * nrow(dsep)
	if (MCX2) {
	 pv<- MCX2(model.df=df, n.obs=nrow(data), model.chi.square=X2)[[1]]
	}else{
	 pv<- 1 - pchisq(q = X2, df = df)
	}
	if (verbose) print(data.frame(C_test=X2, df=df, pvalue=round(pv,6)))
			
	return( list(dag=dag, dsep=dsep, ctest=c(X2, df, pv)) )
}

dsep.test <- function(dag, S, n, limit, ...)
{
 	# d-sep (basis set) testing of a DAG
	A <- ifelse(as_adj(as.undirected(dag), sparse=FALSE) == 1, 0, 1)
	ug <- graph_from_adjacency_matrix(A, mode="undirected", diag=FALSE)
	M <- attr(E(ug), "vnames")
	
	local<- function(x) {
	 s <- strsplit(x,"\\|")
	 ed <- c(s[[1]][1], s[[1]][2])
	 pa.r <- V(dag)$name[SEMgraph::parents(dag, ed[1])]
	 pa.s <- V(dag)$name[SEMgraph::parents(dag, ed[2])]
	 dsep <- union(pa.r, pa.s)
     dsep <- setdiff(dsep, ed)
	 B <- c(ed, dsep)
	 if(length(B) > (n-3)) return(rep(NA,4))
	 p.value <- pcor.test(S, B, n, H0=0.05)
	 set <- paste(B[-c(1:2)], collapse=",")
	 return(data.frame(X=B[1], Y=B[2], SET=set, p.value))
	}

	message("d-separation test (basis set) of ", length(M), " edges...")
	op<- pbapply::pboptions(type = "timer", style = 2)
	df<- vcount(dag)*(vcount(dag)-1)/2 - ecount(dag)
	if ( df > limit ){
	 n_cores <- parallel::detectCores(logical = FALSE)
	 cl <- parallel::makeCluster(n_cores)
	 parallel::clusterExport(cl, c("local", "dag", "S", "n"),
							 envir = environment())
	 SET<- pbapply::pblapply(M, local, cl=cl)
	 parallel::stopCluster(cl)
	}else{
	 SET<- pbapply::pblapply(M, local, cl=NULL)
	}
	SET<- do.call(rbind, lapply(SET, as.data.frame))

	return(SET = na.omit(SET))
}

pcor.test<- function(S, B, n, H0, ...)
{
	#Set objects
	k <- solve(S[B,B])
    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
	q <- length(B)-2
	if( H0 == 0 ) {
	#Test null H0: r=abs(r(X,Y|Z))=0
	 df <- n - 2 - q
	 tval <- r * sqrt(df)/sqrt(1 - r * r)
	 pval <- 2 * pt(-abs(tval), df)
	}else{
	#Test of not-close fit, H0: r=abs(r(X,Y|Z)) vs. r<.05
	 z <- atanh(r)
	 se <-  1/sqrt(n - 3 - q)
	 pval <- pchisq((z/se)^2, df=1, ncp=(atanh(H0)/se)^2, lower.tail=FALSE)
	}
	return(pval)
}

MCX2 <- function (model.df, n.obs, model.chi.square, n.sim = 10000, ...)
{
	#Monte Carlo Chi-square simulator (Author: Bill Shipley) from:
	#devtools::install_github("BillShipley/CauseAndCorrelation")
	# All rights reserved.  See the file COPYING for license terms.
	x <- (-1 + sqrt(1 + 8 * model.df))/2
	if ((x - as.integer(x)) == 0)
	v <- x
	if ((x - as.integer(x)) > 0 & (x - as.integer(x)) < 1) 
	v <- as.integer(x) + 1
	if ((x - as.integer(x)) > 1)return(message("ERROR: check model df !"))
	c.value <- v * (v + 1)/2 - model.df
	MCX2 <- rep(NA, n.sim)
	for (i in 1:n.sim) {
	 dat <- matrix(rnorm(n.obs * v), ncol = v)
	 obs.VCV <- var(dat)
	 model.VCV <- diag(v)
	 diag(model.VCV)[1:c.value] <- diag(obs.VCV)[1:c.value]
	 MCX2[i] <- (n.obs - 1) * (log(det(model.VCV)) + sum(diag(obs.VCV) * 
							(1/diag(model.VCV))) - log(det(obs.VCV)) - v)
	}
	MCprob <- sum(MCX2 >= model.chi.square)/n.sim
	x <- seq(0, max(MCX2))
	theoretical.prob <- dchisq(x, model.df)
	MLprob<- pchisq(model.chi.square, model.df, lower.tail=FALSE)
	
	return(list(MCprob = MCprob, MLprob = MLprob))
}

#' @title Estimate the optimal DAG from an input graph
#'
#' @description Extract the optimal DAG from an input graph, using
#' topological order or top-down order search and LASSO-based algorithm, implemented in
#' \code{\link[glmnet]{glmnet}}.
#'
#' @param graph An igraph object.
#' @param data A matrix whith n rows corresponding to subjects, and p columns
#' to graph nodes (variables).
#' @param LO character for linear order method. If LO="TO" the topological
#' order of the input DAG is enabled (default), while LO="TD" the data-driven
#' top-down minimum conditional variance method is performed.
#' @param beta Numeric value. Minimum absolute LASSO beta coefficient for
#' a new interaction to be retained in the final model. By default,
#' \code{beta = 0}.
#' @param lambdas A vector of regularization LASSO lambda values. If lambdas is
#' NULL, the \code{\link[glmnet]{glmnet}} default using cross-validation lambdas
#' is enabled. If lambdas is NA (default), the tuning-free scheme is enabled by
#' fixing lambdas = sqrt(log(p)/n), as suggested by Janková and van de Geer
#' (2015) and many others. This will both reduce computational time and provide
#' the same result at each run.
#' @param penalty A logical value. Separate penalty factors can be applied to
#' each coefficient. This is a number that multiplies lambda to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. If TRUE (default) weights are based
#' on the graph edges: 0 (i.e., edge present) and 1 (i.e., missing edge) ensures
#' that the input edges will be retained in the final model. If FALSE the
#' \code{\link[glmnet]{glmnet}} default is enabled (all weights equal to 1). Note:
#' the penalty factors are internally rescaled to sum p (the number of variables).
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen.
#' @param ... Currently ignored.
#'
#' @details The optimal DAG is estimated using the order search approach. First
#' a linear order of p nodes is determined, and from this sort, the DAG can be
#' learned using successive penalized (L1) regressions (Shojaie and Michailidis,
#' 2010). The estimate linear order are obtained from \emph{a priori} graph
#' topological order (TO), or with a data-driven high dimensional top-down (TD)
#' approach (best subset regression), assuming a SEM whose error terms have equal
#' variances (Peters and Bühlmann, 2014; Chen et al, 2019). If the input graph is
#' not acyclic, a warning message will be raised, and a cycle-breaking algorithm
#' will be applied (see \code{\link[SEMgraph]{graph2dag}} for details).
#' Output DAG edges will be colored in gray, if they were present in the
#' input graph, and in green, if they are new edges generated by LASSO
#' screening.
#'
#' @return A list of 3 igraph objects:
#' \enumerate{
#' \item "dag", the estimated DAG;
#' \item "dag.new", new estimated connections;
#' \item "dag.old", connections preserved from the input graph.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @seealso \code{\link[SEMgraph]{modelSearch}}
#'
#' @references
#'
#' Tibshirani R, Bien J, Friedman J, Hastie T, Simon N, Taylor J,
#' Tibshirani RJ (2012). Strong rules for discarding predictors in
#' lasso type problems. Royal Statistical Society: Series B
#' (Statistical Methodology), 74(2): 245-266.
#' <https://doi.org/10.1111/j.1467-9868.2011.01004.x>
#' 
#' Shojaie A, Michailidis G (2010). Penalized likelihood methods for
#' estimation of sparse high-dimensional directed acyclic graphs.
#' Biometrika, 97(3): 519-538. <https://doi.org/10.1093/biomet/asq038>
#'
#' Jankova J, van de Geer S (2015). Confidence intervals for high-dimensional
#' inverse covariance estimation. Electronic Journal of Statistics,
#' 9(1): 1205-1229. <https://doi.org/10.1214/15-EJS1031>
#'
#' Peters J, Bühlmann P (2014). Identifiability of Gaussian structural equation
#' models with equal error variances. Biometrika, 101(1):219–228.
#' 
#' Chen W, Drton M, Wang YS (2019). On Causal Discovery with an Equal-Variance
#' Assumption. Biometrika, 106(4): 973-980.
#'
#' @examples
#'
#' # DAG estimation
#' G <- SEMdag(graph = sachs$graph, data = log(sachs$pkc), beta = 0.05)
#'
#' # Model fitting
#' sem <- SEMrun(graph = G$dag, data = log(sachs$pkc), group = sachs$group)
#'
#' # Graphs
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(1,4))
#' plot(sachs$graph, layout=layout.circle, main="input graph")
#' plot(G$dag, layout=layout.circle, main = "Output DAG")
#' plot(G$dag.old, layout=layout.circle, main = "Inferred old edges")
#' plot(G$dag.new, layout=layout.circle, main = "Inferred new edges")
#' par(old.par)
#'
SEMdag<- function(graph, data, LO="TO", beta=0, lambdas=NA, penalty=TRUE, verbose=FALSE, ...)
{
	# Set DAG objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	if (!is_dag(ig) & LO == "TO"){
	 cat("WARNING: input graph is not acyclic !\n")
	 cat(" Applying graph -> DAG conversion...\n")
	 dag<- graph2dag(ig, data) #del cycles & all <->
	}else{ dag<- ig }
	#X <- scale(data[,V(dag)$name])
	X<- as.matrix(data[,V(dag)$name])

	# Estimate DAG using linear ordering (LO) approach:
	x<- DAG_HD_TD(graph=dag, X=X, LO=LO, beta=beta,
	 lambdas=lambdas, penalty=penalty, verbose=verbose)
	colnames(x$adj)<- rownames(x$adj)<- colnames(X)
	
	# Mapping DAG edges on input graph:
	if (sum(x$adj) == 0) return(message("DAG with 0 edges, decrease beta threshold !"))
	ig1<- graph_from_adjacency_matrix(x$adj, mode="directed")
	ig2<- quiet(properties(ig1)[[1]])
	E1<- attr(E(ig2), "vnames")
	E0<- attr(E(ig), "vnames")
	E(ig2)$color<- ifelse(E1 %in% E0, "gray", "green")
	ig3<- ig2-E(ig2)[which(E(ig2)$color == "gray")]
	ig3<- ig3-vertices(V(ig3)$name[igraph::degree(ig3) == 0])
	ig4<- ig2-E(ig2)[which(E(ig2)$color == "green")]
	ig4<- ig4-vertices(V(ig4)$name[igraph::degree(ig4) == 0])
	
	if (verbose) {
	 gplot(ig2)
	 fit<- SEMrun(ig2, X, algo = "ricf", n_rep = 0)
	}
	
	return( list(dag=ig2, dag.new=ig3, dag.old=ig4) )
}

DAG_HD_TD<- function(graph, X, LO, beta, lambdas, penalty, verbose, ...)
{
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	if (LO == "TO"){
	 TO <- names(igraph::topo_sort(graph))
	 rr <- rev(match(TO, colnames(X)))
	 l <- sqrt(log(p)/n)
	}else if (LO == "TD"){
	 J <- 3
	 rr <- rev(getOrdering(X, J))
	 l <- sqrt(log(p)/n)
	 #l <- (2/sqrt(n))*qnorm(1 - 0.05/(2*p*(p-1)))
	}
	cat("Node Linear Ordering with", LO, "setting\n\n")
	if (verbose) {print(colnames(X)[rev(rr)]);cat("\n")}
	
	result <- matrix(0, p, p)
	X <- scale(X)
	A <- as_adj(graph, sparse = FALSE)
	for (ii in 1:(p - 1)) {
		now <- rr[ii]
		this <- sort(rr[(ii + 1):p])
		if (penalty) {
		 pw <- 1 - A[this,now]
		}else{
		 pw <- rep(1, length(this))
		}
		if (sum(pw) == 0) break
		if (length(this) > 1) {
			if (!is.null(lambdas)) {
				 lassom <- glmnet::glmnet(X[, this], X[, now], 
				  lambda = l, penalty.factor = pw)
				 bfit <- coefficients(lassom)[-1]
			}
			else{
				 lassom <- glmnet::cv.glmnet(X[, this], X[, now],
				  lambda = NULL, penalty.factor = pw)
				 bfit <- coefficients(lassom)[-1]
			}
			for (jj in 1:length(this)) {
				if (abs(bfit[jj]) > beta)
				result[this[jj], now] <- 1
			}
		}
		else {
			#lmod <- summary(RcppEigen::fastLm(X[, now] ~ X[, this]))
			lmod <- summary(lm(X[, now] ~ X[, this]))$coefficients
			if (lmod[2, 4] < 0.05 & abs(lmod[2, 1]) > beta) {
			#if (lmod$coef[2, 4] < 0.05) {
				result[this, now] <- 1
			}
		}
	}
	return(list(adj = result, TO = rev(rr)))
}

getOrdering<- function (Y, J, ...)
{
	# Copyright (c) 2020  Wenyu Chen [email ?]
	# https://github.com/WY-Chen/EqVarDAG
	# All rights reserved.  See the file COPYING for license terms.
	p <- dim(Y)[2]
	variances <- apply(Y, MARGIN = 2, sd)
	Theta <- rep(0, p)
	Theta[1] <- which.min(variances)
	out <- sapply(setdiff(1:p, Theta[1]), function(z) {
	  #sum(resid(RcppEigen::fastLm(Y[, z] ~ Y[, Theta[1], drop = F]))^2)
	  sum(resid(lm(Y[, z] ~ Y[, Theta[1], drop = F]))^2)
	 })
	Theta[2] <- setdiff(1:p, Theta[1])[which.min(out)]
	for (i in 3:p) {
	 out <- suppressWarnings(lapply(setdiff(1:p, Theta),
			function(jj) subsets(jj, Y, Theta[seq(i - 1)], J)))
	 nextRoot <- which.min(sapply(out, function(x) {
		min(x$rss)
	 }))
	 Theta[i] <- setdiff(1:p, Theta)[nextRoot]
	}
	return(Theta)
}

subsets<- function (z, Y, Theta, J, mtd = "seqrep") 
{
  #mtd = c("exhaustive","backward", "forward", "seqrep")
  leaps::regsubsets(x = Y[, Theta, drop = F], y = Y[, z, drop = F], 
                    method = mtd, nbest = 1, nvmax = min(J, sum(Theta > 0)), 
                    really.big = TRUE)
}

#' @title Interactome-assisted graph re-seizing
#'
#' @description An input directed graph is re-sized, removing edges
#' or adding edges/nodes. This function takes three input graphs: the
#' first is the input causal model (i.e., a directed graph), and the
#' second can be either a directed or undirected graph, providing a set
#' of connections to be checked against a directed reference network
#' (i.e., the third input) and imported to the first graph.
#'
#' @param g A list of two graphs as igraph objects, g=list(graph1, graph2).
#' @param gnet External directed network as an igraph object. The reference
#' network should have weighted edges, corresponding to their interaction
#' p-values, as an edge attribute \code{E(gnet)$pv}. Then, connections in
#' \code{graph2} will be checked by known connections from the reference network,
#' intercepted by the minimum-weighted shortest path found among the equivalent
#' ones by the Dijkstra algorithm, as implemented in the \pkg{igraph} function
#' \code{all_shortest_paths()}.
#' @param d An integer value indicating the maximum geodesic distance between
#' two nodes in the interactome to consider the inferred interaction between
#' the same two nodes in \code{graph2} as validated, otherwise the edges are
#' removed. For instance, if \code{d = 2}, two interacting nodes must either
#' share a direct interaction or being connected through at most one mediator
#' in the reference interactome (in general, at most \code{d - 1} mediators are
#' allowed). Typical \code{d} values include \code{2} (at most one mediator), or
#' \code{mean_distance(gnet)} (i.e., the average shortest path length for
#' the reference network). Setting d = 0, is equivalent to \code{gnet = NULL}.
#' @param v A logical value. If TRUE (default) new nodes and edges on the
#' validated shortest path in the reference interactome will be added in the
#' re-sized graph.
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen, saving execution time (for large graphs)
#' @param ... Currently ignored.
#'
#' @details Typically, the first graph is an estimated causal graph (DAG),
#' and the second graph is the output of \code{\link[SEMgraph]{SEMdag}} or
#' an external covariance graph. In the former we use the new inferred
#' causal structure stored in the \code{dag.new} object. In the latter, we
#' use the new inferred covariance structure stored in the covariance graph
#' object. Both directed (causal) edges and covariances (i.e., bidirected
#' edges) highlight emergent hidden topological proprieties, absent in the
#' input graph. Estimated directed edges between nodes X and Y
#' are interpreted as either direct links or direct paths mediated by hidden
#' connector nodes. Covariances between any two bow-free nodes X and Y may
#' hide causal relationships, not explicitly represented in the current model.
#' Conversely, directed (or bi-directed) edges could be redundant or artifact,
#' specific to the observed data and could be deleted.
#' Function \code{resizeGraph()} leverage on these concepts to extend/reduce a
#' causal model, importing new connectors or deleting estimated edges, if they are
#' present or absent in a given reference network. The whole process may lead to
#' the discovery of new paths of information flow, and cut edges not corroborate
#' by a validated network. Since added nodes can already be present in the causal
#' graph, network resize may create cross-connections between old and new paths
#' and their possible closure into circuits.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
#'
#' @return "Ug", the re-sized graph, the graph union of the causal graph \code{graph1}
#' and the re-sized graph \code{graph2}
#'
#' @examples
#'
#' \donttest{
#'
#' # Extract the "Protein processing in endoplasmic reticulum" pathway:
#'
#' g <- kegg.pathways[["Protein processing in endoplasmic reticulum"]]
#' G <- properties(g)[[1]]; summary(G)
#'
#' # Extend a graph using new inferred DAG edges (dag+dag.new):
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' dag <- SEMdag(graph = G, data = als.npn, beta = 0.1)
#' gplot(dag$dag)
#' ext <- resizeGraph(g=list(dag$dag, dag$dag.new), gnet = kegg, d = 2)
#' gplot(ext) 
#'
#' # Create a directed graph from correlation matrix, using
#' # i) an empty graph as causal graph,
#' # ii) a covariance graph,
#' # iii) KEGG as reference:
#'
#' corr2graph<- function(R, n, alpha=5e-6, ...)
#' {
#' 	Z <- qnorm(alpha/2, lower.tail=FALSE)
#'	thr <- (exp(2*Z/sqrt(n-3))-1)/(exp(2*Z/sqrt(n-3))+1)
#' 	A <- ifelse(abs(R) > thr, 1, 0)
#' 	diag(A) <- 0
#' 	return(graph_from_adjacency_matrix(A, mode="undirected"))
#' }
#'
#' v <- which(colnames(als.npn) %in% V(G)$name)
#' selectedData <- als.npn[, v]
#' G0 <- make_empty_graph(n = ncol(selectedData))
#' V(G0)$name <- colnames(selectedData)
#' G1 <- corr2graph(R = cor(selectedData), n= nrow(selectedData))
#' ext <- resizeGraph(g=list(G0, G1), gnet = kegg, d = 2, v = TRUE)
#' 
#' #Graphs
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2), mar=rep(1,4))
#' plot(G1, layout = layout.circle)
#' plot(ext, layout = layout.circle)
#' par(old.par)
#'
#' }
#'
resizeGraph<- function(g=list(), gnet, d = 2, v = TRUE, verbose = FALSE, ...)
{
	# Set graph objects (gnet, ig, guu):
	if (!is_directed(gnet)) {
	 return(message(" ERROR: Reference graph is NOT a directed graph !"))
	}
	if (!is_directed(g[[1]])) {
	 return(message(" ERROR: First input graph is NOT a directed graph !"))
	}
	ig<- g[[1]]
	if(!is.null(V(ig)$color)) ig<- delete_vertex_attr(ig, "color")
	if(!is.null(E(ig)$color)) ig<- delete_edge_attr(ig, "color")
	guu<- getNetEdges(g[[2]], gnet, d, yes = is.directed(g[[2]]))
	if (ecount(guu) == 0) {
	 message(" no edges u->u (or u--u) found !")
	 return(graph = ig)
	}
	
	if (v == TRUE) {
	 # all_shortest_paths calculates ALL shortest paths from=x to=y
	 # shortest_paths calculates a SINGLE shortest path from=x to=y
	 ftm<- as_edgelist(guu)
	 vpath<- ftmuv<-  NULL
	 for (i in 1:nrow(ftm)) {
		mode<- ifelse(is.directed(guu) & is.directed(gnet), "out", "all")
		if(distances(gnet,ftm[i,1],ftm[i,2],mode=mode,weights=NA) == Inf) next
		if(is.null(E(gnet)$pv)){
		suppressWarnings(path<- shortest_paths(gnet, ftm[i,1], ftm[i,2],
												mode = mode, weights = NA)$vpath)
		}else{
		suppressWarnings(path<- all_shortest_paths(gnet, ftm[i,1], ftm[i,2],
													mode = mode, weights = NA)$res)
		}
		if(length(path) >1){
		 fX2<- NULL
		 for( k in 1:length(path)) {
			pathk<- induced_subgraph(gnet, V(gnet)$name[path[[k]]])
			fX2[k]<- -2*sum(log(E(pathk)$pv))
		 }
		 path<- path[[which(fX2 == max(fX2))[1]]]
		}else{
		 path<- path[[1]]
		}
		V<- V(gnet)$name[path]
		vpath<- c(vpath, V[-c(1,length(V))])
		for(h in 1:(length(V)-1)) ftmuv<- rbind(ftmuv, c(V[h], V[h+1]))
	 }
	
	 ftmuv<- na.omit(ftmuv[duplicated(ftmuv) != TRUE,])
	 guv<- graph_from_edgelist(ftmuv, directed=is.directed(guu))
	 #guv<- guv - E(guv)[which_mutual(guv)]
	 vv<- V(guv)$name[-which(V(guv)$name %in% V(ig)$name)]
	 V(guv)$color[V(guv)$name %in% vv]<- "green"
	}else{
	 guv<- guu
	}
	
	if (!is.directed(guv)) guv <- orientEdges(guv, gnet)
	Ug<- graph.union(g=list(ig,guv))
	Ug<- quiet(properties(Ug)[[1]])
	E1 <- attr(E(Ug), "vnames")
	E0 <- attr(E(ig), "vnames")
	E(Ug)$color <- ifelse(E1 %in% E0, "gray", "green")
	if (verbose) gplot(Ug, main="Resized Graph (Ug)")
	
	return(graph = Ug)
}

getNetEdges<- function(graph, gnet, d, yes, ...) 
{
	# External validation of discovery edges from a reference interactome
	SET1<- as_edgelist(graph)
	if( nrow(SET1) == 0 ){
	 return(message("n.interactions of input graph = 0 !"))
	}	
	ftm1<- NULL
	for(j in 1:nrow(SET1)) { #j=14
	 cat("\r","edge set=", j, "of", nrow(SET1))
	 #Sys.sleep(0.01)
	 flush.console()
	 a<- SET1[j,1]
	 b<- SET1[j,2]
	 ftm1<- rbind(ftm1, c(a,b))
	 v<- which(V(gnet)$name %in% c(a,b))
	 if (length(v) == 2) {
		if(yes == FALSE) sp<- distances(gnet, a, b, mode="all", weights=NA)
		if(yes == TRUE) sp<- distances(gnet, a, b, mode="out", weights=NA)
		if(sp <= d) ftm1[j,]<- c(a,b) else ftm1[j,]<- c(NA,NA)
	 }else{ ftm1[j,]<- c(NA,NA) }
	} 
	#ftm1<- na.omit(ftm1)
	ftm1<- rbind(na.omit(ftm1), SET1[SET1[,1] == "group",])
	cat("\n\n", "n. edges to be evaluated:", nrow(SET1),
		"\n", "n. edges selected from interactome:", nrow(ftm1), "\n\n")
	guu<- graph_from_edgelist(ftm1, directed=yes) #gplot(guu)
	
	return(guu)
}

#' @title Tree-based structure learning methods
#'
#' @description Four tree-based structure learning methods are implemented
#' with graph and data-driven algorithms. A tree ia an acyclic graph with p vertices
#' and p-1 edges. The graph method refers to the  Steiner Tree (ST), a tree from an
#' undirected graph that connect "seed" with additional nodes in the
#' "most compact" way possible. The data-driven methods propose fast and scalable
#' procedures based on Chu-Liu–Edmonds’ algorithm (CLE) to recover a tree from a full
#' graph. The first method, called Causal Additive Trees (CAT) uses pairwise mutual
#' weights as input for CLE algorithm to recover a directed tree (an "arborescence").
#' The second one applies CLE algorithm for skeleton recovery and extends the skeleton
#' to a tree (a "polytree") represented by a Completed Partially Directed Acyclic Graph
#' (CPDAG). Finally, the Minimum Spanning Tree (MST) connecting an undirected graph (or
#' on an undirected full graph) with minimal edge weights can be identified.
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param seed A vector of seed nodes.  
#' @param type Tree-based structure learning method. Four algorithms 
#' are available:
#' \itemize{
#' \item "ST". Steiner Tree (ST) identification via fast Kou's algorithm 
#' (Kou et al, 1981) connecting a set of seed nodes (called Terminal vertices)
#' with connector nodes (called Steiner vertices) from input graph as defined
#' in \code{graph} with minimal total distance on its edges. By default the edge
#' weights are based on the pairwise correlation, 1-abs(cor(j,k)). If input
#' graph has E(graph)$weight=1, and \code{eweight = "custom"}, ST seeks a minimum
#' subtree (i.e. subtree with minimal number of edges).
#' \item "CAT" (default). Causal additive trees (CAT) algorithm as in Jakobsen et al. 
#' (2022). The argument \code{graph} is set to NULL (i.e., no input graph is needed).
#' In the first step, a (univariate) generalized additive model (GAM) is employed
#' to estimate the residual variances, var(X(j) - [X(j)|X(k)]) for all j != k,
#' then use these to construct edge weights as inputs to the Chu–Liu–Edmonds’
#' algorithm (Chow and Liu, 1968) to recover the arborescence. Argument \code{seed}
#' must be specified to analyse a subset of nodes (variables) of interest.
#' \item "CPDAG". CLE algorithm for Skeleton Recovery and CPDAG
#' estimation as in Lou et al. (2021). Together with "CAT" algorithm, "CPDAG" is 
#' data-driven and the argument \code{graph} is set to NULL.
#' The key idea is to first recover the skeleton of the polytree by applying 
#' the CLE algorithm to the pairwise sample correlations of the data matrix.
#' After the skeleton is recovered, the set of all v-structures can be correctly
#' identified via a simple thresholding approach to pairwise sample correlations.
#' CPDAG can be found applying iteratively only Rule 1 of Meek (1995).
#' Argument \code{seed} must be specified to analyse a subset of nodes
#' (variables) of interest.
#' \item "MST". Minimum Spanning Tree (MST) identification via Prim's algorithm
#' (Prim, 1957). The latter finds the subset of edges that includes every vertex
#' of the graph (as defined in \code{graph}) such that the sum of the weights 
#' of the edges can be minimized. The argument \code{seed} is set to NULL (i.e.,
#' no seed nodes are needed), or if argument \code{seed} is not NULL, argument
#' \code{graph} is set to NULL to recover the MST of the full graph.}
#' @param eweight Edge weight type for igraph object derived from
#' \code{\link[SEMgraph]{weightGraph}} or from user-defined distances. 
#' This option determines the weight-to-distance transform.
#' If set to "NULL" (default), edge weights will be internally computed
#' equal to 1 - abs(cor), i.e., 1 - abs(pairwise Pearson's correlation).
#' If \code{eweight = "kegg"}, repressing interactions (-1) will be set 
#' to 2 (maximum distance), neutral interactions (0) will be set to 1, 
#' and activating interactions (+1) will be set to 1 (minimum distance).
#' If \code{eweight = "zsign"}, all significant interactions will be set 
#' to 1 (minimum distance), while non-significant ones will be set to 2.
#' If \code{eweight = "pvalue"}, weights (p-values) will be transformed 
#' to the inverse of negative base-10 logarithm. 
#' If \code{eweight = "custom"}, the algorithm will use the distance 
#' measure specified by the user as "weight" edge attribute.
#' @param alpha Threshold for rejecting a pair of node being independent in 
#' "CPDAG" algorithm. The latter implements a natural v-structure identification 
#' procedure by thresholding the pairwise sample correlations over all adjacent 
#' pairs of edges with some appropriate threshold. By default, 
#' \code{alpha = 0.05}.
#' @param verbose If TRUE, it shows the output tree (not recommended for large graphs).
#' @param ... Currently ignored.
#'
#' @details If the input graph is a directed graph, ST and MST undirected trees are
#' converted in directed trees using the \code{\link[SEMgraph]{orientEdges}} function.
#'
#' @export
#'
#' @return An \code{igraph} object. If \code{type = "ST"}, seed nodes are 
#' colored in "aquamarine" and connectors in "white". If \code{type = "ST"} and
#' \code{type = "MST"}, edges are colored in "green" if not present in the input,
#' graph. If \code{type = "CPDAG"}, bidirected edges are colored in "black"
#' (if the algorithm is not able to establish the direction of the relationship
#' between x and y).
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#' 
#' @references
#'
#' Kou, L., Markowsky, G., Berman, L. (1981). A fast algorithm for Steiner trees. 
#' Acta Informatica 15, 141–145. <https://doi.org/10.1007/BF00288961>
#'
#' Prim, R.C. (1957). Shortest connection networks and some generalizations Bell
#' System Technical Journal, 37 1389–1401. 
#'
#' Chow, C.K. and Liu, C. (1968). Approximating discrete probability distributions with 
#' dependence trees. IEEE Transactions on Information Theory, 14(3):462–467.
#' 
#' Meek, C. (1995). Causal inference and causal explanation with background knowledge.
#' In Proceedings of the Eleventh conference on Uncertainty in artificial intelligence,
#' 403–410.
#'
#' Jakobsen, M, Shah, R., Bühlmann, P., Peters, J. (2022). 
#' Structure Learning for Directed Trees. arXiv:
#' <https://doi.org/10.48550/arxiv.2108.08871>.
#'
#' Lou, X., Hu, Y., Li, X. (2022). Linear Polytree Structural Equation Models:
#' Structural Learning and Inverse Correlation Estimation. arXiv:
#' <https://doi.org/10.48550/arxiv.2107.10955>
#'
#' @examples
#' 
#' \donttest{
#'
#' library(huge)
#' data <- huge.npn(alsData$exprs)
#' graph <- alsData$graph
#'
#' # graph-based trees
#' seed <- V(graph)$name[sample(1:vcount(graph), 10)]
#' tree1<- SEMtree(graph, data, seed=seed, type="ST", verbose=TRUE)
#' tree2<- SEMtree(graph, data, seed=NULL, type="MST", verbose=TRUE)
#'
#' # data-driven trees
#' V <- colnames(data)[colnames(data) %in% V(graph)$name]
#' tree3<- SEMtree(NULL, data, seed=V, type="CAT", verbose=TRUE)
#' tree4<- SEMtree(NULL, data, seed=V, type="CPDAG", alpha=0.05, verbose=TRUE)
#'
#' }
#'
SEMtree <- function(graph, data, seed, type = "CAT", eweight = NULL, alpha = 0.05, verbose = FALSE, ...)
{
	# Set data and graph objects:
	if (!is.null(graph)) {
	 nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	 ig <- induced_subgraph(graph, vids = V(graph)$name %in% nodes)
	 X <- data[,nodes]
	 if (!is.null(eweight)) {
	  if (eweight == "kegg") E(graph)$weight <- 2 - E(graph)$weight
	  else if (eweight == "zsign") E(graph)$weight <- 2 - abs(E(graph)$zsign)
	  else if (eweight == "pvalue") E(graph)$weight <- 1/(-log10(E(graph)$pv))
	  else if (eweight == "custom") E(graph)$weight <- E(graph)$weight
	 }
	 else if (is.null(eweight)) {
	  A <- (1-abs(cor(X)))*as_adj(ig, sparse=FALSE)[nodes,nodes]
	  d_u <- ifelse(is.directed(ig), "directed", "undirected")
	  graph <- graph_from_adjacency_matrix(A, mode=d_u, weighted=TRUE, diag=FALSE)
	 }
	 # SteinerTree(ST) or MinimumSpanningTree(MST):
	 if (!is.null(seed) & type == "ST") {
	   T <- SteinerTree(graph, seed, eweight = E(graph)$weight)
	 }else if(is.null(seed) & type == "MST") {
	   eattr <- list(weight="mean", "ignore")
	   ug <- as.undirected(graph, edge.attr.comb = eattr) 
	   T <- mst(ug, weights = E(ug)$weight, algorithm = "prim")
	 }
	 if (is.directed(graph) & !is.directed(T)) T <- orientEdges(ug=T, dg=graph) 
	 T <- quiet(properties(T)[[1]])
	 V(T)$color <- ifelse(V(T)$name %in% seed, "aquamarine", "white")
	 E1 <- attr(E(T), "vnames")
	 E0 <- attr(E(graph), "vnames")
	 E(T)$color <- ifelse(E1 %in% E0, "gray60", "green")
	 E(T)$color <- ifelse(which_mutual(T), "black", E(T)$color)
	 E(T)$width <- ifelse(which_mutual(T), 2, 1)

	# Causal Addittive Tree(CAT) or CPDAG Tree or MST:
	}else if(is.null(graph)) {
	  nodes <- colnames(data)[colnames(data) %in% seed]
	  X <- data[, nodes]
	  if (type == "CAT") T <- CAT.R(data = data.frame(X))
	  if (type == "CPDAG") T <- CPDAG(X, alpha, verbose=TRUE)
	  if (type == "MST") {
	   A <- 1-abs(cor(X))
	   gA <- graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
	   T <- mst(gA, weights = E(gA)$weight, algorithm = "prim")
	  }
	}

	if (verbose) {
	 gplot(T)
	 fit <- SEMrun(T, X, algo = "ricf", n_rep = 0)
	}
	
	return(Tree = T)
}

SteinerTree<- function(graph, seed, eweight, ...)
{ 
	# step 0) Define  distance matrix:
	seed<- seed[which(seed %in% V(graph)$name)]
	D<- distances(graph, v=seed, to=seed, mode="all", weights=eweight)
	 	 
	# step 1) Complete undirected distance graph Gd for terminal nodes:
	Gd<- graph_from_adjacency_matrix(D, mode="undirected", weighted=TRUE)
	Gd<- delete_edges(Gd, E(Gd)[which(E(Gd)$weight == Inf)])
	#Gd<- Gd-igraph::edges(E(Gd)[which(E(Gd)$weight == Inf)])

	# step 2) MST T1 of the complete distance graph Gd:
	T1<- mst(Gd, weights=NULL, algorithm="prim")

	# step 3) For each edge in T1, replace it with the shortest path in ig:
	edge_list<- as_edgelist(T1)
	N<- nrow(edge_list)
	subgraph<- vector()

	for (n in 1:N) { 
	 i <- edge_list[n,1]
	 j <- edge_list[n,2]
	 # extract from ig all nodes of the shortest paths between edges of T1:
	 path<- shortest_paths(graph, from=V(graph)[i], to=V(graph)[j], mode="all", weights=eweight, output="both")
	 vpath<- V(graph)$name[path$vpath[[1]]]
	 subgraph<- union(subgraph, vpath)	
	}

	# step 4) MST Ts of the extracted (induced) sub-graph Gs of ig:
	Gs<- induced_subgraph(as.undirected(graph), unique(subgraph))
	Ts<- mst(Gs, weights=NULL, algorithm="prim")

	# step 5) Pruning non-seed genes with degree=1 (one at time) from Ts:
	St<- Ts
	idx<- ifelse(V(St)$name %in% seed == TRUE, FALSE, TRUE)

	i<-1
	I<-length(V(St)[idx])+1
	while( i < I ) {
	 K<- igraph::degree(St, v=V(St), mode="all")
	 todel<- names(which(K == 1))
	 todel<- todel[which(!todel %in% seed)]
	 if( length(todel) > 0 ) {
	  St<- delete.vertices(St, todel)
	 }
	 i<-i+1
	}
	
	return(St)
}

CAT.R<- function(data, limit = 100, ...)
{   	
	local<- function(x){
	 form <- formula(paste0("X",x[[2]],"~","s(X",x[[1]],",bs='ts')"))
	 Cond_exp <- mgcv::gam(form, data = data)
	 #form <- formula(paste0("X",x[[2]],"~","X",x[[1]]))
	 #Cond_exp <- lm(form, data = data)
	 form <- paste0("data$X",x[[2]]," - predict(Cond_exp, newdata=data)")
	 Residual <- eval(parse(text=form))
	 var(Residual)
	}

	# Saving original column names and setting standard column names
	colNames <- gsub("X", "", colnames(data))
	colnames(data) <- paste0("X",seq(1,ncol(data),1))

	# Compute Gaussian Edge Weights:
	if (is.data.frame(data)) {
	 ig0 <- make_full_graph(ncol(data), directed = TRUE)
	 Edges <- as_data_frame(ig0)
	 x <- split(Edges, f = seq(nrow(Edges)))
	 message("Score weighting of ", length(x), " edges...")
	 op <- pbapply::pboptions(type = "timer", style = 2)

	 if (ncol(data) > limit){
	  n_cores <- parallel::detectCores(logical = FALSE)
	  cl <- parallel::makeCluster(n_cores)
	  parallel::clusterExport(cl, c("local"),
	   envir = environment())
	  Vr <- pbapply::pblapply(x, local, cl=cl)
	  parallel::stopCluster(cl)
	 }else{
	  Vr <- pbapply::pblapply(x, local, cl=NULL)
	 }

	 Vr <- cbind(Edges, weight=do.call(rbind,Vr))
	 Vx <- apply(data, 2, var)
	 gW <- graph_from_data_frame(Vr)
	 W <- as_adj(gW, attr="weight", sparse=FALSE)
	 diag(W) <- Vx
	}else{
	 Vx <- rep(1, ncol(data))
	 W <- 1 - data^2
	}

	# Chu-Liu-Edmonds's Algorithm
	NegW <- -log(W%*%diag(1/Vx))
	NegW <- NegW + 2*abs(min(NegW))
	rownames(NegW) <- colnames(NegW) <- colNames
	gx <- graph_from_adjacency_matrix(NegW, mode="directed", weighted=TRUE, diag=FALSE)
	ax <- RBGL::edmondsOptimumBranching(as_graphnel(gx))
	ax <- data.frame(t(ax$edgeList),t(ax$weight))
	OptimalTree <- graph_from_data_frame(ax)

	return(OptimalTree)
}

CPDAG <- function(data, alpha, verbose = FALSE, ...)
{
	# a) MST on the full graph (with zero edges)
	R<- cor(data)
	MI<- -nrow(data)*log(1-R^2 + 1e-16)
	A<- ifelse(MI <= qchisq(1-alpha, df=1), 0, abs(R))
	gA<- graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
	MST<- mst(gA, weights = 1-E(gA)$weight, algorithm = "prim")
			
	# b) CPDAG recovery
	g <- as_adjacency_matrix(MST, sparse=FALSE)
	g <- ifelse(abs(g) != 0, 1, 0)
	pdag <- g
	ind <- which(g == 1, arr.ind = TRUE)
	
	# 1) v-structures for all node triplets i--k--j
	for (i in seq_len(nrow(ind))) {
	 x <- ind[i, 1]
	 y <- ind[i, 2]
	 allZ <- setdiff(which(g[y, ] == 1), x)
	 for (z in allZ) {
		if (g[x, z] == 0 & A[x, z] == 0) {
		 if (verbose) {
		  cat("\n", x, "->", y, "<-", z, "\n")
		}
		 pdag[x, y] <- pdag[z, y] <- 1
		 pdag[y, x] <- pdag[y, z] <- 0
		}
	 }
	}
	if (sum(g) == sum(pdag)) {
	 message(" WARNING: none v-structures are recovery, CPDAG=MST !\n")
	 return(MST)
	}
	
	# 2) apply Rule 1 in the resulting PDAG:
	dagy <- graph2dagitty(pdag, graphType = "pdag")
	# plot( dagitty::graphLayout(dagy) ) 
	cpdag <- dagitty::orientPDAG(dagy)
	# plot( dagitty::graphLayout(cpdag))
	CPDAG <- dagitty2graph(cpdag)
	# gplot(CPDAG)
	CPDAG <- quiet(properties(CPDAG)[[1]])
	E(CPDAG)$color <- ifelse(which_mutual(CPDAG), "black", "gray60")
	E(CPDAG)$width <- ifelse(which_mutual(CPDAG), 2, 1)
	
	return(CPDAG)
}

#' @title Optimal model search strategies
#'
#' @description Four model search strategies are implemented combining
#' \code{SEMdag()}, \code{SEMbap()}, and \code{resizeGraph()} functions.
#' All strategies estimate a new graph by 1) adjusting (BAP deconfounding) the
#' the data matrix and 2) re-sizing the output DAG.
#' @param graph Input graph as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param gnet Reference directed network used to validate and import
#' nodes and interactions.
#' @param d Maximum allowed geodesic distance for directed or undirected
#' shortest path search. A distance \code{d = 0} disables shortest path
#' search (fixed in \code{search = "basic"}), while \code{d = 1} 
#' (fixed in \code{search = "direct"}) only search for directed links
#' (i.e., no mediators are allowed).
#' A distance \code{d} > 1 (defaults to \code{d = 2} for "outer" and
#' "inner" strategies), will search for shortest paths with at most
#' \code{d} - 1 mediators between nodes sharing a significant estimated
#' interaction.
#' Connectors are imported from the reference interactome, as specified
#' by the argument \code{gnet}. If the edges of the reference interactome are
#' weighted by P-value, as defined by the \code{E(gnet)$pv} attribute,
#' the shortest path with the smallest sum of weights will be chosen (e.g.,
#' see \code{\link[SEMgraph]{weightGraph}} for graph weighting options).
#' @param search Search strategy. Four model search strategies are available:
#' \itemize{
#' \item "outer". The estimated DAG is re-sized using
#' \code{\link[SEMgraph]{resizeGraph}} to find new indirect paths (i.e.,
#' inferred directed connections that may hide new mediators). New
#' interactions and connectors will be searched and imported from the
#' reference network (argument \code{gnet}, see above). Both DAG and
#' extended graph complexity can be controlled with \code{beta} > 0 and
#' \code{d} > 1 arguments, respectively. The term "outer" means that new
#' model mediator variables are imported from an external resource (i.e.,
#' the reference network).
#' \item "inner". This strategy is analogous to the "outer" one,
#' but disables external mediator search. In other words, new indirect
#' paths are generated by adding new interactions of the input model, so
#' that mediators will be nodes already present in the input graph. The
#' reference network is still used to validate new model paths. Also in
#' this case, \code{beta} > 0 and \code{d} > 1 are used.
#' \item "direct". The input graph structure is improved through direct
#' (i.e., adjacent) link search, followed by interaction validation and
#' import from the reference network, with no mediators (i.e., \code{d = 1}).
#' \item "basic" (default). While the previous strategies rely on the
#' input graph and the reference network to integrate knowledge to the
#' final model, the "basic" strategy is data-driven. The input graph is
#' needed to define the topological order. The argument \code{gnet} is
#' set to NULL (i.e., no reference network is needed) and argument
#' \code{d = 0}. Model complexity can be still controlled by setting
#' \code{beta} > 0.
#' }
#' @param beta Numeric value. Minimum absolute LASSO beta coefficient for
#' a new interaction to be retained in the estimated DAG backbone. Lower
#' \code{beta} values correspond to more complex DAGs.
#' By default, \code{beta} is set to 0 (i.e., maximum complexity).
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, \code{method}
#' is set to "BH" (i.e., Benjamini-Hochberg multiple test correction).
#' @param alpha Significance level for false discovery rate (FDR) used
#' for local d-separation tests. This argument is used to
#' control data de-correlation. A higher \code{alpha} level includes more
#' hidden covariances, thus considering more sources of confounding.
#' If \code{alpha = 0}, data de-correlation is disabled.
#' By default, \code{alpha = 0.05}.
#' @param limit An integer value corresponding to the number of missing
#' edges of the extracted acyclic graph. Beyond this limit, multicore
#' computation is enabled to reduce the computational burden.
#' @param verbose If TRUE, it shows intermediate graphs during the
#' execution (not recommended for large graphs).
#' @param ... Currently ignored.
#'
#' @details Search strategies can be ordered by decreasing conservativeness
#' respect to the input graph, as: "direct", "inner", "outer", and "basic".
#' The first three strategies are knowledge-based, since they require an
#' input graph and a reference network, together with data, for
#' knowledge-assisted model improvement. The last one does not require
#' any reference and the output model structure will be data-driven.
#' Output model complexity can be limited using arguments \code{d} and
#' \code{beta}.
#' While d is fixed to 0 or 1 in "basic" or "direct", respectively;
#' we suggest starting with \code{d = 2} (only one mediator)
#' for the other two strategies.
#' For knowledge-based strategies, we suggest to to start with
#' \code{beta = 0}. Then, beta can be relaxed (0 to < 0.1) to improve
#' model fitting, if needed. Since data-driven models can be complex,
#' we suggest to start from \code{beta = 0} when using the "basic" strategy.
#' The \code{beta} value can be relaxed until a good model fit is obtained.
#' Argument alpha determines the extent of data adjustment: lower alpha
#' values for FDR correction correspond to a smaller number of significant
#' confounding factors, hence a weaker correction
#' (default \code{alpha = 0.05}).
#'
#' @export
#'
#' @return The output model as well as the adjusted dataset are returned
#' as a list of 2 objects:
#' \itemize{
#' \item "graph", the output model as an igraph object;
#' \item "data", the adjusted dataset.
#' }
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' \donttest{
#'
#' # Comparison among different model estimation strategies
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' # Models estimation
#' m1 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
#'       search = "direct", beta = 0, alpha = 0.05)
#' m2 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
#'       d = 2, search = "inner", beta = 0, alpha = 0.05)
#' m3 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
#'       d = 2, search = "outer", beta = 0, alpha = 0.05)
#' m4 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = NULL,
#'       search = "basic", beta = 0.1, alpha = 0.05)
#'
#' # Graphs
#' #old.par <- par(no.readonly = TRUE)
#' #par(mfrow=c(2,2), mar= rep(1,4))
#' gplot(m1$graph, main = "direct graph")
#' gplot(m2$graph, main = "inner graph")
#' gplot(m3$graph, main = "outer graph")
#' gplot(m4$graph, main = "basic graph")
#' #par(old.par)
#'
#' }
#'
modelSearch<- function(graph, data, gnet = NULL, d = 2, search = "basic",
                       beta = 0, method = "BH", alpha = 0.05,
                       limit = 30000, verbose = FALSE, ...)
{
	# Step by step search:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	Zt<- as.matrix(data[,nodes])
	Gt<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	cat("Step1: BAP deconfounding...\n")
	Zt1<- quiet(SEMbap(Gt, Zt, method=method, alpha=alpha, limit=limit))
	cat("Step2: DAG estimation...\n")
	if(is.null(Zt1)) {
	 Zt1$data <- Zt
	 Gt1<- quiet(SEMdag(Gt, Zt1$data, beta=beta)$dag)
	}else{
	 Gt1<- quiet(SEMdag(Zt1$dag, Zt1$data, beta=beta)$dag)
	}
	if (is.null(Gt1)) return(NULL)
	E1<- attr(E(Gt1), "vnames")
	E0<- attr(E(Gt), "vnames")
	E(Gt1)$color<- ifelse(E1 %in% E0, "gray", "green")
	Gt1.new<- Gt1 - E(Gt1)[which(E(Gt1)$color == "gray")]
	if(ecount(Gt1.new) == 0){
	 return(message("DAG with 0 new edges, decrease beta threshold !"))
	}
	#gplot(Gt1); gplot(Gt1.new); plot(Zt1$guu); plot(Zt1$dag)
	
	cat("Step3: DAG resize (remove edges/add nodes)...\n\n")
	if (search == "basic"){
	 cat("None DAG resize for basic search !", "\n\n")
	 Gt2<- Gt1
	 dataZ<- Zt1$data
	}
	if (search == "direct"){
	 Gt2<- resizeGraph(g=list(Gt,Gt1.new), gnet, d=1, v=FALSE, verbose=FALSE)
	 dataZ<- Zt1$data
	}
	if (search == "inner") {
	 Gt2<- resizeGraph(g=list(Gt,Gt1.new), gnet, d=d, v=FALSE, verbose=FALSE)
	 dataZ<- Zt1$data
	}
	if (search == "outer") {
	 Gt2<- resizeGraph(g=list(Gt,Gt1.new), gnet, d=d, verbose=FALSE)
	 green<- V(Gt2)$name[V(Gt2)$color == "green"]
	 dataZ<- cbind(Zt1$data, data[,which(nodes %in% green)])
	}
	cat("Done.\n\n")
	C_test<- Shipley.test(Gt2, dataZ, verbose=TRUE)
	if (verbose) {
	 gplot(Gt2, main="Estimated Extended Graph")
	 cat("\n")
	 fit<- SEMrun(Gt2, dataZ, algo = "ricf", n_rep = 0)
	}
	
	return(list(graph = Gt2, data = dataZ))
}
