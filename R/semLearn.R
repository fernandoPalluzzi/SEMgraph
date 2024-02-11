#  SEMgraph library
#  Copyright (C) 2019-2024 Mario Grassi; Fernando Palluzzi; Barbara Tarantino 
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
#' the data matrix, Y or (iii) spectral transformation of Y. 
#' 
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and columns to
#' graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If \code{NULL} (default), confouding within group
#' will not be considered.
#' @param dalgo Deconfounding method. Four algorithms are available: 
#' \itemize{
#' \item "cggm" (default). The algorithm make: (i) exhaustive search
#' of bow-free significant covariances (see details) through  
#' \code{\link[SEMgraph]{Shipley.test}} function; (ii) estimation of the inverse
#' of the selected covariance matrix (i.e. the precision matrix, W) through
#' \code{\link[ggm]{fitConGraph}} function; (iii) obtain the de-correlated data
#' matrix, Z by multiplying the data matrix, Y rightward by the square root of
#' the estimated precision matrix, Z=YW^(1/2) as suggested by Grassi, Palluzzi
#' and Tarantino (2022).
#' \item "glpc". The algorithm first makes an exhaustive search of bow-free
#' significant covariances through \code{\link[SEMgraph]{Shipley.test}} function.
#' Once obtained the adjacency matrix, Graph-Laplacian PCA (gLPCA) algorithm
#' (Jiang et al., 2013) learns a low dimensional representation of the observed
#' data matrix that incorporates bow-free structure. Then, the DAG is extended
#' by including the confounding proxies, i.e. LVs, as additional source nodes
#' defined by last q principal component scores of gLPCA and these LV scores are
#' added to the data matrix, Z=cbind(LV,Y).
#' \item "pc". The procedure add additional source nodes to DAG as in "glpc"
#' algorithm, but confounding proxies are the q principal component scores
#' extracted by Spectral decomposition (SVD) selecting only graph nodes and
#' without graph edge information and bow-free covariance search.
#' \item "trim". Ćevid et al. (2020) suggest multiplying the data
#' matrix, Y leftward by a well selected spectrum transformation matrix, T
#' which modifies the singular values of Y, while keeping its singular vectors
#' intact, Z=TY. Trim transform limits all singular values to be at most some
#' costant (t), where t = median of the singuar values.
#' }
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, \code{method}
#' is set to "BH" (i.e., Benjamini-Hochberg multiple test correction).
#' @param alpha Significance level for false discovery rate (FDR) used
#' for d-separation test. This argument is used to
#' control data de-correlation. A higher \code{alpha} level includes more
#' hidden covariances, thus considering more sources of confounding.
#' If \code{alpha = 0}, data de-correlation is disabled. By default,
#' \code{alpha = 0.05}.
#' @param hcount The number of latent (or hidden) variables. By default
#' \code{hcount="auto"}, the hidden count is determined with a
#' permutation method (see details). Currently ignored if (dalgo ="cggm"
#' or "trim").
#' @param cmax Maximum number of parents set, C. This parameter can be
#' used to perform only those tests where the number of conditioning
#' variables does not exceed the given value. High-dimensional conditional
#' independence tests can be very unreliable. By default, cmax = Inf.
#' @param limit An integer value corresponding to the graph size (vcount)
#' tolerance. Beyond this limit, the precision matrix is estimated by
#' "glasso" algorithm (FHT, 2008) to reduce the computational burden of the
#' exaustive BAP search of the \code{\link[SEMgraph]{Shipley.test}} procedure.
#' By default, \code{limit = 200}.
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
#' of Y1 and Y2 (Shipley, 2000).
#' A new bow-free covariance is added if there is a significant (Y1, Y2)
#' association at a significance level \code{alpha}, after multiple testing
#' correction. The selected covariance between pairs of nodes is interpreted
#' as the effect of a latent variable (LV) acting on both nodes; i.e., the LV
#' is an unobserved confounder. BAP-based algorithms adjust (or de-correlate)
#' the observed data matrix by conditioning out the latent triggers responsible
#' for the nuisance edges.
#' For "pc" algorithm the number of hidden proxies, q is determined by a permutation
#' method. It compares the singular values to what they would be if the variables
#' were independent, which is estimated by permuting the columns of the data matrix,
#' Y and selects components if their singular values are larger than those of the
#' permuted data (for a review see Dobriban, 2020).
#' While for "glpc" algorithm, q is determined by the number of clusters by
#' spectral clustering through \code{\link[igraph]{cluster_leading_eigen}} function.
#' If the input graph is not acyclic, a warning message will be raised, and a
#' cycle-breaking algorithm will be applied (see \code{\link[SEMgraph]{graph2dag}}
#' for details).
#'
#' @return A list of four objects:
#' \itemize{
#' \item "dag", the directed acyclic graph (DAG) extracted from input graph.
#' If (dalgo = "glpc" or "pc"), the DAG also includes LVs as source nodes.
#' \item "guu", the bow-free covariance graph, BAP = dag + guu. If (dalgo =
#' "pc" or "trim"), guu is equal to NULL
#' \item "adj", the adjacency matrix of selected bow-free covariances; i.e, the 
#' missing edges selected after multiple testing correction. If (dalgo = "pc"
#' or "trim"), adj matrix is equal to NULL.
#' \item "data", the adjusted (de-correlated) data matrix or if (dalgo = "glpc",
#' or "pc"), the combined data matrix, where the first columns represent LVs
#' scores and the other columns are the raw data.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38(20), 4829–4830.
#' <https://doi.org/10.1093/bioinformatics/btac567>
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs.
#' Structural Equation Modeling, 7(2), 206-218.
#' <https://doi.org/10.1207/S15328007SEM0702_4>
#'
#' Jiang B, Ding C, Bin L, Tang J (2013). Graph-Laplacian PCA: 
#' Closed-Form Solution and Robustness. IEEE Conference on Computer
#' Vision and Pattern Recognition, 3492-3498. 
#' <https://doi.org/10.1109/CVPR.2013.448>
#' 
#' Ćevid D,  Bühlmann P, Meinshausen N (2020). Spectral deconfounding via
#' perturbed sparse linear models. J. Mach. Learn. Res, 21(232), 1-41.
#' <http://jmlr.org/papers/v21/19-545.html>
#'
#' Dobriban E (2020). Permuatation methods for Factor Analysis and PCA.
#' Ann. Statist. 48(5): 2824-2847
#' <https://doi.org/10.1214/19-AOS1907>
#'
#' Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse covariance
#' estimation with the graphical lasso. Biostatistics, 9(3), 432-441.
#' <https://doi.org/10.1093/biostatistics/kxm045>
#'
#' @examples
#'
#' #Set function param
#' graph <- sachs$graph
#' data <- log(sachs$pkc)
#' group <-sachs$group
#'
#' # BAP decounfounding with CGGM (default)
#' bap <- SEMbap(graph, data, verbose = TRUE)
#'
#' # SVD decounfounding with trim method
#' svd <- SEMbap(graph, data, dalgo = "trim")
#'
#' # Model fitting (with node-perturbation)
#' sem1 <- SEMrun(graph, data, group)
#' bap1 <- SEMrun(bap$dag, bap$data, group)
#' svd1 <- SEMrun(svd$dag, svd$data, group)
#'
SEMbap <- function(graph, data, group = NULL, dalgo = "cggm",
					method = "BH", alpha = 0.05, hcount = "auto",
					cmax = Inf, limit = 200, verbose = FALSE, ...)
{
	# Set graph and data objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	dag<- graph2dag(graph, data, bap=FALSE) #del cycles & all <->
	dataY<- data[, V(dag)$name]
		
	if (dalgo == "cggm" | dalgo == "glpc") {
	 if (vcount(dag) <= limit) {
	  cat("Bow-free covariances search. Use method:", dalgo, "...\n")
	  Z <- estimatePSI(dag, dataY, group, dalgo, method, alpha, cmax)
	 }else{
	  cat("Bow-free covariances search. Use method: glasso ...\n")
	  Z <- estimateGGM(dag, dataY, group, dalgo)
	 }
	 if (is.null(Z$guu)) {
	  return(message("NULL covariance graph, data de-correlation is not performed !"))
	 }else{
	  guu<- Z$guu; adj<- Z$adj; Z<- Z$data
	  cls <- length(cluster_leading_eigen(guu))
	  #cls <- count_components(guu)
	  cat("Number of clusters / number of nodes:", cls, "/", vcount(guu), "\n\n")
	}

	 # Covariance and latent variables graphs (guu, gLV)
	 if (verbose) {
	  ftm<- as_edgelist(as.undirected(guu))
	  ftmLV<- NULL
	  V(guu)$color<- "white"
	  for (i in 1:nrow(ftm)) ftmLV<- rbind(ftmLV, cbind(rep(paste0("L",i), 2),ftm[i,]))
	  gLV<- graph_from_data_frame(ftmLV, directed=TRUE)
	  V(gLV)$color<- ifelse(substr(V(gLV)$name,1,1)=="L","yellow","white")
	 
	  plot(guu, main="bow-free covariance graph (guu)")
	  Sys.sleep(3)
	  plot(gLV, main="bow-free latent variables graph (gLV)")
	  #Sys.sleep(3)
	  #gplot(graph.union(g=list(dag,gLV)), main="BAP graph (dag+gLV)")
	 }
	 
	 if (dalgo == "glpc") {
	  beta <- ifelse(cls > 3, 1, 0.75)
	  Z <- estimateLV(adj, dataY, group, hcount=cls, beta=beta)
	  dag <- map_hidden_dag(dag, Z)
	 }
	}

	else if (dalgo == "pc" | dalgo == "fa"){
	 Z <- estimateFX(dataY, group, dalgo=dalgo, hcount=hcount, verbose=verbose)
	 dag <- map_hidden_dag(dag, Z)
	 guu <- adj<- NULL
	}else{
	 Z <- estimateFX(dataY, group, dalgo=dalgo, hcount=hcount, verbose=verbose)
	 guu <- adj<- NULL
	}

	# SEM fitting with adjusted bow-free covariances:
	if (verbose) fit<- SEMrun(dag, Z, algo = "ricf", n_rep = 0)

	return( list(dag=dag, guu=guu, adj=adj, data=Z) )
}

estimatePSI <- function(dag, data, group, dalgo, method, alpha, cmax, ...)
{
	# Set data objects
	if (is.null(group)){
	 Y <- data
	}else{
	 Y_0 <- data[group == 0, ]
	 Y_1 <- data[group == 1, ]
	}
	
	estimate_ggm <- function(dag, X, dalgo, method, alpha, cmax, ...)
	{
	  # Set the bow-free covariance graph
	  dsep <- quiet(Shipley.test(dag, X, cmax=cmax, limit=100, verbose=FALSE)$dsep)
	  d_sep <- subset(dsep, p.adjust(dsep[,5 ], method) < alpha)
	  guu <- graph_from_data_frame(d_sep[,c(1,2)], directed = FALSE)
	  if (ecount(guu) > 0) {
	   adj <- as_adj(guu, type = "both", attr = NULL, sparse = FALSE)
	   s <- round(sqrt(nrow(X))/log(ncol(X)))
	   cat("Number of bow-free covariances / df :", nrow(d_sep),"/",nrow(dsep),"\n")
	   cat("Max parent set(S) / Sparsity idx(s) :", max(dsep[, 3]),"/",s,"\n")
	  }else{
	   return(NULL)
	  }

	  # Set the complete adjacency matrix
	  idx <- which(V(dag)$name %in% rownames(adj) == FALSE)
	  if (length(idx) > 0) {
	    R <- matrix(0, length(idx), ncol(adj))
	    C <- matrix(0, nrow(adj), length(idx))
	    D <- matrix(0, length(idx), length(idx))
	    adj <- rbind(cbind(D,R), cbind(C,adj))
	    rownames(adj)[1:length(idx)] <- V(dag)$name[idx]
	    colnames(adj)[1:length(idx)] <- V(dag)$name[idx]
	  }
      adj<- adj[colnames(X), colnames(X)]
	  
	  if (dalgo == "glpc") return(list(guu = guu, adj = adj, Z = X))
	 
	  # Fitting a concentration graph model with HTF procedure
	  suppressWarnings(
	   w <- tryCatch(ggm::fitConGraph(adj, S=cor(X), n=nrow(X))$Shat,
                     error = function(err) cor(X)*adj)
	  )
	  colnames(w) <- rownames(w) <- colnames(X)
	  	  
	  if (!corpcor::is.positive.definite(w)){
	   w <- corpcor::cov.shrink(w, verbose = FALSE)
	   #w<- corpcor::cor.shrink(w, verbose = TRUE)
	   #w<- corpcor::make.positive.definite(w)
	  }
	  E<- eigen(w) # Eigenvalues-eigenvectors of w
	  R<- E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
	  #sum(solve(w) - R %*% R)
	  Z <- as.matrix(X)%*%R
	  colnames(Z) <- colnames(X)
	  
	  return( list(guu = guu, adj = adj, Z = Z) )
	}

	if (is.null(group)){
     Z <- estimate_ggm(dag, Y, dalgo, method, alpha, cmax)
	}else{
	 Z_0 <- estimate_ggm(dag, Y_0, dalgo, method, alpha, cmax)
	 Z_1 <- estimate_ggm(dag, Y_1, dalgo, method, alpha, cmax)
	 group<- c(rep(0,nrow(Z_0$Z)),rep(1,nrow(Z_1$Z)))
	 Z01 <- cbind(group, rbind(Z_0$Z, Z_1$Z))
	 guu <- Z_0$guu %u% Z_1$guu
	 adj <- as_adj(guu, type = "both", sparse = FALSE)
	 Z <- list(guu=guu, adj=adj, Z=Z01)
	}

	return(list(guu = Z$guu, adj = Z$adj, data = Z$Z))
}

estimateGGM <- function(dag, data, group, dalgo, ...)
{
	# Set data objects
	if (is.null(group)){
	 Y <- data
	}else{
	 Y_0 <- data[group == 0, ]
	 Y_1 <- data[group == 1, ]
	}
	
	estimate_glasso <- function(dag, X, dalgo, ...)
	{
	  # graphical LASSO of large-scale data, rho=sqrt(log(p)/n)
	  rho <- sqrt(log(ncol(X))/nrow(X))
	  ug <- as.undirected(dag, mode="collapse")
	  V(ug)$name <- 1:vcount(ug)
	  ftm <- igraph::as_data_frame(ug)[1:2]
	  zero <- cbind(as.numeric(ftm[,1]),as.numeric(ftm[,2]))
	  wi <- glasso::glasso(s=cor(X), rho=rho, zero=zero)$wi
	  rownames(wi) <- colnames(wi) <- colnames(X)
	  
	  ai <- ifelse(abs(wi) > 0, 1, 0)
	  guu <- graph_from_adjacency_matrix(ai, mode="undirected", diag=FALSE)
	  if (ecount(guu) > 0) {
	   df0 <- vcount(dag)*(vcount(dag)-1)/2 - ecount(dag)
	   s <- round(sqrt(nrow(X))/log(ncol(X)))
	   d <- max(igraph::degree(guu))
	   cat("Number of bow-free covariances / df :", ecount(guu), "/", df0, "\n")
	   cat("Max node degree(d) / Sparsity idx(s):", d,"/",s,"\n")
	  }else{
	   return(NULL)
	  }
	  
	  if (dalgo == "glpc") return(list(guu = guu, adj = ai, Z = X))

	  # Eigenvalues-eigenvectors of wi
	  E <- eigen(wi*ai)
	  R <- E$vectors%*%diag(sqrt(E$values))%*%t(E$vectors)
	  #sum(wi - R %*% R)
	  Z <- as.matrix(X)%*%R
	  colnames(Z) <- colnames(X)
	  
	  return( list(guu = guu, adj = ai, Z = Z) )
	}

	if (is.null(group)){
     Z <- estimate_glasso(dag, Y, dalgo)
	}else{
	 Z_0 <- estimate_glasso(dag, Y_0, dalgo)
	 Z_1 <- estimate_glasso(dag, Y_1, dalgo)
	 group<- c(rep(0,nrow(Z_0$Z)),rep(1,nrow(Z_1$Z)))
	 Z01 <- cbind(group, rbind(Z_0$Z, Z_1$Z))
	 guu <- Z_0$guu %u% Z_1$guu
	 adj <- as_adj(guu, type = "both", sparse = FALSE)
	 Z <- list(guu=guu, adj=adj, Z=Z01)
	}

	return(list(guu = Z$guu, adj = Z$adj, data = Z$Z))
}

estimateLV <- function(adj, data, group, hcount, beta, ...)
{
	# Set data objects
	if (is.null(group)){
	 Y <- data[, colnames(adj)]
	}else{
	 Y_0 <- data[group == 0, colnames(adj)]
	 Y_1 <- data[group == 1, colnames(adj)]
	}

	# Estimate latent variables via gLPCA

	estimate_LV <- function(X, q, b)
	{
	 X <- scale(X)
	 n <- nrow(X)
	 p <- ncol(X)

	 S <- cov(X)*(n-1)/n
	 E <- eigen(S)
	 e1 <- E$values[1]
	 
	 adj<- S*adj #signed Adjacency
	 dj <- rowSums(abs(adj), na.rm = FALSE)
	 D <- diag(dj)
	 L <- D - adj #signed Laplacian
	 E <- eigen(L)
	 f1 <- E$values[1]

	 one <- rep(1,p)
	 I <- diag(one)
	 if (sum(D) == 0) {
	  b <- 0; f1 <- 1
	 }

	 G <- (1-b)*(I-S/e1) + b*(L/f1 + one%*%t(one)/p)
	 E <- svd(G, nu = 0)
	 W <- E$v #sum(round(t(W)%*%W, 6))
	 U <- X%*%W

	 #ev<- sort(1-E$d^2, decreasing = TRUE)
	 #q <- screePlot(ev, method="knee")
	 LV <- scale(U[,(p-q+1):p])
	 colnames(LV) <- paste0("LV", seq_len(q))
	 return(cbind(LV,X))
	}
	
	if (is.null(group)){
	 Z <- estimate_LV(Y, q=hcount, b=beta)
	}else{
	 Z_0 <- estimate_LV(Y_0, q=hcount, b=beta)
	 Z_1 <- estimate_LV(Y_1, q=hcount, b=beta)
	 group<- c(rep(0,nrow(Z_0)),rep(1,nrow(Z_1)))
	 Z <- cbind(group, rbind(Z_0, Z_1))
	}

	return(data = Z)
}

estimateFX <- function(data, group, dalgo, hcount, verbose, ...)
{	
	# Set data objects
	if (is.null(group)){
	 Y <- data
	 if (is.character(hcount)) {
	  hcount<- estimate_latent_count(Y, method=hcount, verbose=verbose)
	 }
	}else{
	 Y_0 <- data[group == 0, ]
	 Y_1 <- data[group == 1, ]
	 if (is.character(hcount)) {
	  K_0 <- estimate_latent_count(Y_0, method=hcount, verbose=verbose)
	  K_1 <- estimate_latent_count(Y_1, method=hcount, verbose=verbose)
	  hcount <- (K_0+K_1)/2
	 }
	}

	# Estimate deconfounding data via svd(X, nu, nv)

	estimate_FX <- function(X, q)
	{
	 X <- scale(X)
	 n <- nrow(X)
	 #p <- ncol(X)
	 
	 if (dalgo == "pc") {
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  LV <- sqrt(n-1)*as.matrix(x$u[,1:q]) #round(t(LV)%*%LV,6)
	  colnames(LV) <- paste0("LV", seq_len(q))
	  return(cbind(LV,X))
	 }

	 if (dalgo == "fa") {
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  LV <- factor.analysis(X, r=q, method="ml")$Z
	  LV <- scale (LV)
	  colnames(LV) <- paste0("LV", seq_len(q))
	  return(cbind(LV,X))
	 }

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
	  #E <- eigen(cov(X))
      #V <- E$vectors[,1:q]
      #Y<- X%*%(diag(p)-V%*%t(V))
      #colnames(Y)<- colnames(X)
      #return(Y)
      #if (p < n) {
      # V <- svd(X, nu = 0, nv = q)$v # round(t(V)%*%V,6)
      #  Y <- X%*%(diag(p)-V%*%t(V))
      #}else{
      #U <- svd(X, nu = q, nv = 0)$u # round(t(U)%*%U,6)
      # Y <- (diag(n)-U%*%t(U))%*%X
      #}
	  #colnames(Y) <- colnames(X)
      #return(Y)
	  x <- svd(X, nv = 0)
	  if (q == 0) return(X)
	  U <- sqrt(n)*x$u #round(t(U)%*%U,6)
	  d <- c(rep(0,q),x$d[-c(1:q)])/x$d
	  return( U%*%diag(d/n)%*%t(U)%*%X )
	 }
	}
	
	if (is.null(group)){
	 Z <- estimate_FX(Y, q=hcount)
	}else{
	 Z_0 <- estimate_FX(Y_0, q=hcount)
	 Z_1 <- estimate_FX(Y_1, q=hcount)
	 group<- c(rep(0,nrow(Z_0)),rep(1,nrow(Z_1)))
	 Z <- cbind(group, rbind(Z_0, Z_1))
	}

	return(data = Z)
}

estimate_latent_count <- function(X, method, verbose, ...)
{
	n <- nrow(X)
	p <- ncol(X)
	r <- min(n-1,p)
	svdX <- svd(scale(X), nu=0, nv=r)
	seed <- 1
		
	if (method == "auto") {
	# Parallel analys by N permutation (Dobriban, 2020)
	 permutation_thresholding <- function(X, quantile = 0.99, seed = 1)
	 {
		N <- 50
		X <- scale(X)
		X <- X[, !is.na(apply(X, 2, sum))]
		X <- X[, sort(apply(X, 2, sd),
					index.return = TRUE,
					decreasing = TRUE)$ix[seq_len(min(1000, ncol(X)))]]
		r <- min(dim(X))
		
		evals <- matrix(0, nrow = N, ncol = r)
		set.seed(seed)
		for (i in seq_len(N)) {
		 X_perm <- apply(X, 2, function(xx) sample(xx))
		 evals[i, ] <- svd(X_perm, nu = 0, nv = 0)$d[seq_len(r)]
		}
		thresholds <- apply(
		 evals, 2,
		 function(xx) quantile(xx, probs = quantile)
		)

		# limit number of confounders to at most 10 LVs
		limit <- ifelse(r > 10, 10, r)
		# last which crosses the threshold
		crosses <- (svd(X, nu = 0, nv = 0)$d > thresholds)[1:limit]
		return(max(which(c(TRUE, crosses)) - 1))
	 }
	 q <- permutation_thresholding(X, seed = seed)
	 idx <- ceiling(q)
	}
		
	if (method == "knee") {
	# Screening at the knee point (Raiche et al, 2013)
	 ev <- svdX$d / sqrt(n)
	 ev <- ev[seq_len(round(length(ev) / 2))]
	 values <- seq(length(ev))

	 d1 <- diff(ev) / diff(values) # first derivative
	 d2 <- diff(d1) / diff(values[-1]) # second derivative
	 idx <- which.max(abs(d2))
	}
	
	if (method == "ED") {
	# Screening at the Empirical Distribution (ED) (Onatsky, 2010) 
	 ev <- svdX$d^2 / n
	 rmax <- round(3 * sqrt(r))
	 j <- rmax + 1
     diffs <- ev - c(ev[-1], 0)

	 for (i in 1:10) { #i=1
		y <- ev[j:(j+4)]
		x <- ((j-1):(j+3))^(2/3)
		lm.coef <- lm(y ~ x)
		delta <- 2 * abs(lm.coef$coef[2])
		idx <- which(diffs[1:rmax] > delta)
		hatr <- ifelse(length(idx) == 0, 0, max(idx))
		newj <- hatr + 1
		if (newj == j) break
		j <- newj
	}
	 idx <- hatr
    }
	
	if (method == "IC") {
	# Screening at min(Information Criteria)(Bai and Ng, 2003)
	 V <- svdX$v
	 M <- ifelse(r > 10, 10, r)
	 K_losses <- c()
	 for (K in 1:M){ #K=1
	  Vk <- V[,1:K]
	  #css <- t(Vk %*% (t(Vk) %*% t(X)))
	  css <- X%*%Vk%*%t(Vk)
	  loss <- log(sum(X - css)^2/(n*p)) # log likelihood
	  loss <- loss + K*((p+n)/(p*n))*log((p*n)/(p+n)) # penalizer
	  K_losses <- rbind(K_losses, loss)
	 }
	 idx <- which(K_losses == min(K_losses))
	}
	
	# visualize the scree plot
	ev <- svdX$d^2/n
	if (verbose) {
	 plot(c(1:r), ev[1:r], type = "p", ylab = "Eigenvalue",
	  xlab= "Number of principal components", cex = 1.1,
	  cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
	 abline(h = ev[idx] - 0.005, lty = 2, col = "red")
	}
	pve <- round(cumsum(ev[1:idx])[idx]/sum(ev),2)
	#logger::log_info("Estimated {idx} latent confounders")
	message(paste0("Estimated ", idx," (", 100*pve,"%) latent confounders"))

	return(idx)
}

map_hidden_dag <- function(graph, data, cg=NULL, verbose=FALSE, ...)
{
	VH <- colnames(data)[grepl("LV",colnames(data))]
	gH <- graph + igraph::vertices(VH)
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

#' @title Estimate a DAG from an input (or empty) graph
#'
#' @description Two-step extraction of the optimal DAG from an input (or empty)
#' graph, using in step 1) graph topological order or bottom-up search order,
#' and in step 2) parent recovery with the LASSO-based algorithm (FHT, 2010),
#' implemented in \code{\link[glmnet]{glmnet}}.
#'
#' @param graph An igraph object or a graph with no edges (make_empty_graph(n=0)).
#' @param data A matrix whith n rows corresponding to subjects, and p columns
#' to graph nodes (variables).
#' @param LO character for linear order method. If LO="TO" or LO="TL" the
#' topological order (resp. level) of the input graph is enabled, while LO="BU"
#' the data-driven bottom-up search of vertex (resp. layer) order is performed
#' using the vertices of the empty graph. By default \code{LO = "TO"}.
#' @param beta Numeric value. Minimum absolute LASSO beta coefficient for
#' a new direct link to be retained in the final model. By default,
#' \code{beta = 0}.
#' @param eta Numeric value. Minimum fixed eta threshold for bottom-up search
#' of vertex (eta = 0) or layer (eta > 0) ordering. Use eta = NULL, for estimation
#' of eta adaptively with half of the sample data. By default, \code{eta = 0}.
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
#' @details The extracted DAG is estimated using the two-step order search approach.
#' First a vertex (node) or level (layer) order of p nodes is determined, and from
#' this sort, the DAG can be learned using in step 2) penalized (L1) regressions
#' (Shojaie and Michailidis, 2010). The estimate linear order are obtained from
#' \emph{a priori} graph topological vertex (TO) or level (TL) ordering, or with a
#' data-driven Bottom-up (BU) approach, assuming a SEM whose error terms have equal
#' variances (Peters and Bühlmann, 2014). The BU algorithm first estimates the last
#' element (the terminal vertex) using the diagonal entries of the inverse covariance
#' matrix with: t = argmin(diag(Omega)), or the terminal layer (> 1 vertices) with
#' d = diag(Omega)- t < eta. And then, it determines its parents with L1 regression.
#' After eliminating the last element (or layer) of the ordering, the algorithm applies
#' the same procedure until a DAG is completely estimated. In high-dimensional data 
#' (n < p), the inverse covariance matrix is computed by glasso-based algorithm
#' (FHT, 2008), implemented in \code{\link[glasso]{glasso}}. If the input graph is
#' not acyclic, in TO or TL, a warning message will be raised, and a cycle-breaking
#' algorithm will be applied (see \code{\link[SEMgraph]{graph2dag}} for details).
#' Output DAG will be colored: vertices in cyan, if they are source nodes, and in
#' orange, if they are sink nodes, and edges in gray, if they were present in the
#' input graph, and in green, if they are new edges generated by LASSO screening.
#'
#' @return A list of 3 igraph objects plus the vertex ordering:
#' \enumerate{
#' \item "dag", the estimated DAG;
#' \item "dag.new", new estimated connections;
#' \item "dag.old", connections preserved from the input graph;
#' \item "LO", the estimated vertex ordering.
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
#' Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse covariance
#' estimation with the graphical lasso. Biostatistics, 9(3), 432-441.
#' <https://doi.org/10.1093/biostatistics/kxm045>
#'
#' Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for
#' Generalized Linear Models via Coordinate Descent.
#' Journal of Statistical Software, Vol. 33(1), 1-22.
#' <https://doi.org/10.18637/jss.v033.i01>
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
#' <https://doi.org/10.1093/biomet/ast043>
#'
#' @examples
#'
#' #Set function param
#' ig <- sachs$graph
#' X <- log(sachs$pkc)
#' group <- sachs$group
#'
#' # DAG estimation (default values)
#' dag0 <- SEMdag(ig, X)
#' sem0 <- SEMrun(ig, X, group)
#'
#' # Graphs
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(1,4))
#' plot(sachs$graph, layout=layout.circle, main="input graph")
#' plot(dag0$dag, layout=layout.circle, main = "Output DAG")
#' plot(dag0$dag.old, layout=layout.circle, main = "Inferred old edges")
#' plot(dag0$dag.new, layout=layout.circle, main = "Inferred new edges")
#' par(old.par)
#'
#' # Four DAG estimation
#' dag1 <- SEMdag(ig, X, LO="TO")
#' dag2 <- SEMdag(ig, X, LO="TL")
#' dag3 <- SEMdag(ig, X, LO="BU", eta=0)
#' dag4 <- SEMdag(ig, X, LO="BU", eta=NULL)
#' 
#' unlist(dag1$LO)
#' dag2$LO
#' unlist(dag3$LO)
#' dag4$LO
#' 
#' # Graphs
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(2,4))
#' gplot(dag1$dag, main="TO")
#' gplot(dag2$dag, main="TL")
#' gplot(dag3$dag, main="BU")
#' gplot(dag4$dag, main="TLBU")
#' par(old.par)
#' 
SEMdag<- function(graph, data, LO="TO", beta=0, eta=NULL, lambdas=NA, penalty=TRUE, verbose=FALSE, ...)
{
	# Set DAG objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	if (!is_dag(ig) & (LO == "TO" | LO == "TL")){
	 cat("WARNING: input graph is not acyclic !\n")
	 cat(" Applying graph -> DAG conversion...\n")
	 dag<- graph2dag(ig, data) #del cycles & all <->
	}else{ dag<- ig }
	#X<- as.matrix(scale(data[, V(dag)$name]))
	X<- scale(data[, V(dag)$name])

	# Estimate DAG using linear ordering (LO) approach:
	x<- getParents(dag, X, LO, beta, eta, lambdas, penalty, verbose)
	#if (LO == "TO") x$L<- unlist(x$L)
	
	# Mapping DAG edges on input graph:
	ig1<- x$gxy
	ig2<- quiet(properties(ig1)[[1]])
	E1<- attr(E(ig2), "vnames")
	E0<- attr(E(ig), "vnames")
	E(ig2)$color<- ifelse(E1 %in% E0, "gray", "green")
	dout <- igraph::degree(ig2, mode= "out") #sink
	din <- igraph::degree(ig2, mode= "in") #source
	V(ig2)$color[dout == 0]<- "orange"
	V(ig2)$color[din == 0]<- "cyan"
	ig3<- ig2-E(ig2)[which(E(ig2)$color == "gray")]
	ig3<- ig3-vertices(V(ig3)$name[igraph::degree(ig3) == 0])
	ig4<- ig2-E(ig2)[which(E(ig2)$color == "green")]
	ig4<- ig4-vertices(V(ig4)$name[igraph::degree(ig4) == 0])
	
	if (verbose) {
	 gplot(ig2)
	 fit<- SEMrun(ig2, X, algo = "ricf", n_rep = 0)
	}

	return( list(dag=ig2, dag.new=ig3, dag.old=ig4, LO=x$L) )
}

getParents<- function(dag, X, LO, beta, eta, lambdas, penalty, verbose, ...)
{
	# STEP 1: Variables ordering recovery
	cat("Node Linear Ordering with", LO, "setting\n\n")
	if (LO == "TO"){
	 TO <- names(igraph::topo_sort(dag))
	 L <- as.list(rev(TO))
	}else if (LO == "TL"){
	 L <- buildLevels(dag)
	}else if (LO == "BU"){
	 L <- buildLayers(X, rho = NULL, eta = eta)
	}
	X <- as.matrix(X)
	n <- nrow(X)
	p <- length(L)
	if (is.na(lambdas)) lambdas <- sqrt(log(p)/n)
	adj <- as_adj(dag, sparse = FALSE)

	# STEP 2: DAG recovery with BU ordering
	ftm<- NULL
	for (i in 1:(p - 1)) { #i=1
	  yy<- which(colnames(X) %in% L[[i]])
	  ll<- unique(unlist(L[(i+1):length(L)]))
	  xx<- which(colnames(X) %in% ll)
	  y<- as.matrix(X[,yy])
	  x<- as.matrix(X[,xx]) 
	  if (ncol(y) == 1) colnames(y)<- colnames(X)[yy]
	  if (ncol(x) == 1) colnames(x)<- colnames(X)[xx]
	  if (penalty) {
		 pw <- matrix(1,ncol(x),ncol(y)) - adj[xx, yy]
	  }else{
		 pw <- matrix(1,ncol(x),ncol(y))
	  }
	  ftmi<- glmnet.set(x, y, beta, lambdas, pw)
	  ftm<- rbind(ftm, ftmi)
	}

	if (is.null(ftm))
 	 return(message("DAG with 0 edges, change input values !"))
	del<- which(duplicated(ftm) == TRUE)
	if (length(del) > 0) ftm<- ftm[-del,]
	gxy<- graph_from_data_frame(ftm, directed=TRUE)
		
	return(list(gxy = gxy, L = L))
}

glmnet.set<- function(x, y, beta, lambdas, pw, ...)
{ 
	set.seed(1324)
	ftm <- NULL
	for (k in 1:ncol(y)) {
	   if (ncol(x) > 1) {
		if (!is.null(lambdas)) {
			fit <- glmnet::glmnet(x=x, y=y[,k], family="gaussian",
					#lambda = sqrt(log(ncol(y)+ncol(x))/nrow(y)))
					lambda = lambdas,
					penalty.factor = pw[,k])
			b <- coef(fit, s = NULL)[-1,]
		}else{
			fit <- glmnet::cv.glmnet(x=x, y=y, family="gaussian",
					lambda = NULL,
					penalty.factor = pw[,k])
			b <- coef(fit, s = fit$lambda.min)[-1,]
		}
	   }else{
			fit <- summary(lm(y[,k] ~ x))$coefficients[2,4]
			b <- ifelse(fit < 0.05, 1, 0)
			names(b) <- colnames(x)
	   }
	   from <- names(b)[abs(b) > beta]
	   to <- rep(colnames(y)[k], length(from))
	   ftm <- rbind(ftm, cbind(from, to))
	}

	return(ftm)
}

buildLevels <- function(dag, ...)
{
	Leaf_removal <- function(dag)
	{
	 levels <- list()
	 level <- 1
	 repeat {
	  leaves <- igraph::degree(dag, mode= "out")
	  levels[[level]]<- names(leaves)[leaves == 0]
  	  dag <- delete_vertices(dag, names(leaves)[leaves == 0])
	  level <- level+1
	  if (vcount(dag)==0 | ecount(dag)==0) break
	 }
	 levels[[level]] <- V(dag)$name
	 names(levels)<- 1:level
	 return(levels)
	}

	# leaf-removal(dag)
	l1<- Leaf_removal(dag)
	if (length(l1) == 2) return(l1)
	# leaf removal(dagT)
	adj <- as_adj(dag, sparse=FALSE)
	dagT <- graph_from_adjacency_matrix(t(adj), mode="directed")
	l2 <- Leaf_removal(dagT)
	l2 <- rev(l2)
	# number-of-layers 
	L <- max(length(l1), length(l2))

	# combine BU-ordering (dag+dagT)
	l3 <- list()
	l3[[1]] <- l1[[1]] #sink
	l3[[L]] <- l2[[L]] #source
	for (k in 2:(L-1)){
	 lk <- unique(c(l1[[k]], l2[[k]]))
	 Lk <- unlist(l3[c(1:(k-1),L)])
	 l3[[k]] <- setdiff(lk, Lk)
	}

	return(l3)
}

buildLayers <- function(X, rho = NULL, eta = NULL, eta.scaler = 1)
{
	### Set the required parameters from data matrix
	X <- as.matrix(X)
	p <- ncol(X)
	n <- nrow(X)
	#if (is.null(rho)) rho <- sqrt(log(p)/n)

	### Layer Estimation
	layers <- list()
	l <- 1
	RemNode <- colnames(X) ## initial remaining nodes

	while (length(RemNode) > 0) {
	 ### Estimate residual variance for each RemNode
	 X1 <- X[, RemNode]
	 if (is.null(rho)) rho <- sqrt(log(ncol(X1))/n)
	 if (ncol(as.matrix(X1)) > 1) {
	  if (n > p) {
	   wi <- solve(cov(X1))
	  }else{
	   wi <- glasso::glasso(s=cov(X1), rho=rho)$wi
	   rownames(wi) <- colnames(wi) <- colnames(X1)
	  }
	  resvars <- diag(wi)
	  diff.resvars <- resvars - min(resvars)
	  #Recover layers using X2 to determine eta adaptively 
	  if (is.null(eta)) {
	   X2 <- X[sample(1:n, n/2), RemNode]
	   wi <- glasso::glasso(s=cov(X2), rho=sqrt(2)*rho)$wi
	   rownames(wi) <- colnames(wi) <- colnames(X2)
	   resvars2 <- diag(wi)
	   eta_hat <- min(abs(resvars - resvars2)) * eta.scaler
	   layers[[l]] <- RemNode[diff.resvars <= eta_hat]
	  }else{
	  #Recover layers with fixed eta    
	   layers[[l]] <- RemNode[diff.resvars <= eta]
	  }
	 }else{
	  layers[[l]] <- RemNode
	 }

	 ### Updating the remaining node ###
	 RemNode<- RemNode[-which(RemNode %in% layers[[l]])]

	 l <- l + 1 # layer indicator
	}

	return(layers)
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
#' and the second graph is the output of either \code{\link[SEMgraph]{SEMdag}}
#' or \code{\link[SEMgraph]{SEMbap}}. Alternatively, the first graph is an
#' empthy graph, and the second graph is a external covariance graph.
#' In the former we use the new inferred causal structure stored in the
#' \code{dag.new} object. In the latter, we use the new inferred covariance
#' structure stored in the \code{guu} object. Both directed (causal) edges
#' inferred by \code{SEMdag()} and covariances (i.e., bidirected edges)
#' added by \code{SEMbap()}, highlight emergent hidden topological
#' proprieties, absent in the input graph. Estimated directed edges between
#' nodes X and Y are interpreted as either direct links or direct paths
#' mediated by hidden connector nodes. Covariances between any two bow-free
#' nodes X and Y may hide causal relationships, not explicitly represented
#' in the current model. Conversely, directed edges could be redundant or
#' artifact, specific to the observed data and could be deleted.
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
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
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
#' with graph and data-driven algorithms.
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param seed A vector of seed nodes.  
#' @param type Tree-based structure learning method. Four algorithms 
#' are available:
#' \itemize{
#' \item "ST"(default). Steiner Tree (ST) identification via fast Kou's algorithm 
#' (Kou et al, 1981) connecting a set of seed nodes (called Terminal vertices)
#' with connector nodes (called Steiner vertices) from input graph as defined
#' in \code{graph} with minimal total distance on its edges. By default the edge
#' weights are based on the pairwise correlation, 1-abs(cor(j,k)). If input
#' graph has E(graph)$weight=1, and \code{eweight = "custom"}, ST seeks a minimum
#' subtree (i.e., the subtree with minimal number of edges).
#' \item "CAT". Causal additive trees (CAT) algorithm as in Jakobsen et al. 
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
#' no seed nodes are needed).
#' }
#' @param eweight Edge weight type for igraph object can be externally derived
#' using \code{\link[SEMgraph]{weightGraph}} or from user-defined distances. 
#' This option determines the weight-to-distance transform. If set to:
#' \itemize{
#' \item "NULL" (default), edge weights will be internally computed
#' equal to 1 - abs(pairwise Pearson's correlation).
#' \item "kegg", repressing(-1), neutral(0) and activating(+1) kegg
#' interactions will be multiplied by "zsign" attributes, and positive
#' (i.e., concordant) values will be set to 1 (minimum distance), while
#' negative (i.e., discordant) values will be set to 2.
#' \item "zsign", all significant interactions (abs(zsign) > 0) will be
#' set to 1 (minimum distance), while non-significant (zsign=0) ones will
#' be set to 2.
#' \item "pvalue", edge p-value atributes will be transformed to the
#' inverse of negative base-10 logarithm, 1/(-log(E(graph)$pv)).
#' \item "custom", the algorithm will use the distance measure specified
#' by the user as "weight" edge attribute in the input graph.
#' }
#' @param alpha Threshold for rejecting a pair of node being independent in 
#' "CPDAG" algorithm. The latter implements a natural v-structure identification 
#' procedure by thresholding the pairwise sample correlations over all adjacent 
#' pairs of edges with some appropriate threshold. By default, 
#' \code{alpha = 0.05}.
#' @param verbose If TRUE, it shows the output tree (not recommended for large graphs).
#' @param ... Currently ignored.
#'
#' @details A tree ia an acyclic graph with p vertices and p-1 edges. The graph method
#' refers to the  Steiner Tree (ST), a tree from an undirected graph that connect "seed"
#' with additional nodes in the "most compact" way possible. The data-driven methods
#' propose fast and scalable procedures based on Chu-Liu–Edmonds’ algorithm (CLE) to
#' recover a tree from a full graph. The first method, called Causal Additive Trees (CAT)
#' uses pairwise mutual weights as input for CLE algorithm to recover a directed tree
#' (an "arborescence"). The second one applies CLE algorithm for skeleton recovery and
#' extends the skeleton to a tree (a "polytree") represented by a Completed Partially
#' Directed Acyclic Graph (CPDAG). Finally, the Minimum Spanning Tree (MST) connecting
#' an undirected graph with minimal edge weights can be identified.
#' To note, if the input graph is a directed graph, ST and MST undirected trees are
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
#' Grassi M, Tarantino B (2023). SEMtree: tree-based structure learning methods
#' with structural equation models. 
#' Bioinformatics, 39 (6), 4829–4830 <https://doi.org/10.1093/bioinformatics/btad377>
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
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
#'
#' # graph-based trees
#' graph <- alsData$graph
#' seed <- V(graph)$name[sample(1:vcount(graph), 10)]
#' tree1 <- SEMtree(graph, als.npn, seed=seed, type="ST", verbose=TRUE)
#' tree2 <- SEMtree(graph, als.npn, seed=NULL, type="MST", verbose=TRUE)
#'
#' # data-driven trees
#' V <- colnames(als.npn)[colnames(als.npn) %in% V(graph)$name]
#' tree3 <- SEMtree(NULL, als.npn, seed=V, type="CAT", verbose=TRUE)
#' tree4 <- SEMtree(NULL, als.npn, seed=V, type="CPDAG", alpha=0.05, verbose=TRUE)
#'
#' }
#'
SEMtree <- function(graph, data, seed, type = "ST", eweight = NULL, alpha = 0.05, verbose = FALSE, ...)
{
	# Set data and graph objects:
	if (!is.null(graph)) {
	 nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	 ig <- induced_subgraph(graph, vids = V(graph)$name %in% nodes)
	 X <- data[,nodes]
	 if (!is.null(eweight)) {
	  if (eweight == "kegg") E(graph)$weight <- 2 - E(graph)$weight*E(graph)$zsign
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

	# Causal Addittive Tree(CAT) or CPDAG Tree:
	}else if(is.null(graph)) {
	  nodes <- colnames(data)[colnames(data) %in% seed]
	  X <- data[, nodes]
	  if (type == "CAT") T <- CAT.R(data = data.frame(X))
	  if (type == "CPDAG") T <- CPDAG(X, alpha, verbose=TRUE)
	  #if (type == "MST") {
	  # A <- 1-abs(cor(X))
	  # gA <- graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
	  # T <- mst(gA, weights = E(gA)$weight, algorithm = "prim")
	  #}
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
	colNames <- sub(".", "", colnames(data))
	colnames(data) <- paste0("X",seq(1,ncol(data),1))

	# Compute Gaussian Edge Weights:
	if (is.data.frame(data)) {
	 ig0 <- make_full_graph(ncol(data), directed = TRUE)
	 Edges <- igraph::as_data_frame(ig0)
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
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
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
						verbose = FALSE, ...)
{
	# Step by step search:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	Zt<- as.matrix(data[,nodes])
	Gt<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	cat("Step1: BAP deconfounding...\n")
	Zt1<- quiet(SEMbap(Gt, Zt, method=method, alpha=alpha))
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
	 Gt2<- quiet(resizeGraph(g=list(Gt,Gt1.new), gnet, d=1, v=FALSE, verbose=FALSE))
	 dataZ<- Zt1$data
	}
	if (search == "inner") {
	 Gt2<- quiet(resizeGraph(g=list(Gt,Gt1.new), gnet, d=d, v=FALSE, verbose=FALSE))
	 dataZ<- Zt1$data
	}
	if (search == "outer") {
	 Gt2<- quiet(resizeGraph(g=list(Gt,Gt1.new), gnet, d=d, verbose=FALSE))
	 green<- V(Gt2)$name[which(V(Gt2)$color == "green")]
	 dataZ<- cbind(Zt1$data, data[,which(colnames(data) %in% green)])
	}
	cat("Done.\n")
	if (verbose) {
	 gplot(Gt2, main="Estimated Extended Graph")
	 cat("\n")
	 fit<- SEMrun(Gt2, dataZ, algo = "ricf", n_rep = 0)
	 C_test<- Shipley.test(Gt2, dataZ, verbose=TRUE)
	}
	
	return(list(graph = Gt2, data = dataZ))
}
