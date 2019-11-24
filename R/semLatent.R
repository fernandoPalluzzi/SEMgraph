# SEMgraph 0.3.3
#
# This is the SEMgraph package for statistical graph analysis using
# Structural Equation Models
# Typically, SEMgraph requires three input elements:
# 1) Interactome (e.g. KEGG singaling pathways or STRING PPI)
# 2) Quantitative data (e.g. GWAS, DNA-methylation arrays, RNA-seq)
# 3) Case/Control identifiers
#

ggm.set <- function(x, y, alpha) {
	# P-values via Large-Scale (glasso) Gaussian Graphical Model
	Z <- cbind(x, y)
	#Z[1:6,1:6]
	t <- ncol(x)*ncol(y)
	wi <- SILGGM::SILGGM(Z, method = "D-S_GL")$p_precision
	rownames(wi) <- colnames(wi) <- colnames(Z)
	#wi[1:6,1:6]
	rxy <- 1*as.matrix(wi[(ncol(x) + 1):ncol(Z), 1:ncol(x)] < alpha/t)
	return(rxy)
}

psi2guv <- function(guu, ig, gnet, alpha, verbose, ...)
{
	# Adding new nodes to undirected graph (guu) from the interactome  
	vids <- which(V(guu)$name %in% V(gnet)$name)
	guu <- induced_subgraph(graph = guu, vids = vids)
	if (verbose) {
		plot(guu, main = "Covariance network (guu)")
		Sys.sleep(5)
	}
	ftm <- as_edgelist(guu)
	vpath <- ftmuv <- NULL
	
	for (i in 1:nrow(ftm)) { #i=1
		if(distances(gnet, ftm[i, 1], ftm[i, 2], weights = NA) == Inf) next
		path <- shortest_paths(gnet, ftm[i, 1], ftm[i, 2],
		                       mode = "all",
		                       weights = NA)$vpath[[1]]
		V <- V(gnet)$name[path]
		#V; V(guu)$name
		vpath <- c(vpath, V[-c(1, length(V))])
		#vpath
		for(h in 1:(length(V) - 1)) ftmuv <- rbind(ftmuv, c(V[h], V[h + 1]))
	}
	
	# Graph with added nodes from interactome (guv)
	ftmuv <- stats::na.omit(ftmuv[duplicated(ftmuv) != TRUE,])
	if (nrow(ftmuv) > 0) {
		guv <- graph_from_data_frame(ftmuv, directed = FALSE)
		guv <- simplify(guv, remove.loops = TRUE)
		vv <- V(guv)$name[-which(V(guv)$name %in% V(ig)$name)]
		uv <- V(ig)$name[which(V(ig)$name %in% unique(vpath))]
		
		V(guv)$color[V(guv)$name %in% V(guu)$name] <- "lightblue"
		V(guv)$color[V(guv)$name %in% vv] <- "yellow"
		V(guv)$color[V(guv)$name %in% uv] <- "green"
		#E(guv)$weight <- 1
		if(verbose) {
			plot(guv, main = "Extended connector network (guv)")
			Sys.sleep(5)
		}
	} else {
		guv <- NA
		cat("\nno edges u--v found !\n\n")
	}
	
	return(guv)
}

psi2gUX <- function(duv, ig, gnet, data, alpha, verbose, ...)
{
	# Select the source nodes from the directed graph guv
	degree.in <- igraph::degree(ig, v = V(ig)$name, mode = "in")
	xx <- V(ig)$name[degree.in == 0]
	cxx <- V(duv)$name[V(duv)$name %in% xx]
	
	if (length(cxx) > 1) {
		AN <- unique(unlist(neighborhood(gnet, nodes = cxx, mode = "in")))
		#AN <- parents(gnet, nodes = cxx)
		#AN <- ancestors(gnet, nodes = cxx)
		AN <- V(gnet)$name[AN]
		#length(AN)
		if (length(AN) == 0) {
			cat("\nno AN(x) found !\n\n")
			return(NA)
		}
		Y <- as.matrix(data[, cxx])
		#head(Y)
		idx <- colnames(data) %in% AN[-which(AN %in% cxx)]
		X <- as.matrix(data[, idx])
		#head(X)
		
		Bxx <- ggm.set(x = as.matrix(X), y = as.matrix(Y), alpha = alpha)
		if (sum(Bxx) == 0) {
			cat("\nno AN(x) found !\n\n")
			return(NA)
		}
		ftmX <- NULL
		for (i in 1:nrow(Bxx)) { #i=1
			from <- colnames(Bxx)[Bxx[i,] == 1]
			to <- rownames(Bxx)[i]
			ftmX <- rbind(ftmX, cbind(from, rep(to, length(from))))
		}
		
		gUX <- graph_from_data_frame(ftmX, directed = TRUE)
		V(gUX)$color <- "yellow"
		V(gUX)$color[V(gUX)$name %in% cxx] <- "orange"
		if (verbose) {
			plot(gUX, main = "Extended sources network (gUX)")
			Sys.sleep(5)
		}
	} else {
		cat("\nno psi(x,x) found !\n\n")
		return(NA)
	}
	
	return(gUX)
}

psi2gUY <- function(duv, ig, gnet, data, alpha, verbose, ...)
{
	# Select the sink nodes from the directed graph guv
	degree.out <- igraph::degree(duv, v = V(duv)$name, mode = "out")
	yy <- V(ig)$name[degree.out == 0]
	cyy <- V(duv)$name[V(duv)$name %in% yy]
	
	if (length(cyy) > 1) {
		DE <- unique(unlist(neighborhood(gnet, nodes = cyy, mode = "out")))
		#DE <- siblings(gnet, nodes = cyy)
		#DE <- descendants(gnet, nodes = cyy)
		DE <- V(gnet)$name[DE]
		if (length(DE) == 0) {
			cat("\nno DE(y) found !\n\n")
			return(NA)
		}
		Y <- as.matrix(data[, cyy])
		#head(Y)
		idx <- colnames(data) %in% DE[-which(DE %in% cyy)]
		X <- as.matrix(data[, idx])
		#head(X)
		Byy <- ggm.set(x = as.matrix(X), y = as.matrix(Y), alpha = alpha)
		if (sum(Byy) == 0) {
			cat("\nno DE(y) found !\n\n")
			return(NA)
		}
		ftmY <- NULL
		for (i in 1:nrow(Byy)) { #i=1
			from <- rownames(Byy)[i]
			to <- colnames(Byy)[Byy[i,] == 1]
			ftmY <- rbind(ftmY, cbind(rep(from, length(to)), to))
		}
		
		gUY <- graph_from_data_frame(ftmY, directed = TRUE)
		V(gUY)$color <- "yellow"
		V(gUY)$color[V(gUY)$name %in% cyy] <- "orange"
		if (verbose) {
			plot(gUY, main = "Extended target network (gUY)")
			Sys.sleep(5)
		}
	} else {
		cat("\nno psi(y,y) found !\n\n")
		return(NA)
	}
	
	return(gUY)
}

#' @title Network community plotting utility
#'
#' @description Merge and plot network communities of a graph as single 
#' nodes.
#' @param g An igraph object.
#' @param membership a vector of node memberships.
#' @param l igraph layout option.
#' @param global Logical value. If TRUE, the plot of the input graph 
#' (coloured by cluster membership) will be generated alongwith independent 
#' module plots. If the input graph is very large, plotting could be 
#' computationally intensive (by default, global = FALSE).
#'
#' @import igraph
#' @importFrom graph nodes edgeNames isDirected nodeRenderInfo edgeRenderInfo graphRenderInfo
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @export
#'
#' @examples
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' membership <- clusterGraph(graph, type = "tahc", size = 10)
#' G <- cplot(graph, membership, global = TRUE)
#'
cplot <- function(g, membership, l = layout.auto, global = FALSE) {
	
	# Cluster visualization
	V(g)$M[which(V(g)$name %in% names(membership))] <- membership
	if(is.character(V(g)$M)) V(g)$M <- substr(V(g)$M, 3, 10L)
	V(g)$M[is.na(V(g)$M)] <- 999
	V(g)$color <- as.numeric(V(g)$M) + 1
	try(weight <- E(g)$weight)
	E(g)$weight <- 1
	if(global) {
		gplot(g)
		Sys.sleep(5)
	}
	
	# Within cluster visualization
	M <- names(table(V(g)$M))
	K <- length(table(V(g)$M))
	vcol <- as.numeric(M) + 1
	
	HM <- lapply(1:K, function(x) induced_subgraph(g, V(g)$name[V(g)$M == M[x]]))
	names(HM) <- paste0("HM", M)
	d <- igraph::degree(g, mode = "all")*2 + 1
	glv <- lapply(1:K, function(x) {
		plot(HM[[x]], vertex.color = vcol[x],
		     vertex.size = d[V(HM[[x]])$name],
		     main = paste0("Hidden Module ", M[x]))
		Sys.sleep(3)
	})
	try(E(g)$weight <- weight)
	
	return(invisible(c(list(g = g), HM)))
}

#' @title Hidden modules building
#'
#' @description Graph clustering utility.
#' @param graph An igraph object.
#' @param HM Hidden model type. This parameter is only required by 
#' \code{\link[SEMgraph]{SEMfsr}}), and it is automatically set to NULL 
#' when clusterGraph is used as a stand alone function.
#' For each defined hidden module: 
#' (i) if HM = "LV", a latent variable (LV) will be defined as common 
#' unknown cause acting on cluster nodes; (ii) if HM = "CV", cluster nodes 
#' will be considered as regressors of a latent composite variable (CV); 
#' (iii) if HM = "UV", an unmeasured variable (UV) is defined, where source 
#' nodes of the module (i.e., in-degree = 0) act as common regressors 
#' influencing the other nodes via an unmeasured variable (see also 
#' \code{\link[SEMgraph]{SEMfsr}}).
#' @param type Network clustering method. If type = "ebc" (default), network 
#' modules are generated using edge betweenness clustering method (see 
#' \code{\link[igraph]{cluster_edge_betweenness}}). Other clustering methods 
#' include: "fgc" (fast greedy method), "lbc" (label propagation method), 
#' "lec" (leading eigenvector method), "loc" (multi-level optimization), 
#' "opc" (optimal communiy structure), "sgc" (spinglass statistical 
#' mechanics), "wtc" (walktrap community structure with short random walks), 
#' and "tahc" (tree agglomerative hierarchical clustering).
#' @param size Minimum number of nodes per hidden module. By default, a 
#' minimum number of 3 nodes is required.
#' @param verbose A logical value. If FALSE (default), the processed graphs 
#' will not be plotted to screen, saving execution time (they will be 
#' returned anyway).
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @importFrom stats qnorm cov2cor cor hclust as.dist cutree
#' @importFrom graphics abline
#' @export
#'
#' @references
#'
#' Fortunato S, Hric D. Community detection in networks: A user guide (2016). 
#' Phys Rep; 659: 1-44. http://dx.doi.org/10.1016/j.physrep.2016.09.002
#' 
#' Yu M, Hillebrand A, Tewarie P, Meier J, van Dijk B, Van Mieghem P, 
#' Stam CJ (2015). Hierarchical clustering in minimum spanning trees. 
#' Chaos 25(2): 023107.  https://doi.org/10.1063/1.4908014
#'
#' @return If HM is not "none" a list of 3 objects is returned:
#' \enumerate{
#' \item "gHM", subgraphs containing latent variable or composite variable modules;
#' \item "M", cluster membership vector for each node;
#' \item "gHC", the list of modules produced by clusterGraph.
#' }
#' If HM is "none", only the cluster membership vector is returned.
#'
#' @seealso \code{\link[SEMgraph]{SEMfsr}}
#' @examples
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' membership <- clusterGraph(graph, type = "tahc", size = 10, verbose = TRUE)
#'
clusterGraph <- function(graph, type, HM = "none", size = 3, verbose = FALSE, ...)
{
	# Clustering by cluster_edge_betweenness() or by TAHC()
	if (!is_directed(graph)) {
		ug <- graph
	} else {
		ug <- as.undirected(graph,
		                    mode = "collapse",
		                    edge.attr.comb = "ignore")
	}
	
	if (type == "tahc") {
		# Tree Agglomerative Hierarchical Clustering (TAHC)
		mst <- minimum.spanning.tree(ug,
		                             weights = NULL,
		                             algorithm = NULL)
		G <- distances(mst, v = V(mst), to = V(mst),
		               mode = "all",
		               weights = NA)
		D <- 1 - cor(x = G, method = "spearman")
		#D[1:10, 1:10]
		hMST <- hclust(as.dist(D), method = "average")
		tahc <- cutree(hMST, h = 0.2)
		#print(table(tahc)); cat("\n")
		cnames <- as.numeric(names(table(tahc)))[table(tahc) >= size]
		membership <- tahc[tahc %in% cnames]
		if(verbose) {
			plot(hMST, labels = FALSE, xlab = "", sub = "")
			abline(h = 0.2, col = "red")
			Sys.sleep(5)
		}
	} else {
		#type= "wtc"
		if (type == "ebc") cls <- cluster_edge_betweenness(ug, weights = NULL)
		if (type == "fgc") cls <- cluster_fast_greedy(ug, weights = NULL)
		if (type == "lbc") cls <- cluster_label_prop(ug, weights = NA)
		if (type == "lec") cls <- cluster_leading_eigen(ug, weights = NA)
		if (type == "loc") cls <- cluster_louvain(ug, weights = NA)
		#if (type == "opc") cls <- cluster_optimal(ug, weights = NA)
		if (type == "sgc") cls <- cluster_spinglass(ug, weights = NA)
		if (type == "wtc") cls <- cluster_walktrap(ug, weights = NULL)
		cat("modularity =", modularity(cls), "\n\n")
		print(sort(sizes(cls)))
		cat("\n")
		cnames <- as.numeric(names(sizes(cls)[sizes(cls) >= size]))
		membership <- membership(cls)[membership(cls) %in% cnames]
		if(verbose) {
			plot(cls, ug)
			Sys.sleep(5)}
	}
		
	K <- length(cnames)
	if (K == 0) return(cat("No communities with size >", size, "found!\n"))
	
	if (HM == "UV") {
		V(graph)$M <- 999
		V(graph)$M[which(V(graph)$name %in% names(membership))]<- membership
		M <- names(table(V(graph)$M))
		K <- length(table(V(graph)$M))
		gHC <- lapply(1:K, function(x) induced_subgraph(graph, V(graph)$name[V(graph)$M==M[x]]))
		names(gHC) <- paste0("HM", M)
		gLM <- NULL
		
	} else if (HM == "LV") {
		ftm <- data.frame(from = c(paste0("LX", membership)),
			              to = names(membership))
		gLM <- graph_from_data_frame(ftm, directed = TRUE)
		#V(gLM)$name
		V(gLM)$LV <- 0
		V(gLM)$LV[1:K] <- 1
		V(gLM)$color <- ifelse(V(gLM)$LV == 1, "lightblue", "yellow")
		# plot(gLM)
		gHC <- NULL
		
	} else if (HM == "CV") {
		ftm <- data.frame(from = names(membership),
			              to = c(paste0("CY", membership)))
		gLM <- graph_from_data_frame(ftm, directed = TRUE)
		#V(gLM)$name
		V(gLM)$LV <- 0
		V(gLM)$LV[(vcount(gLM) - K + 1):vcount(gLM)] <- 1
		V(gLM)$color <- ifelse(V(gLM)$LV == 1, "lightblue", "green")
		# plot(gLM)
		gHC <- NULL
		
	} else if (HM == "none") {
		return(membership)
	}
	
	if (verbose & HM != "UV") {
		plot(gLM)
		Sys.sleep(5)
	}
	#readline(prompt = "Press [enter] to continue")
	
	return(list(gHM = gLM, M = membership, gHC = gHC))
}

#' @title Network extension through covariance analysis
#'
#' @description SEMextend extends an input network using significant 
#' residual covariances (derived from \code{\link[SEMgraph]{SEMbap}}) 
#' between three types of nodes: (i) connectors, being nodes with in- and 
#' out-degree greater than 0, (ii) sources, with in-degree = 0, and 
#' (iii) targets, with out-degree = 0. Network extension is performed in 
#' three ways: (i) new external nodes connecting pairs of nodes with 
#' significant residual covariances; (ii) new external parent nodes 
#' connecting source nodes; (iii) new external child nodes connecting 
#' target nodes. Parent and child nodes are found using de-sparsified 
#' graphical lasso method (D-S_GL), from \code{\link[SILGGM]{SILGGM}} 
#' function. Added nodes (connectors, regressors, and response variables) 
#' can already be present in the input graph. The extended output network 
#' will be fitted using \code{\link[SEMgraph]{SEMfit}}.
#' @param fit A fitted SEM object of class \code{\link{lavaan}} produced 
#' by \code{\link[SEMgraph]{SEMbap}}.
#' @param data A matrix with rows corresponding to subjects, and columns 
#' to graph nodes.
#' @param gnet External interaction network as an igraph object. Interaction 
#' data from this network will be used to integrate additional interaction 
#' information inside the model. Two preset databases are available: 
#' (i) kegg, for KEGG signaling pathways (directed), and (ii) string, 
#' for STRING protein interactions (undirected). Note that, if the input 
#' graph is undirected, an undirected reference interactome is expected. 
#' In case a directed interactome is used with an undirected input, a 
#' directed output network will be enforced.
#' @param B Node-node interaction fixed weight. If B is NULL (default), 
#' beta coefficients will be estimated by MLE. If B is numeric, it will 
#' be used as a scaling factor for the edge weights in the graph object 
#' (graph attribute E(graph)$weight). Since SEMgraph scales data before 
#' model fitting, we suggest a grid search for the optimal B value in the 
#' interval [0, 0.3]. As a rule of thumb, to our experience, B = 0.1 
#' performs well on any network.
#' @param perm Number of permutations. By default, perm is set to 0 and 
#' conventional standard errors will be computed. If perm > 1, P-values 
#' will be computed from a moment-based chi-squared approximation derived 
#' from the empirical distribution of permuted data (Larson and Owen, 2015). 
#' To reduce computational time costs per permutation, we suggest perm = 500 
#' (this will leave P-values precision almost unaltered). If perm = 1, 
#' no P-values are calculated.
#' @param alpha Significance level used for GGM search of common regressors 
#' and common response variables. By default, alpha is set to 0.05. 
#' As a general rule, to limit the number of added nodes with large networks, 
#' we suggest to set alpha as 1/df, where the degrees of freedom 
#' df = vcount(graph)*(vcount(graph)-1)/2 - ecount(graph). This will limit 
#' the number of new nodes proportionally to the graph sparsity.
#' @param verbose A logical value. If FALSE (default), the processed graphs 
#' will not be plotted to screen, saving execution time (they will be 
#' returned anyway).
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @import lavaan
#' @import SILGGM
#' @importFrom cate factor.analysis
#' @export
#'
#' @references
#' 
#' Grassi M & Palluzzi F (in preparation). SEMgraph: An R Package for 
#' Pathway and Network Analysis of Genomics Data with Structural Equation 
#' Models (SEM). Journal of Statistical Software (xxxx) 
#' 
#' Zhang R, Ren Z, Chen W (2018). SILGGM: An extensive R package for 
#' efficient statistical inference in large-scale gene networks. 
#' PLoS Comput. Biol., 14(8): e1006369. 
#' https://doi.org/10.1371/journal.pcbi.1006369
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "fit", SEM fitted lavaan object;
#' \item "gest", group effect estimates and P-values on graph nodes;
#' \item "model", SEM model as a string;
#' \item "graph", list of four igraph objects, contianing:
#' \itemize{
#' \item "duv", the extended directed (or bidirected) network with new 
#' nodes from the external interactome,
#' \item "gUX", extended source network with unmeasured variables,
#' \item "gUY", extended target network with unmeasured variables,
#' \item "Ug", union of duv, gUX, and gUY graphs;
#' }
#' \item "dataXY", input data subset mapping graph nodes, plus group at 
#' the first column (if no group is specified, this column will take NA 
#' values).
#' }
#' @seealso \code{\link{igraph}}, \code{\link[lavaan]{lavaan}}, 
#' \code{\link{SILGGM}}
#' 
#' @examples
#' # Specifying data groups
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert 
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' 
#' # Network fitting
#' fit <- SEMfit(graph, data, group, B = NULL, perm = 10000)
#' # Network degrees of freedom
#' ndf <- vcount(graph)*(vcount(graph) - 1)/2 - ecount(graph)
#' 
#' # Bow-free interaction search through GGMs
#' ggm <- SEMggm(fit = fit, gnet = kegg, d = 2, perm = 10000,
#'               alpha = 1/ndf,
#'               verbose = FALSE)
#' 
#' # Network extension
#' ext <- SEMext(fit = ggm, gnet = kegg, data = data, B = NULL, d = 2,
#'               perm = 10000,
#'               alpha = 1/ndf,
#'               verbose = TRUE)
#' 
#' # Results summary
#' summary(ext$gest)
#' pval <- ext$gest@res$pchisq[-c(1:3)]
#' length(which(p.adjust(pval, method = "BH") < 0.1))
#' 
#' # Plot output graphs
#' par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
#' W <- ext$graph$guv; E(W)$weight <- 1       # connectors
#' X <- ext$graph$gux                         # sources
#' Y <- ext$graph$guy                         # targets
#' # U = input network + C + X + Y            # union
#' U <- ext$graph$Ug; E(U)$weight <- 1
#' plot(W, main = "Extended connector network")
#' plot(X, main = "Extended source network")
#' plot(Y, main = "Extended target network")
#' plot(U, main = "Extended input network")
#' 
#' # SEM fitting of the extracted graphs (guv, gux, guy)
#' sem.w <- SEMfit(graph = W, data = data, group = group, B = NULL, perm = 0)
#' summary(sem.w$fit)
#' sem.x <- SEMfit(graph = X, data = data, group = group, B = NULL, perm = 0)
#' summary(sem.x$fit)
#' sem.y <- SEMfit(graph = Y, data = data, group = group, B = NULL, perm = 0)
#' summary(sem.y$fit)
#'
SEMext<- function(fit, data, gnet, d, B=NULL, perm=0, alpha=0.05, verbose=FALSE, ...)
{
	# Set SEM (group) and graph (ig, guu, gnet) objects :
	group<- fit$dataXY[,1]
	if( is.na(group[1]) ) group<- NULL
	ig<- fit[[4]]$ig
	guu<- fit[[4]]$guu
	vids<- which(V(gnet)$name %in% colnames(data))
	gnet<- induced_subgraph(graph=gnet, vids= vids)
	SET1<- as_edgelist(guu)
	if( nrow(SET1) == 0 ){
	 cat("NULL SEM fitting: No.covariances=0 !","\n\n")
	 return(list(fit=NULL, gest=NULL, model=NULL, graph=NULL, dataXY=NULL))
	}	
	
	# Search of bow-free acyclic covariances from interactome
	if( is.directed(gnet) ) ugnet<- as.undirected(gnet, mode="collapse")
	ftm1<- NULL
	for(j in 1:nrow(SET1)) { #j=14
	 cat("\r","edge=", j, "of", nrow(SET1))
     #Sys.sleep(0.01)
	 flush.console()
	 a<- SET1[j,1]
	 b<- SET1[j,2]
	 ftm1<- rbind(ftm1, c(a,b))
	 v<- which(V(ugnet)$name %in% c(a,b))
	 if ( length(v) == 2 ) {
	  sp<- distances(ugnet, a, b, mode="all", weights=NA)
	  if(sp <= d) ftm1[j,]<- c(a,b) else ftm1[j,]<- c(NA,NA)
	 }else{ ftm1[j,]<- c(NA,NA) }
	} 
	ftm1<- na.omit(ftm1)
	cat("\n", "n.selected covariances:", nrow(SET1),
	    "n.imported from interactome:", nrow(ftm1), "\n\n")
			
	# Search of external nodes from interactome
	guv<- psi2guv(guu=guu, ig=ig, gnet=gnet, alpha, verbose=verbose)
	if( is.directed(gnet) ){
	 duv<- mergeGraph(g=list(guv), gref=ig, gnet=gnet)
	 gUX<- psi2gUX(duv=duv, ig=ig, gnet=gnet, data, alpha, verbose=verbose)
	 gUY<- psi2gUY(duv=duv, ig=ig, gnet=gnet, data, alpha, verbose=verbose)
	}else{ gUX<- gUY<- NA }
	
	# SEM fitting of the merged (expanded) graph Ug
	cat("\n")
	if( is.directed(gnet) ){
	 g<- list(ig, duv, gUX, gUY)
	 Ug<- mergeGraph(g=g[!is.na(g)], gref=ig, gnet=gnet)
	}else{
	 if( is.directed(ig) & !is.na(guv)[1] ) guv<- as.directed(guv, mode="mutual")
	 Ug<- graph.union(g=list(ig, guv)[!is.na(list(ig, guv))])
	}
	fit2<- SEMfit(graph=Ug, data=data, group=group, B=B, perm=perm)
	#summary(fit2$fit)
	graph<- list(guv=guv, gux=gUX, guy=gUY, Ug=Ug)
	
	return( list(fit=fit2$fit, gest=fit2$gest, model=fit2$model, graph=graph, dataXY=fit2$dataXY) )
}

#' @title Hidden module scoring and fitting via Factor Score Regression (FSR)
#'
#' @description Generate factor scores, principal component scores, or 
#' projection scores of latent, composite, and unmeasured variable modules, 
#' respectively, and fit them in a SEM with exogenous group effect.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and 
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the number 
#' of subjects. Each vector element must be 1 for cases and 0 for control 
#' subjects.
#' @param HM Hidden model type. For each defined hidden module: 
#' (i) if HM = "LV", a latent variable (LV) will be defined as common 
#' unknown cause acting on cluster nodes; (ii) if HM = "CV", cluster nodes 
#' will be considered as regressors of a latent composite variable (CV); 
#' (iii) if HM = "UV", an unmeasured variable (UV) model will be generated 
#' for each latent module, where source nodes (i.e., in-degree = 0) act 
#' as common regressors influencing the other nodes via an unmeasured variable.
#' @param size Minimum number of nodes per hidden module. By default, a 
#' minimum number of 3 nodes is required.
#' @param type Graph clustering method, from \code{\link{igraph}} R package 
#' (see also \code{\link[SEMgraph]{clusterGraph}}).
#' @param verbose A logical value. If FALSE (default), the processed graphs 
#' will not be plotted to screen, saving execution time (they will be 
#' returned anyway).
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @import lavaan
#' @importFrom cate factor.analysis
#' @import SILGGM
#' @export
#'
#' @references
#' 
#' Hoshino, T., & Bentler, P.M. (2013). Bias in factor score regression 
#' and a simple solution. In de Leon, A.R., & Chough, K.C. (Eds.). 
#' Analysis of Mixed Data: Methods & Applications. 
#' New York: Chapman and Hall/CRC
#' 
#' Bai, J; Li, K (2012). Statistical analysis of factor models of high 
#' dimension. Ann. Statist. 40 (2012), no. 1, 436-465. 
#' doi:10.1214/11-AOS966.
#' 
#' Davies PT & Tso M K-S (1982). Procedures for Reduced-Rank Regression. 
#' Journal of the Royal Statistical Society. Series C (Applied Statistics), 
#' Vol. 31, No. 3, pp. 244-255
#' 
#' @return A list of 3 objects:
#' \enumerate{
#' \item "fit", hidden module fitting as a lavaan object;
#' \item "M", hidden module nodes membership;
#' \code{\link[SEMgraph]{clusterGraph}} function;
#' \item "dataHM", hidden module data matrix.
#' }
#' @seealso \code{\link[cate]{factor.analysis}}
#' @examples
#' # Data loading
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#' 
#' # Module finding and Factor Score Regression
#' fsr.uv <- SEMfsr(graph = graph, data = data, group = group,
#'                  type = "ebc",
#'                  HM = "LV",
#'                  size = 15,
#'                  verbose = TRUE)
#' 
#' # Output summary
#' summary(fsr.uv$fit)
#' # Group membership
#' table(fsr.uv$M)
#' # Module scores
#' head(fsr.uv$dataHM)
#' 
#' # Hidden modules
#' gHC <- fsr.uv$gHC
#' gHC
#' gplot(x = gHC$g)
#' 
#' # Hidden modules reduction and fitting pipeline
#' fit <- SEMfit(graph, data, group, B = NULL, perm = 10000)
#' ggm <- SEMggm(fit = fit, gnet = kegg, d = 2,
#'               perm = 10000,
#'               alpha = 5E-05,
#'               verbose = FALSE)
#' ext <- SEMext(fit = ggm, gnet = kegg, data = data, B = NULL, d = 2,
#'               perm = 10000,
#'               alpha = 5E-05,
#'               verbose = FALSE)
#' cg <- mergeNodes(graph = ext$graph$Ug, membership = fsr.uv$M)
#' properties(cg)
#' 
#' # Adding hidden module data to the original data
#' data2 <- cbind(fsr.uv$dataHM[, -1], data)
#' head(data2)[, 1:10]
#' fit.hm <- SEMfit(graph = cg, data = data2, group = group, perm = 10000)
#' summary(fit.hm$fit)
#' summary(fit.hm$gest)
#'
SEMfsr <- function(graph, data, group, type, HM, size = 3, verbose = FALSE, ...)
{  
	# Set SEM objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	ig <- simplify(ig, remove.loops = TRUE)
	
	# Hidden modules LX -> Y
	if (HM == "LV") {
		LX <- clusterGraph(graph = ig, type = type,
		               HM = "LV",
		               size = size,
		               verbose = verbose)
		if (length(LX) == 0) return(list(fit = NA, M = NA, dataHM = NA))
		gLM <- LX[[1]]
		membership <- LX[[2]]
		gLC <- LX[[3]]
		LX <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "L"]
		
		# Latent Variables(LV) model
		K <- as.numeric(names(table(membership)))
		LV <- NULL
		for(k in 1:length(LX)) { #k=1
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[, which(colnames(dataY) %in% Xk)])
			fa1 <- cate::factor.analysis(Y = Y, r = 1, method = "ml")$Z
			LV <- cbind(LV, fa1)
		}
		M <- paste0("LM", membership)
		names(M) <- names(membership)
		colnames(LV) <- gsub("X", "M", LX)
		rownames(LV) <- rownames(dataY)
		dataLC <- cbind(group, LV)
		#head(dataLC)
		# group mean differences effects
		model <- paste0(colnames(LV), "~group")
	}
	#cat(model)
	
	# Hidden modules X -> LY
	if (HM == "CV") {
		LY <- clusterGraph(graph = ig, type = type,
		               HM = "CV",
		               size = size,
		               verbose = verbose)
		if (length(LY) == 0) return(list(fit = NA, M = NA, dataHM = NA))
		gLM <- LY[[1]]
		membership <- LY[[2]]
		gLC <- LY[[3]]
		LY <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "C"]
		
		# Composite Variables(CV) model
		K <- as.numeric(names(table(membership)))
		CV <- NULL
		for(k in 1:length(LY)) { #k=1
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[,which(colnames(dataY) %in% Xk)])
			pc1 <- cate::factor.analysis(Y = Y, r = 1, method = "pc")$Z
			CV <- cbind(CV, pc1)
	}
	M <- paste0("CM", membership)
	names(M) <- names(membership)
	colnames(CV) <- gsub("CY", "CM", LY)
	rownames(CV) <- rownames(dataY)
	dataLC <- cbind(group, CV)
	#head(dataLC)
	# group mean differences effects
	model <- paste0(colnames(CV), "~group")
	}
	#cat model
	
	# Hidden modules X -> UV -> Y
	if (HM == "UV") {
		if (!is.directed(graph)) {
			return(cat("\nUV is not applicable with udirected graph\n"))
		}
		LXY <- clusterGraph(graph = ig, type = type,
		                HM = "UV",
		                size = size,
		                verbose = verbose)
		if(length(LXY) == 0) return(list(fit = NA, M = NA, dataHM = NA))
		gLM <- NULL
		membership <- LXY[[2]]
		gLC <- LXY[[3]]
		LXY <- paste0("HM", names(table(membership)))
		
		# Unmeasured Variables(UV) model
		UV <- na <- NULL
		for(k in 1:length(LXY)) { #k=1
			gk <- gLC[[which(names(gLC) %in% LXY)[k]]]
			#gplot(gk)
			d <- igraph::degree(gk, mode = "in")
			idx <- which(colnames(dataY) %in% V(gk)$name[d == 0])
			Xk <- as.matrix(dataY[, idx])
			#head(Xk)
			if(ncol(Xk) > nrow(Xk)) {
				na <- c(na, k)
				next
			}
			idy <- which(colnames(dataY) %in% V(gk)$name[d > 0])
			Yk <- as.matrix(dataY[, idy])
			#head(Yk)
			Uk <- Xk%*%solve(t(Xk)%*%Xk)%*%t(Xk)%*%Yk
			#head(Uk)
			spc1 <- cate::factor.analysis(Y = as.matrix(Uk), r = 1,
			                              method = "pc")$Z
			UV <- cbind(UV, spc1)
			#head(UV)
	 }
	 if(length(na) == 0) colnames(UV) <- LXY else colnames(UV) <- LXY[-na]
	 M <- paste0("UM", membership)
	 names(M) <- names(membership)
	 colnames(UV) <- gsub("HM", "UM", LXY)
	 rownames(UV) <- rownames(dataY)
	 dataLC <- cbind(group, UV)
	 # head(dataLC)
	 # group mean differences effects
	 model <- paste0(colnames(UV), "~group")
	}
	#cat model
	
	if (length(group) > 0) {
		fsr <- sem(model, data = dataLC, se = "standard", fixed.x = TRUE)
		if (fsr@Fit@converged == TRUE) {
			srmr <- fitMeasures(fsr, c("srmr"))
			cat("Model converged:", fsr@Fit@converged, "srmr:", srmr, "\n\n")
		} else {
			cat("Model converged:", fsr@Fit@converged, "srmr:", NA, "\n\n")
			fsr<- NULL
		}
		#summary(fsr)
	} else if (length(group) == 0) {
		fsr <- NULL
		dataLC <- cbind(group = rep(NA, nrow(dataY)), dataLC)
	}
	
	return(list(fit = fsr, M = M, dataHM = dataLC))
}

