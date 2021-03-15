#  SEMgraph library
#  Copyright (C) 2019-2021 Fernando Palluzzi; Mario Grassi
#  e-mail: <fernando.palluzzi@gmail.com>
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
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #


diagonalizePsi <- function(g = list(graph, guu), data, ...)
{
	# Set graph and data objects
	graph <- g[[1]]
	V <- colnames(data)[colnames(data) %in% V(graph)$name]
	Y <- scale(data[, V])
	graph <- induced_subgraph(graph, vids = which(V(graph)$name %in% V))
	A0 <- as_adj(as.undirected(graph), type = "both", sparse = FALSE)[V, V]
	
	# Precision fitting of guu -> wi
	guu <- g[[2]]
	adj <- as_adj(guu, sparse = FALSE)
	idx <- which(rownames(A0) %in% rownames(adj) == FALSE)
	if (length(idx) > 0) {
		R <- matrix(0, length(idx), ncol(adj))
		C <- matrix(0, nrow(adj), length(idx))
		I <- diag(length(idx))
		adj <- rbind(cbind(I, R), cbind(C, adj))
		rownames(adj)[1:length(idx)] <- rownames(A0)[idx]
		colnames(adj)[1:length(idx)] <- rownames(A0)[idx]
	}
	Sigma <- cor(Y[, colnames(adj)])
	wi <- GGMncv::constrained(Sigma, adj)$Theta
	colnames(wi) <- rownames(wi) <- colnames(adj)
	if (!corpcor::is.positive.definite(wi)) {
		wi <- corpcor::cor.shrink(wi, verbose = FALSE)
		#wi <- corpcor::cov.shrink(wi, verbose = TRUE)
		#wi <- corpcor::make.positive.definite(wi)
	}
	E <- eigen(wi) # Eigenvalues and eigenvectors of w
	R <- E$vectors%*%diag(sqrt(E$values))%*%t(E$vectors)
	#sum(wi - R %*% R)
	Y <- Y[,colnames(wi)]
	YR <- as.matrix(Y)%*%R
	colnames(YR) <- colnames(Y)
	
	return(data = YR[, V])
}

psi2guv <- function(guu, ig, gnet, verbose, ...)
{
	# Adding new nodes to undirected graph (guu) from the interactome
	vids <- which(V(guu)$name %in% V(gnet)$name)
	guu <- induced_subgraph(graph = guu, vids = vids)
	if(verbose) {
		plot(guu, main = "direct(or covariance) graph (guu) in gnet")
		Sys.sleep(3)
	}
	ftm <- as_edgelist(guu)
	vpath <- ftmuv <- NULL
	
	for (i in 1:nrow(ftm)) {
		
		mode <- ifelse(is.directed(guu) & is.directed(gnet), "out", "all")
		
		if (distances(gnet, ftm[i, 1], ftm[i, 2], mode = mode,
		              weights = NA) == Inf) next
		
		if (is.null(E(gnet)$pv)) {
			suppressWarnings(path <- shortest_paths(gnet, ftm[i, 1],
			                                        ftm[i, 2],
			                                        mode = mode,
			                                        weights = NA)$vpath)
		} else {
			path <- all_shortest_paths(gnet, ftm[i, 1], ftm[i, 2],
			                           mode = mode,
			                           weights = NA)$res
		}
		
		if (length(path) > 1) {
			fX2 <- NULL
			for (k in 1:length(path)) {
				pathk <- induced_subgraph(gnet, V(gnet)$name[path[[k]]])
				fX2[k] <- -2*sum(log(E(pathk)$pv))
			}
			path <- path[[which(fX2 == max(fX2))[1]]]
		} else {
			path <- path[[1]]
		}
		V <- V(gnet)$name[path]
		vpath <- c(vpath, V[-c(1, length(V))])
		for(h in 1:(length(V) - 1)) ftmuv <- rbind(ftmuv, c(V[h], V[h + 1]))
	}
	
	# Graph with added nodes from interactome (guv)
	ftmuv <- na.omit(ftmuv[duplicated(ftmuv) != TRUE,])
	ftmuv <- matrix(ftmuv, ncol = 2)
	
	if (nrow(ftmuv) > 0) {
		mode <- ifelse(is.directed(guu) & is.directed(gnet), TRUE, FALSE)
		guv <- graph_from_data_frame(ftmuv, directed = mode)
		guv <- simplify(guv, remove.loops = TRUE)
		vv <- V(guv)$name[-which(V(guv)$name %in% V(ig)$name)]
		uv <- V(ig)$name[which(V(ig)$name %in% unique(vpath))]
		V(guv)$color[V(guv)$name %in% V(guu)$name] <- "deepskyblue"
		V(guv)$color[V(guv)$name %in% vv] <- "gold"
		V(guv)$color[V(guv)$name %in% uv] <- "green2"
		
		if(verbose) {
			plot(guv, main = "Extended connector graph (guv)")
			Sys.sleep(3)
		}
	
	} else {
		cat("\n", "no edges u->u (or u--v) found !", "\n\n")
		guv <- make_empty_graph(n = 0)
	}
	return(guv)
}

#' @title Subgraph mapping
#'
#' @description Map groups of nodes onto an input graph, based on a 
#' membership vector.
#' @param graph An igraph object.
#' @param membership Cluster membership vector for each node.
#' @param l graph layout. One of the \code{\link{igraph}} layouts. 
#' If this argument is ignored, an automatic layout will be applied.
#' @param map A logical value. Visualize cluster mapping over the input 
#' graph. If FALSE (default), visualization will be disabled. For large 
#' graphs, visualization may take long.
#' @param verbose A logical value. If FALSE (default), the processed 
#' graphs will not be plotted to screen, saving execution time (they will 
#' be returned in output anyway).
#' @param ... Currently ignored. 
#'
#' @import igraph
#' @importFrom stats qnorm cov2cor cor hclust as.dist cutree
#' @importFrom graphics abline
#' @export
#'
#' @return The list of clusters and cluster mapping as igraph objects.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @seealso \code{\link[SEMgraph]{clusterGraph}}, 
#' \code{\link[SEMgraph]{clusterScore}}
#'
#' @examples
#' library(SEMdata)
#' G <- kegg.pathways$"Amyotrophic lateral sclerosis (ALS)"
#' # Largest connected component
#' G <- properties(G)[[1]]
#' membership <- clusterGraph(graph = G, type = "wtc")
#' cplot(G, membership, map = TRUE)
#' 
#' \dontrun{
#' cplot(G, membership, map = FALSE, verbose = TRUE)
#' }
#'
cplot <- function(graph, membership, l = layout.auto, map = FALSE,
                  verbose = FALSE, ...)
{
	# Overall cluster visualization
	V(graph)$M <- 9999
    V(graph)$M[which(V(graph)$name %in% names(membership))] <- membership
	if (map) {
		V(graph)$color <- V(graph)$M + 1
		gplot(graph)
		Sys.sleep(3)
	}
	
	# Within cluster visualization
	M <- names(table(V(graph)$M))
	K <- length(table(V(graph)$M))
	vcol <- as.numeric(M) + 1
	
	HM <- lapply(1:K, function(x) induced_subgraph(graph,
	             V(graph)$name[V(graph)$M == M[x]]))
	
	names(HM) <- paste0("HM", M)
	d <- igraph::degree(graph, mode = "all")*2 + 1
	
	if (verbose) {
		glv <- lapply(1:K, function(x) {
			          E(HM[[x]])$weight <- 1
			          plot(HM[[x]],
			          vertex.color = vcol[x],
			          vertex.size = d[V(HM[[x]])$name],
			          layout = l,
			          main = paste0("Hidden Module ", M[x]))
		Sys.sleep(3)})
	}
	
	return(invisible(c(list(graph = graph), HM)))
}

#' @title Graph nodes merging by a user-defined membership attribute
#'
#' @description Merge groups of graph nodes using a custom membership
#' attribute (e.g., cluster membership).
#' @param graph Network as an igraph object.
#' @param membership Cluster membership. A vector of cluster membership
#' identifiers, where vector names correspond to graph node names.
#' Topological graph clustering can be done using 
#' \code{\link[SEMgraph]{clusterGraph}}.
#' @param HM Hidden model label. If HM = "LV", a latent variable (LV) 
#' will be defined as common unknown cause acting on cluster nodes.
#' If HM = "CV", cluster nodes will be considered as regressors of a 
#' latent composite variable (CV). Finally, if HM = "UV", an unmeasured 
#' variable (UV) is defined, where source nodes of the module (i.e., 
#' in-degree = 0) act as common regressors influencing the other nodes 
#' via an unmeasured variable.
#' @param ... Currently ignored.
#'
#' @import igraph
#' @importFrom graph combineNodes
#' @export
#'
#' @return A network with merged nodes as an igraph object.
#' @seealso \code{\link[SEMgraph]{clusterGraph}}
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#' library(SEMdata)
#' G <- kegg.pathways$"Amyotrophic lateral sclerosis (ALS)"
#' # Largest connected component
#' G <- properties(G)[[1]]
#' membership <- clusterGraph(graph = G, type = "wtc")
#' M <- mergeNodes(G, membership, HM = "LV")
#' gplot(M)
#'
mergeNodes <- function(graph, membership, HM, ...)
{
	# Set membership object
	if (is.numeric(membership)) {
		nodes <- names(membership)
		membership <- paste0(HM, membership)
		names(membership) <- nodes
	}

	LM <- NULL
	for (i in 1:length(table(membership))) { #i=1
		m <- names(table(membership))[i]
		LMi <- V(graph)$name[which(V(graph)$name %in% names(membership)[membership == m])]
		LM <- c(LM, list(LMi))
	}
	names(LM) <- names(table(membership))
	
	# Visualize graph object
	gLM <- as_graphnel(graph)
	for (i in 1:length(LM)) {
		gLMi <- graph::combineNodes(LM[[i]], gLM, names(LM)[i], mean)
		gLM <- gLMi
	}
	
	ig <- graph_from_graphnel(gLM)
	if (length(V(ig)$color) == 0) V(ig)$color <- "white"
	V(ig)$color[substr(V(ig)$name, 2, 2) == "V"] <- "orange"
	vcol <- V(ig)$color
	names(vcol) <- V(ig)$name
	gplot(ig)
	
	return(gLM = ig)
}

#' @title Topological graph clustering
#'
#' @description Topological graph clustering methods.
#' @param graph An igraph object.
#' @param type Topological clustering methods. If type = "tahc", network 
#' modules are generated using the tree agglomerative hierarchical 
#' clustering method (Yu et al., 2015).
#' Other non-tree clustering methods from igraph package include: "wtc" 
#' (default value; walktrap community structure with short random walks), 
#' "ebc" (edge betweeness clustering), "fgc" (fast greedy method), "lbc" 
#' (label propagation method), "lec" (leading eigenvector method), "loc" 
#' (multi-level optimization), "opc" (optimal communiy structure), "sgc" 
#' (spinglass statistical mechanics).
#' @param HM Hidden model type. Enables the visualization of the hidden 
#' model. If set to "none" (default), no HM is visualized.
#' For each defined hidden module:
#' (i) if HM = "LV", a latent variable (LV) will be defined as common
#' unknown cause acting on cluster nodes; (ii) if HM = "CV", cluster nodes
#' will be considered as regressors of a latent composite variable (CV);
#' (iii) if HM = "UV", an unmeasured variable (UV) is defined, where source
#' nodes of the module (i.e., in-degree = 0) act as common regressors
#' influencing the other nodes via an unmeasured variable (see also
#' \code{\link[SEMgraph]{clusterScore}}).
#' @param size Minimum number of nodes per module. By default, a minimum 
#' number of 5 nodes is required. 
#' @param verbose A logical value. If FALSE (default), the processed graphs 
#' will not be plotted to screen, saving execution time (they will be 
#' returned in output anyway). 
#' @param ... Currently ignored. 
#'
#' @import igraph
#' @importFrom stats qnorm cov2cor cor hclust as.dist cutree
#' @importFrom graphics abline
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
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
#' \item "gHM", subgraph containing hidden modules as an igraph object;
#' \item "membership", cluster membership vector for each node;
#' \item "gHC", the list of modules as igraph objects.
#' }
#' If HM is "none", only the cluster membership vector is returned.
#'
#' @seealso \code{\link[SEMgraph]{clusterScore}}, \code{\link[SEMgraph]{cplot}}
#'
#' @examples
#' library(SEMdata)
#' G <- kegg.pathways$"Amyotrophic lateral sclerosis (ALS)"
#' # Largest connected component
#' G <- properties(G)[[1]]
#' membership <- clusterGraph(graph = G, type = "wtc", HM = "LV", verbose = TRUE)
#'
clusterGraph <- function(graph, type = "wtc", HM = "none", size = 5,
                         verbose = FALSE, ...)
{
    # Set undirected igraph object
	if (!is_directed(graph)) {
		ug <- graph
	} else {
		ug <- as.undirected(graph, mode = "collapse",
		                    edge.attr.comb = "ignore")
	}
	
	if (type == "tahc") {
		# Tree Agglomerative Hierarchical Clustering (TAHC)
		mst <- minimum.spanning.tree(ug, weights = NULL, algorithm = NULL)
		G <- distances(mst, v = V(mst), to = V(mst), mode = "all",
		               weights = NA)
		D <- 1 - cor(x = G, method = "spearman")
		hMST <- hclust(as.dist(D), method = "average")
		tahc <- cutree(hMST, h = 0.2)
		cnames <- as.numeric(names(table(tahc)))[table(tahc) >= size]
		membership <- tahc[tahc %in% cnames]
		if(verbose) {
			plot(hMST, labels = FALSE, xlab = "", sub = "")
			abline(h = 0.2, col = "red")
			Sys.sleep(3)
		}
	
	} else {
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
			Sys.sleep(3)
		}
	}
	
	K <-  length(cnames)
	if (K == 0) return(cat("WARNING: no communities with size >=", size, "\n"))
	if (HM == "UV") {
		gHC <- cplot(graph, membership = membership, map = FALSE,
		             verbose = FALSE)[-1]
		ftm <- Vxx <- NULL
	
	for (i in 1:K) {
		d <- igraph::degree(gHC[[i]], mode = "in")
		Vx <- V(gHC[[i]])$name[d == 0]
		Vy <- V(gHC[[i]])$name[d != 0]
		ftm <- rbind(ftm, cbind(Vx, rep(paste0("UV", i), length(Vx))))
		ftm <- rbind(ftm, cbind(rep(paste0("UV", i), length(Vy)), Vy))
		Vxx <- c(Vxx, Vx)
	}
	gLM <- graph_from_data_frame(ftm, directed = TRUE)
	V(gLM)$color <- "yellow"
	V(gLM)$color[substr(V(gLM)$name, 1, 1) == "U"] <- "lightblue"
	V(gLM)$color[V(gLM)$name %in% Vxx] <- "green"
	
	} else if (HM == "LV") {
		ftm <- data.frame(from = c(paste0("LX", membership)),
		                  to = names(membership))
		gLM <- graph_from_data_frame(ftm, directed = TRUE)
		V(gLM)$LV <- 0
		V(gLM)$LV[1:K] <- 1
		V(gLM)$color <- ifelse(V(gLM)$LV == 1, "lightblue", "yellow")
		gHC <- NULL
	
	} else if (HM == "CV") {
		ftm <- data.frame(from = names(membership),
		                  to = c(paste0("CY", membership)))
		gLM <- graph_from_data_frame(ftm, directed = TRUE)
		V(gLM)$LV <- 0
		V(gLM)$LV[(vcount(gLM) - K + 1):vcount(gLM)] <- 1
		V(gLM)$color <- ifelse(V(gLM)$LV == 1, "lightblue", "green")
		gHC <- NULL
	
	} else if (HM == "none") {
		return( membership )
	}
	
    if (verbose == TRUE) {
		plot(gLM)
	}
	
	return(list(gHM = gLM, membership = membership, gHC = gHC))
}

#' @title Interactome-assisted graph extension
#'
#' @description Extend an input directed graph, importing new interactions 
#' from a second graph. Added interactions will be chosen among those 
#' available in a given reference interactome.
#'
#' @param g A list of two graphs as igraph objects.
#' @param data A matrix with rows corresponding to subjects, and columns
#' to graph nodes.
#' @param gnet External interaction network as an igraph object. Interaction
#' data from this network will be used to integrate additional interaction
#' information inside the graph.
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen, saving execution time (they will be
#' returned anyway).
#' @param ... Currently ignored.
#'
#' @details This function takes two input graphs: the first is the input 
#' causal model (i.e., a directed graph), and the second can be either 
#' a directed or undirected graph, providing a set of connections to be 
#' checked against the reference network and imported to the first graph. 
#' Typically, the second graph is the output of either 
#' \code{\link[SEMgraph]{SEMdag}} or \code{\link[SEMgraph]{SEMbap}}. 
#' In the former we use the new inferred causal structure stored in the 
#' \code{dag.red} object. In the latter, we use the new inferred covariance 
#' structure stored in the \code{guu} object. In both cases, new hidden 
#' directed paths and new nodes (i.e., new mediators) can be revealed.
#'
#' @import igraph
#' @import lavaan
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' Grassi M, Palluzzi F (2021). SEMgraph: An R Package for Causal Network 
#' Analysis of High-Throughput Data with Structural Equation Models. 
#' xxxxx x(x): xxxxx. https://doi.org/xxxxx
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "Ug", the extended graph (union of the input graph and guv);
#' \item "guv", the directed subgraph added to the input graph.
#' }
#'
#' @examples
#' library(SEMdata)
#' G <- kegg.pathways$"Steroid biosynthesis"
#' G <- properties(G)[[1]]
#' 
#' # Extend a graph using new inferred DAG edges
#' 
#' library(SEMdata)
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' dag <- SEMdag(graph = G, data = als.npn, beta = 0.1)
#' ext <- extendGraph(list(dag$dag, dag$dag.red), data = als.npn, gnet = kegg)
#' gplot(ext$Ug)
#' 
#' # Extend a graph using the inferred bow-free path diagram
#' 
#' bap <- SEMbap(graph = G, data = als.npn, gnet = kegg, d = 1, alpha = 0.05)
#' ext <- extendGraph(list(bap$bap, bap$guu), data = als.npn, gnet = kegg)
#' gplot(ext$Ug)
#'  
#' # Create a graph from correlation matrix, using KEGG as reference
#' 
#' v <- which(colnames(als.npn) %in% V(G)$name)
#' selectedData <- als.npn[, v]
#' G0 <- make_empty_graph(n = ncol(selectedData))
#' V(G0)$name <- colnames(selectedData)
#' 
#' G1 <- corr2graph(R = cor(selectedData), n = nrow(selectedData),
#'                  type = "tmfg")
#' ext <- extendGraph(list(G0, G1), data = selectedData, gnet = kegg)
#' par(mfrow=c(1,2), mar=rep(1,4))
#' plot(G1, layout = layout.circle)
#' plot(ext$Ug, layout = layout.circle)
#'
extendGraph <- function(g = list(), data, gnet, verbose = FALSE, ...)
{
	# Set graph (ig, guu, gnet) objects
	graph <- g[[1]]
	if (!is_directed(graph)) {
		cat("ERROR: The first input graph is not a directed graph.\n")
		return(NULL)
	}
	guu <- g[[2]]
	vids <- which(V(gnet)$name %in% colnames(data))
	gnet <- induced_subgraph(graph = gnet, vids = vids)
	vids <- which(V(graph)$name %in% colnames(data))
	graph <- induced_subgraph(graph = graph, vids = vids)
	ig <- graph - E(graph)[E(graph)$color == "red"]
	if (!is.null(E(ig)$weight)) ig <- delete_edge_attr(ig, "weight")
	if (!is.null(E(ig)$color)) ig <- delete_edge_attr(ig, "color")
	if (!is.null(V(ig)$color)) ig <- delete_vertex_attr(ig, "color")
	
	# Search external nodes from interactome
	guv <- psi2guv(guu = guu, ig = ig, gnet = gnet, verbose = verbose)
	if (ecount(guv) == 0) return(list(Ug = ig, guv = guv))
	
	# Union graph
	
	if (is.directed(guv) & is.directed(gnet)) {
		Ug <- graph.union(g = list(ig, guv))
	}
	
	if (!is.directed(guv) & is.directed(gnet)) {
		guv <- orientEdges(ug = guv, dg = gnet, data = NULL)
		Ug <- graph.union(g = list(ig, guv))
	}
	
	if (!is.directed(guv) & !is.directed(gnet)) {
		guv <- orientEdges(ug = guv, dg = NULL, data = data)
		Ug <- graph.union(g = list(ig, guv))
	}
	
	E1 <- attr(E(Ug), "vnames")
	E0 <- attr(E(ig), "vnames")
	E(Ug)$color <- ifelse(E1 %in% E0, "blue", "red")
	
	return(list(Ug = Ug, guv = guv))
}

#' @title Module scoring
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
#' for each module, where source nodes (i.e., in-degree = 0) act as common 
#' regressors influencing the other nodes via an unmeasured variable.
#' @param size Minimum number of nodes per hidden module. By default, a
#' minimum number of 5 nodes is required.
#' By default, HM is set to "LV" (i.e., the latent variable model).
#' @param type Graph clustering method. If type = "tahc", network 
#' modules are generated using the tree agglomerative hierarchical 
#' clustering method (Yu et al., 2015).
#' Other non-tree clustering methods from igraph package include: "wtc" 
#' (default value; walktrap community structure with short random walks), 
#' "ebc" (edge betweeness clustering), "fgc" (fast greedy method), "lbc" 
#' (label propagation method), "lec" (leading eigenvector method), "loc" 
#' (multi-level optimization), "opc" (optimal communiy structure), "sgc" 
#' (spinglass statistical mechanics).
#' By default, the "wtc" method is used.
#' @param verbose A logical value. If TRUE, intermediate graphs will be 
#' displayed during the execution. In addition, a condensed graph with 
#' clusters as nodes will be fitted and showed to screen (see also 
#' \code{\link[SEMgraph]{mergeNodes}}). By default, verbode = FALSE.
#' @param ... Currently ignored.
#'
#' @import igraph
#' @import lavaan
#' @importFrom cate factor.analysis
#' @export
#'
#' @seealso 
#' See \code{\link[SEMgraph]{clusterGraph}} and \code{\link[SEMgraph]{cplot}}
#' for graph clustering, and \code{\link[cate]{factor.analysis}} for 
#' factor analysis.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' Grassi M, Palluzzi F (2021). SEMgraph: An R Package for Causal Network 
#' Analysis of High-Throughput Data with Structural Equation Models. 
#' xxxxx x(x): xxxxx. https://doi.org/xxxxx
#'
#' @return A list of 3 objects:
#' \enumerate{
#' \item "fit", hidden module fitting as a lavaan object;
#' \item "membership", hidden module nodes membership;
#' \code{\link[SEMgraph]{clusterGraph}} function;
#' \item "dataHM", hidden module data matrix with cluster scores.
#' }
#'
#' @examples
#' 
#' library(SEMdata)
#' library(huge)
#' 
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' C <- clusterScore(graph = alsData$graph, data = als.npn,
#'                   group = alsData$group,
#'                   HM = "LV",
#'                   type = "wtc",
#'                   verbose = TRUE)
#' summary(C$fit)
#' head(C$dataHM)
#' table(C$membership)
#'
clusterScore <- function(graph, data, group, HM = "LV", type = "wtc",
                         size = 5, verbose = FALSE, ...)
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
		LX <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "L"]

		# Latent Variables(LV) model
		K <- as.numeric(names(table(membership)))
		LV <- NULL
		for(k in 1:length(LX)) {
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[, which(colnames(dataY) %in% Xk)])
			fa1 <- cate::factor.analysis(Y = Y, r = 1, method = "ml")$Z
			LV <- cbind(LV, fa1)
		}
		colnames(LV) <- gsub("LX", "LV", LX)
		rownames(LV) <- rownames(dataY)
		dataLC <- cbind(group, LV)
		
		# Group mean differences effects
		model <- paste0(colnames(LV), "~group")
	}
	
	# Hidden modules X -> LY
	if (HM == "CV") {
		LY <- clusterGraph(graph = ig, type = type,
		               HM = "CV",
		               size = size,
		               verbose = verbose)
		if (length(LY) == 0) return(list(fit = NA, M = NA, dataHM = NA))
		gLM <- LY[[1]]
		membership <- LY[[2]]
		LY <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "C"]

		# Composite Variables(CV) model
		K <- as.numeric(names(table(membership)))
		CV <- NULL
		for(k in 1:length(LY)) {
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[,which(colnames(dataY) %in% Xk)])
			pc1 <- cate::factor.analysis(Y = Y, r = 1, method = "pc")$Z
			CV <- cbind(CV, pc1)
	}
	colnames(CV) <- gsub("CY", "CV", LY)
	rownames(CV) <- rownames(dataY)
	dataLC <- cbind(group, CV)
	
	# Group mean differences effects
	model <- paste0(colnames(CV), "~group")
	}
	
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
		membership <- LXY[[2]]
		gLC <- LXY[[3]]
		LXY <- paste0("HM", names(table(membership)))
		
		# Unmeasured Variables(UV) model
		UV <- na <- NULL
		for (k in 1:length(LXY)) {
			gk <- gLC[[which(names(gLC) %in% LXY)[k]]]
			d <- igraph::degree(gk, mode = "in")
			idx <- which(colnames(dataY) %in% V(gk)$name[d == 0])
			Xk <- as.matrix(dataY[, idx])
			idy <- which(colnames(dataY) %in% V(gk)$name[d > 0])
			if (ncol(Xk) > nrow(Xk) | length(idx) == 0 | length(idy) == 0) {
				na <- c(na, k)
				next
			}
			Yk <- as.matrix(dataY[, idy])
			Uk <- Xk%*%solve(t(Xk)%*%Xk)%*%t(Xk)%*%Yk
			spc1 <- cate::factor.analysis(Y = as.matrix(Uk), r = 1,
			                              method = "pc")$Z
			UV <- cbind(UV, spc1)
		}
		
		if (length(na) == 0) {
			colnames(UV) <- gsub("HM", "UV", LXY)
		} else {
			colnames(UV) <- gsub("HM", "UV", LXY[-na])
		}
		rownames(UV) <- rownames(dataY)
		dataLC <- cbind(group, UV)
		
		# Group mean differences effects
		model <- paste0(colnames(UV), "~group")
	}
	
	if (length(group) > 0) {
		fsr <- sem(model, data = dataLC, se = "standard", fixed.x = TRUE)
		if (fsr@Fit@converged == TRUE) {
			srmr <- fitMeasures(fsr, c("srmr"))
			cat("Model converged:", fsr@Fit@converged, "\nSRMR:", srmr, "\n\n")
		} else {
			cat("Model converged:", fsr@Fit@converged, "\nSRMR:", NA, "\n\n")
			fsr<- NULL
		}
	
	} else if (length(group) == 0) {
		fsr <- NULL
		dataLC <- cbind(group = rep(NA, nrow(dataY)), dataLC)
	}
	
	if (verbose == TRUE) {
		X <- cbind(dataLC, data)
		gM <- mergeNodes(graph, membership, HM = HM)
		gplot(gM)
		sem1 <- SEMfit(gM, X, group)
	}
	
	return(list(fit = fsr, membership = membership, dataHM = dataLC))
}
