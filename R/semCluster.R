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

#' @title Topological graph clustering
#'
#' @description Topological graph clustering methods.
#' @param graph An igraph object.
#' @param type Topological clustering methods. If \code{type = "tahc"},
#' network modules are generated using the tree agglomerative hierarchical
#' clustering method (Yu et al., 2015). Other non-tree clustering methods
#' from \code{\link{igraph}} package include: "wtc"
#' (default value; walktrap community structure with short random walks),
#' "ebc" (edge betweeness clustering), "fgc" (fast greedy method), "lbc"
#' (label propagation method), "lec" (leading eigenvector method), "loc"
#' (multi-level optimization), "opc" (optimal community structure), "sgc"
#' (spinglass statistical mechanics).
#' @param HM Hidden model type. Enables the visualization of the hidden
#' model, gHM. If set to "none" (default), no gHM igraph object is saved.
#' For each defined hidden module:
#' (i) if \code{HM = "LV"}, a latent variable (LV) will be defined as
#' common unknown cause acting on cluster nodes; (ii) if \code{HM = "CV"},
#' cluster nodes will be considered as regressors of a latent composite
#' variable (CV); (iii) if \code{HM = "UV"}, an unmeasured variable (UV)
#' is defined, where source nodes of the module (i.e., in-degree = 0)
#' act as common regressors influencing the other nodes via an unmeasured
#' variable (see also \code{\link[SEMgraph]{clusterScore}}).
#' @param size Minimum number of nodes per module. By default, a minimum
#' number of 5 nodes is required.
#' @param verbose A logical value. If FALSE (default), the gHM igraph
#' will not be plotted to screen, saving execution time (they will be
#' returned in output anyway).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Fortunato S, Hric D. Community detection in networks: A user guide (2016).
#' Phys Rep; 659: 1-44. <https://dx.doi.org/10.1016/j.physrep.2016.09.002>
#'
#' Yu M, Hillebrand A, Tewarie P, Meier J, van Dijk B, Van Mieghem P,
#' Stam CJ (2015). Hierarchical clustering in minimum spanning trees.
#' Chaos 25(2): 023107. <https://doi.org/10.1063/1.4908014>
#'
#' @return If HM is not "none" a list of 2 objects is returned:
#' \enumerate{
#' \item "gHM", subgraph containing hidden modules as an igraph object;
#' \item "membership", cluster membership vector for each node.
#' }
#' If HM is "none", only the cluster membership vector is returned.
#'
#' @seealso \code{\link[SEMgraph]{clusterScore}}, \code{\link[SEMgraph]{cplot}}
#'
#' @examples
#'
#' # Clustering ALS graph with WTC method and LV model
#' G <- properties(alsData$graph)[[1]]
#' clv <- clusterGraph(graph = G, type = "wtc", HM = "LV")
#' gplot(clv$gHM, l = "fdp")
#' table(clv$membership)
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
		mst <- mst(ug, weights = NULL, algorithm = NULL)
		G <- distances(mst, v = V(mst), to = V(mst), mode = "all", weights = NA)           
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
	if (K == 0) return(message("WARNING: no communities with size >=", size, "."))
	
	if (HM == "UV") {
		gHC <- cplot(graph, membership = membership)[-1]
		ftm <- Vxx <- NULL
	for (i in 1:K) {
		d <- igraph::degree(gHC[[i]], mode = "in")
		Vx <- V(gHC[[i]])$name[d == 0]
		Vy <- V(gHC[[i]])$name[d != 0]
		ftm <- rbind(ftm, cbind(Vx, rep(paste0("UV",cnames[i]), length(Vx))))
		ftm <- rbind(ftm, cbind(rep(paste0("UV",cnames[i]), length(Vy)), Vy))
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

	} else if (HM == "CV") {
		ftm <- data.frame(from = names(membership),
		                  to = c(paste0("CY", membership)))
		gLM <- graph_from_data_frame(ftm, directed = TRUE)
		V(gLM)$LV <- 0
		V(gLM)$LV[(vcount(gLM) - K + 1):vcount(gLM)] <- 1
		V(gLM)$color <- ifelse(V(gLM)$LV == 1, "lightblue", "green")

	} else if (HM == "none") {
		return( membership )
	}

    if (verbose == TRUE) {
		gplot(gLM, l = "fdp")
	}

	return(list(gHM = gLM, membership = membership))
}

#' @title Module scoring
#'
#' @description Generate factor scores, principal component scores, or
#' projection scores of latent, composite, and unmeasured variable modules,
#' respectively, and fit them with an exogenous group effect.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#' @param HM Hidden model type. For each defined hidden module:
#' (i) if \code{HM = "LV"}, a latent variable (LV) will be defined as
#' common unknown cause acting on cluster nodes; (ii) if \code{HM = "CV"},
#' cluster nodes will be considered as regressors of a latent composite
#' variable (CV); (iii) if \code{HM = "UV"}, an unmeasured variable (UV)
#' model will be generated for each module, where source nodes (i.e.,
#' in-degree = 0) act as common regressors influencing the other nodes
#' via an unmeasured variable.
#' By default, HM is set to "LV" (i.e., the latent variable model).
#' @param size Minimum number of nodes per hidden module. By default, a
#' minimum number of 5 nodes is required.
#' @param type Graph clustering method. If \code{type = "tahc"}, network
#' modules are generated using the tree agglomerative hierarchical
#' clustering method (Yu et al., 2015).
#' Other non-tree clustering methods from igraph package include: "wtc"
#' (default value; walktrap community structure with short random walks),
#' "ebc" (edge betweenness clustering), "fgc" (fast greedy method), "lbc"
#' (label propagation method), "lec" (leading eigenvector method), "loc"
#' (multi-level optimization), "opc" (optimal communiy structure), "sgc"
#' (spinglass statistical mechanics).
#' By default, the "wtc" method is used.
#' @param verbose A logical value. If TRUE, intermediate graphs will be
#' displayed during the execution. In addition, a reduced graph with
#' clusters as nodes will be fitted and showed to screen (see also
#' \code{\link[SEMgraph]{mergeNodes}}). By default, \code{verbode = FALSE}.
#' @param ... Currently ignored.
#'
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
#' Palluzzi F, Grassi M (2021). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' <arXiv:2103.08332>
#'
#' @return A list of 3 objects:
#' \enumerate{
#' \item "fit", hidden module fitting as a lavaan object;
#' \item "membership", hidden module nodes membership;
#' \code{\link[SEMgraph]{clusterGraph}} function;
#' \item "dataHM", data matrix with cluster scores in first columns.
#' }
#'
#' @examples
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' C <- clusterScore(graph = alsData$graph, data = als.npn,
#'                   group = alsData$group,
#'                   HM = "LV",
#'                   type = "wtc",
#'                   verbose = FALSE)
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
			return(message("UV is not applicable with udirected graph !"))
		}
		LXY <- clusterGraph(graph = ig, type = type,
		                HM = "UV", size = size,
		                verbose = verbose)
		if(length(LXY) == 0) return(list(fit = NA, M = NA, dataHM = NA))
		membership <- LXY[[2]]
		gLC <- cplot(graph, membership = membership)[-1]
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
	dataHM<- cbind(dataLC, dataY)

	return(list(fit = fsr, membership = membership, dataHM = dataHM))
}

#' @title Cluster extraction utility
#'
#' @description Extract and fit clusters from an input graph.
#'
#' @param graph Input network as an igraph object.
#' @param membership A vector of cluster membership IDs. If NULL, clusters
#' will be automatically generated with \code{\link[SEMgraph]{clusterGraph}}
#' using the edge betweenness clustering ("ebc") algorithm.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. Group specification enables node perturbation
#' testing. By default, \code{group = NULL}.
#' @param map Logical value. If TRUE, the plot of the input graph
#' (coloured by cluster membership) will be generated along with independent
#' module plots. If the input graph is very large, plotting could be
#' computationally intensive (by default, \code{map = FALSE}).
#' @param verbose Logical value. If TRUE, a plot will be showed for each
#' cluster.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return List of clusters as igraph objects and fitting results for
#' each cluster as a lavaan object.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @examples
#'
#' \donttest{
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' adjdata <- SEMbap(alsData$graph, als.npn)$data
#'
#' # Clusters creation
#' clusters <- extractClusters(graph = alsData$graph, data = adjdata)
#' head(parameterEstimates(clusters$fit$HM1))
#' head(parameterEstimates(clusters$fit$HM2))
#' head(parameterEstimates(clusters$fit$HM4))
#' gplot(clusters$clusters$HM2)
#'
#' # Map cluster on the input graph
#' g <- alsData$graph
#' c <- clusters$clusters$HM2
#' V(g)$color <- ifelse(V(g)$name %in% V(c)$name, "gold", "white")
#' gplot(g)
#'
#' }
#'
extractClusters<- function(graph, data, group = NULL, membership = NULL, map = FALSE, verbose = FALSE, ...) 
{
	if (is.null(membership)) {
	 membership<- clusterGraph(graph, type="ebc", HM="none", size=5, verbose=FALSE)
	}
	clusters <- cplot(graph, membership, l=layout.auto, map, verbose)[-1]
	N <- length(clusters)
	if ("HM9999" %in% names(clusters)) N <- N - 1
	res <- NULL
	lav <- list()
	for (i in 1:N) {
		cat("\r","cluster=", i, "of", N)
		flush.console()
		fit<- quiet(SEMrun(clusters[[i]], data, group))
		if(is.null(fit)) next
		if (!is.null(group) & vcount(clusters[[i]]) > 100){ 
		 dev_df <- fit$fit$ricf$dev/fit$fit$ricf$df
		 srmr <- fit$fit$fitIdx[3]
		 pv1<- Brown.test(x=fit$dataXY[,-1], p=fit$gest$pvalue, theta=fit$gest$Stat, tail="positive")
		 pv2<- Brown.test(x=fit$dataXY[,-1], p=fit$gest$pvalue, theta=fit$gest$Stat, tail="negative")
		}else{
		 dev_df <- fitMeasures(fit$fit, "chisq")/fitMeasures(fit$fit, "df")
		 srmr <- fitMeasures(fit$fit, "srmr")
		 pv1<- Brown.test(x=fit$dataXY[,-1], p=fit$gest$pvalue, theta=fit$gest$est, tail="positive")
		 pv2<- Brown.test(x=fit$dataXY[,-1], p=fit$gest$pvalue, theta=fit$gest$est, tail="negative")
		}
		
		dfc<- data.frame(
		 cluster = names(clusters)[i],
		 n.nodes = vcount(clusters[[i]]),
		 n.edges = ecount(clusters[[i]]),
		 dev_df = round(dev_df, 3),
		 srmr = round(srmr, 3),
		 V.pv.act = round(pv1, 6),
		 V.pv.inh = round(pv2, 6))
		 
		res <- rbind(res, dfc)
		lav <- c(lav, list(fit$fit))
	}
	message("\n\nFound ", nrow(res), " clusters with > 5 nodes")
	rownames(res) <- NULL
	names(lav) <- res$cluster
	
	return(list(clusters = clusters, fit = lav, dfc = res))
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
#'
#' \donttest{
#'
#' # Clustering ALS graph with WTC method
#' G <- alsData$graph
#' membership <- clusterGraph(graph = G, type = "wtc")
#' cplot(G, membership, map = TRUE, verbose = FALSE)
#' cplot(G, membership, map = FALSE, verbose = TRUE)
#' # The list of cluster graphs !
#' cg <- cplot(G, membership); cg
#'
#' }
#'
cplot<- function (graph, membership, l = layout.auto, map = FALSE, verbose = FALSE, ...) 
{
    V(graph)$M <- 9999
    V(graph)$M[which(V(graph)$name %in% names(membership))] <- membership
    V(graph)$color <- V(graph)$M + 1
	V(graph)$weight <- 1
    if (map) {
		verbose <- FALSE
	    plot(graph, layout = l)
    }
    M <- names(table(V(graph)$M))
    K <- length(table(V(graph)$M))
    vcol <- as.numeric(M) + 1
    HM <- lapply(1:K, function(x)
	       induced_subgraph(graph, names(membership[membership == M[x]])))
    if ("9999" %in% M) {
     HM[[K]] <- induced_subgraph(graph, V(graph)$name[V(graph)$M == 9999])
	}
	names(HM) <- paste0("HM", M)
    d <- igraph::degree(graph, mode = "all") * 2 + 1
    if (verbose) {
        glv <- lapply(1:K, function(x) {
                plot(HM[[x]], vertex.color = vcol[x], vertex.size = d[V(HM[[x]])$name], 
                 layout = l, main = paste0("Hidden Module ", 
                  M[x]))
            Sys.sleep(3)
        })
    }
    return(invisible(c(list(graph = graph), HM)))
}

#' @title Graph nodes merging by a membership attribute
#'
#' @description Merge groups of graph nodes using hierarchical clustering
#' with prototypes derived from \code{\link[protoclust]{protoclust}} or 
#' custom membership attribute (e.g., cluster membership derived from
#' \code{\link[SEMgraph]{clusterGraph}}).
#' @param graph network as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param h Cutting the minimax clustering at height, h = 1 - abs(cor(j,k)),
#' yielding a merged node (and a reduced data set) in which every node in the
#' cluster has correlation of at least cor(j,k) with the prototype node.
#' By default, \code{h = 0.5}, i.e. cor(j,k) = 0.5.
#' @param membership Cluster membership. A vector of cluster membership
#' identifiers as numeric values, where vector names correspond to graph
#' node names. By default, \code{membership = NULL}.
#' @param HM Hidden cluster label. If membership is derived from clusterGraph:
#' HM = "LV", a latent variable (LV) will be defined as common unknown cause
#' acting on cluster nodes. If HM = "CV", cluster nodes will be considered as regressors of a
#' latent composite variable (CV). Finally, if HM = "UV", an unmeasured
#' variable (UV) is defined, where source nodes of the module (i.e.,
#' in-degree = 0) act as common regressors influencing the other nodes
#' via an unmeasured variable. By default, \code{HM = NULL}
#' @param ... Currently ignored.
#'
#' @details Hierarchical clustering with prototypes (or Minmax linkage) is
#' unique in naturally associating a node (the prototypes) with every
#' interior node of the dendogram. Thus, for each merge we have a single
#' representative data point for the resulting cluster (Bien, Tibshirani, 2011).
#' These prototypes can be used to greatly enhance the interpretability of
#' merging nodes and data reduction for SEM fitting.
#
#' @export
#'
#' @return A list of 2 objects is returned:
#' \enumerate{
#' \item "gLM", A graph with merged nodes as an igraph object;
#' \item "membership", cluster membership vector for each node.
#' }
#'
#' @seealso \code{\link[SEMgraph]{clusterGraph}}
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Bien J, Tibshirani R (2011). Hierarchical Clustering With Prototypes via
#' Minimax Linkage. Journal of the American Statistical Association
#' 106(495): 1075-1084. <doi:10.1198/jasa.2011.tm10183>
#'
#' @examples
#'
#' # Gene memberships with prototypes with h=0.5
#' G <- properties(alsData$graph)[[1]]
#' M <- mergeNodes(G, alsData$exprs, h = 0.5)
#'
#' # Gene memberships with EBC method and size=10
#' m <- clusterGraph(G, type = "ebc", size = 10)
#' M <- mergeNodes(G, alsData$exprs, membership = m, HM = "LV")
#'
#' # Gene memberships defined by user
#' c1 <- c("5894", "5576", "5567", "572", "598")
#' c2 <- c("6788", "84152", "2915", "836", "5530")
#' c3 <- c("5603", "6300", "1432", "5600")
#' m <- c(rep(1,5), rep(2,5), rep(3,4))
#' names(m) <- c(c1, c2, c3)
#' M <- mergeNodes(G, alsData$exprs, membership = m, HM = "CV")
#'
mergeNodes<- function(graph, data, h = 0.5, membership = NULL, HM = NULL, ...)
{
	# Set graph and membership object :
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	graph <- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	if(is.numeric(membership)){
	 nodes <- names(membership)
	 membership <- paste0(HM, membership)
	 names(membership) <- nodes
	}else{
	 membership <- prototype(data=data[,nodes], h=h, size=3)
	}

	LM<- NULL
	for ( i in 1:length(table(membership)) ) {
	 m<- names(table(membership))[i]
	 LMi<- V(graph)$name[which(V(graph)$name %in% names(membership)[membership==m])]
	 LM<- c(LM, list(LMi))
	}
	names(LM)<- names(table(membership))
		
	# visualize graph object :
	gLM<- as_graphnel(graph)
	for ( i in 1:length(LM) ) {
	 gLMi<- graph::combineNodes(LM[[i]], gLM, names(LM)[i], mean)
	 gLM<- gLMi
	}
	
	ig<- graph_from_graphnel(gLM)
	if( length(V(ig)$color) == 0 ) V(ig)$color<- "white"
	if( length(HM) == 1){
	  V(ig)$color[substr(V(ig)$name,2,2)=="V"]<- "pink"
	} else{
	  V(ig)$color[substr(V(ig)$name,1,1)=="p"]<- "green"
	}
	vcol<- V(ig)$color
	names(vcol)<- V(ig)$name
	V(ig)$name <- gsub("p", "", V(ig)$name)
	gplot(ig)
	
	return(list(gLM=ig, membership=membership))
}

prototype<- function(data, h, size, ...)
{
	# hierarchical clustering with prototypes
	D <- as.dist(1-abs(cor(data)))
	hc <- protoclust::protoclust(D)
	plot(hc);abline(h=h, lty=1, col="red")
	cutd <- protoclust::protocut(hc, k=NULL, h=h)
	protos <- hc$labels[cutd$protos]
	cln <- sort(cutd$cl)
	
	# cluster membership with nodes > size
	nrep <- as.numeric(table(cln))
	cl <- unlist(lapply(1:length(protos), function(x) rep(paste0("p",protos[x]), nrep[x])))
	names(cl) <- names(cln)
	csize <- names(table(cl))[table(cl) >= size]
	cl <- cl[cl %in% csize]
		
	return(membership = cl)
}
