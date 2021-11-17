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
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #

#' @title SEM-based gene set analysis
#'
#' @description Gene Set Analysis (GSA) via self-contained test for group
#' effect on signaling (directed) pathways as SEM, evaluating overall pathway
#' perturbation, perturbation emission from source nodes, and perturbation
#' accumulation on target nodes. Approximate randomization test P-values
#' of specific node and aggregated group effects will be computed.
#' For directed graphs, they include: the sum of group effects adjusted
#' by residual variances (D), the sum of the tagret nodes perturbation
#' (i.e., group effect) accumulation from source nodes (A), and the sum
#' of the source nodes perturbation emission towards target nodes (E).
#' For undirected graphs, the sum of group effects (D), adjusted
#' by residual variances, will be estimated.
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
#' @return A list of 2 objects:
#' \enumerate{
#' \item "gsa", A data.frame reporting the following information for each
#' pathway in the input list:
#' \itemize{
#' \item "N.nodes", pathway size (number of nodes);
#' \item "N.DRNs", number of differential regulated genes (DRNs) within the pathway,
#' after multiple test correction with Benjamini-Hochberg method;
#' \item "pD", significance of the sum of group effects, adjusted by the
#' residual variance;
#' \item "pA", significance of the sum of tagret nodes perturbation
#' (i.e., group effect) accumulation from source nodes;
#' \item "pE", significance of the sum of source nodes perturbation
#' (i.e., group effect) emission towards target nodes;
#' \item "pvalue", Fisher's combined P-value of pD, pA, and pE.
#' }
#' \item "DRN", a list with DRNs names per pathways.
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats pchisq p.adjust na.omit
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' \dontrun{
#'
#' # Nonparanormal(npn) transformation

#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' # Selection of FTD pathways from kegg.pathways.Rdata
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
#' GSA$DRN
#' 
#' }
#'
SEMgsa <- function(g = list(), data, group, method = "BH", alpha = 0.05,
                   n_rep = 1000, ...)
{
	# Set SEM objects
	pX2 <- function(x) 1 - pchisq(-2*sum(log(x)), 2*length(x))
	gs <- names(g)
	K <- length(g)
	del <- NULL
	res.tbl <- NULL
	DRN <- list()

	for (k in 1:K) {
		cat(paste("# ", k, ": ", gs[k], "\n", sep = ""))
		quiet(ig <- properties(g[[k]])[[1]])

		# SEM fitting
		fit <- NULL
		err <- paste(" ValueError: none of pathway #", k,
		             " variables are present in the dataset.\n",
		             sep = "")
		tryCatch(quiet(fit <- SEMricf(graph = ig, data, group,
		                              random.x = FALSE,
	                                  n_rep)),
	             error = function(c) cat(err))

		if (length(fit[[1]]) == 0) {
			del <- c(del, k)
			next
		}

		p <- ncol(fit$dataXY)
		B <- (diag(p) - fit$fit$ricf$Bhat)[-1, -1]
		if (sum(B) == 0) {
			del <- c(del, k)
			next
		}

		pval <- fit$gest$pvalue[-c(1:3)]
		genes <- gsub("X", "", rownames(fit$gest))[-c(1:3)]
		genes <- genes[p.adjust(pval, method = method) < alpha]
		DRN <- c(DRN, list(genes))

		# data.frame of SEM results
		df <- data.frame(
		 N.nodes = vcount(ig),
		 #N.edges= ecount(ig),
		 N.DRNs = sum(p.adjust(pval, method = method) < alpha),
		 pD = round(fit$gest[1, 4], 6),
		 pA = round(fit$gest[2, 4], 6),
		 pE = round(fit$gest[3, 4], 6))
		pvalue = round(pX2(x = df[1, 3:5]), 6)
		res.tbl <- rbind(res.tbl, cbind(df, pvalue))
	}

	if(is.null(del)) {
		rownames(res.tbl) <- gs
		names(DRN) <- gs
	} else {
		rownames(res.tbl)<- gs[-del]
		names(DRN)<- gs[-del]
	}

	return(list(gsa = res.tbl, DRN = DRN))
}

#' @title Graph properties summary and graph decomposition
#'
#' @description Produces a summary of network properties and returns
#' graph components (ordered by decreasing size), without self-loops.
#' @param graph Input network as an igraph object.
#' @param data An optional data matrix whith rows corresponding to subjects,
#' and columns to graph nodes (variables). Nodes will be mapped onto
#' variable names.
#' @param ... Currently ignored.
#'
#' @import igraph
#' @export
#'
#' @return List of graph components, ordered by decreasing size (the first
#' component is the giant one), without self-loops.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Extract the "Type II diabetes mellitus" pathway:
#' g <- kegg.pathways[["Type II diabetes mellitus"]]
#' properties(g)
#' summary(g)
#'
properties <- function(graph, data = NULL, ...)
{
    if (!is_igraph(graph)) ig <- graph_from_graphnel(graph) else ig <- graph
    if (!is.null(data)) {
		nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
		ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	}
	ig <- simplify(ig, remove.loops = TRUE)
	gcs <- igraph::decompose.graph(ig)
	vsize <- lapply(1:length(gcs), function(x) vcount(gcs[[x]]))
	ig1 <- gcs[[which.max(vsize)]]
	cat("Frequency distribution of graph components\n\n")
	tab <- table(sapply(gcs, vcount))
	tab <- data.frame(n.nodes = as.numeric(names(tab)),
	                  n.graphs = as.numeric(tab))
	print(tab)
	cat("\nPercent of vertices in the giant component:",
	    round(100*vcount(ig1)/vcount(ig), 1), "%\n\n")

	print(c(is.simple = is_simple(ig1),
	        #is.connected = is_connected(ig1),
	        is.dag = is_dag(ig1),
	        is.directed = is_directed(ig1),
	        is.weighted = is_weighted(ig1)))

	cat("\n")
	print(c(which.mutual = table(which_mutual(ig1))))

	return(invisible(list(ig1 = ig1, gcs = gcs[-which.max(vsize)])))
}

mergeGraph <- function(g = list(), gref = NULL, gnet = NULL, verbose = FALSE,
                       ...)
{
	g1 <- list()
	for (i in 1:length(g)) {
		V(g[[i]])$color <- NA
		g[[i]] <- igraph::remove.vertex.attribute(g[[i]], "color")
		if (!is_igraph(gnet) || !is_directed(gnet)) {
			if (is.directed(g[[i]])) {
				g1 <- c(g1, list(as.undirected(g[[i]], mode = "collapse")))
			} else {
				g1 <- c(g1, list(g[[i]]))
			}
		} else {
			if (!is.directed(g[[i]])) {
				g1 <- c(g1, list(orientEdges(graph = g[[i]], dg = gnet)))
			} else {
				g1 <- c(g1, list(g[[i]]))
			}
		}
	}

	Ug <- igraph::graph.union(g1)

	if (is_igraph(gref)) {
		V(Ug)$color <- ifelse(V(Ug)$name %in% V(gref)$name,
		                      "lightblue",
		                      "yellow")
	} else {
		V(Ug)$color <- "white"
	}
	if (verbose) plot(Ug)

	return(Ug)
}

#' @title Graph plotting with renderGraph
#'
#' @description Wrapper for function renderGraph of the R package
#' Rgraphwiz.
#'
#' @param graph An igraph or graphNEL object.
#' @param l Any layout supported by \code{Rgraphviz}. It can be one among:
#' "dot" (default), "neato", "circo", "fdp", "osage", "twopi".
#' @param main Plot main title (by default, no title is added).
#' @param cex.main Main title size (default = 1).
#' @param font.main Main title font (default = 1). Available options
#' are: 1 for plain text, 2 for bold, 3 for italics, 4 for bold italics,
#' and 5 for symbol.
#' @param color.txt Node text color (default = "black").
#' @param fontsize Node text size (default = 16).
#' @param cex Another argument to control node text size (default = 0.6).
#' @param shape Node shape (default = "circle").
#' @param color Node border color (default = "gray70").
#' @param lty Node border outline (default = 1).
#' Available options include: 0 for blank, 1 for solid line, 2 for dashed,
#' 3 for dotted, 4 for dotdash, 5 for longdash, and 6 for twodash.
#' @param lwd Node border thickness (default = 1).
#' @param h Manual node height (default = "auto").
#' @param w Manual node width (default = "auto").
#' @param psize Automatic node size (default = 80).
#' @param ... Currently ignored.
#'
#' @import igraph
#' @importFrom graph nodes edgeNames isDirected nodeRenderInfo edgeRenderInfo graphRenderInfo
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @export
#'
#' @return gplot returns invisibly the graph object produced by Rgraphviz
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' gplot(sachs$graph, main = "input graph")
#'
#' sem <- SEMrun(sachs$graph, sachs$pkc)
#' gplot(sem$graph, main = "output graph")
#'
gplot <- function(graph, l = "dot", main = "", cex.main = 1, font.main = 1,
                  color.txt = "black", fontsize = 16, cex = 0.6,
                  shape = "circle", color = "gray70", lty = 1, lwd = 1,
                  w = "auto", h = "auto", psize = 80, ...)
{
	# Set graphNEL object
	g <- as_graphnel(graph)

	vcol <- V(graph)$color
	vshape <- V(graph)$shape
	vsize <- V(graph)$size
	vlab <- V(graph)$label
	ecol <- E(graph)$color
	elwd <- E(graph)$width
	elab <- E(graph)$label
	
	if (length(vcol) > 0) {
		names(vcol) <- V(graph)$name
		vcol <- vcol[graph::nodes(g)]
	}
	if (length(vshape) > 0) {
		names(vshape) <- V(graph)$name
		vshape <- vshape[graph::nodes(g)]
	}
	if (length(vsize) > 0) {
		names(vsize) <- V(graph)$name
		vsize <- 10*vsize[graph::nodes(g)]
	} else {
		vsize <- rep(psize, length = vcount(graph))
		names(vsize) <- V(graph)$name
	}
	if (length(vlab) > 0) {
		names(vlab) <- V(graph)$name
		vlab <- vlab[graph::nodes(g)]
	}
	if (length(color) > 1) {
		names(color) <- V(graph)$name
		color <- color[graph::nodes(g)]
	}
	if (length(ecol) > 0) {
		names(ecol) <- gsub("\\|", "~", attr(E(graph), "vnames"))
		ecol <- ecol[graph::edgeNames(g, recipEdges = "distinct")]
	}
	if (length(elwd) > 0) {
		names(elwd) <- gsub("\\|", "~", attr(E(graph), "vnames"))
		elwd <- elwd[graph::edgeNames(g, recipEdges = "distinct")]
	}
	if (length(elab) > 0) {
		names(elab) <- gsub("\\|", "~", attr(E(graph), "vnames"))
		elab <- elab[graph::edgeNames(g,recipEdges = "distinct")]
	} else {
		elab <- rep("", ecount(graph))
		names(elab) <- gsub("\\|", "~", attr(E(graph), "vnames"))
	}

	if (graph::isDirected(g)) {
		g <- Rgraphviz::layoutGraph(g, layoutType = l,
		                            edgeAttrs = list(label = elab))
	} else {
		g <- Rgraphviz::layoutGraph(g, layoutType = "fdp",
		                            edgeAttrs = list(label = elab))
	}

	if (w == "auto") w <- vsize
	if (h == "auto") h <- vsize

	graph::nodeRenderInfo(g)<- list(col = color, fill = vcol, lty = lty,
	                                label = vlab, lwd = lwd,
	                                textCol = color.txt,
	                                fontsize = fontsize, cex = cex,
	                                shape = vshape, width = w, height = h)
	graph::edgeRenderInfo(g) <- list(col = ecol, lty = 1, lwd = elwd)
	graph::graphRenderInfo(g) <- list(main = main, cex.main = cex.main,
	                                  font.main = font.main)
	Rgraphviz::renderGraph(g)
	
	return(invisible(g))
}

#' @title Correlation matrix to graph
#'
#' @description Convert a correlation matrix to an igraph object.
#' @param R Correlation matrix.
#' @param n Sample size (i.e., the number of subjects).
#' @param alpha Significance level used to compute the correlation threshold.
#' By default, \code{alpha = 0.05}.
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default,
#' \code{method = "none"} (i.e., no multiple test correction).
#' See \code{\link[stats]{p.adjust}} for other correction methods.
#' @param type Graph building method. If \code{type} is either
#' \code{"marg"} or \code{"cond"}, marginal or conditional correlation
#' tests will be used, respectively.
#' If \code{type = "mst"}, input correlations are converted to distances
#' and a minimum spanning tree is generated from the distance matrix,
#' using Prim's algorithm (Prim, 1957).
#' If \code{type = "tmfg"}, a triangulate maximally graph is generated
#' from the given correlation matrix (Massara et al., 2016).
#' @param ... Currently ignored.
#'
#' @import lavaan
#' @import igraph
#' @export
#'
#' @references
#'
#' Palluzzi F, Grassi M (2021). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' <arXiv:2103.08332>
#'
#' Massara GP, Di Matteo T and Aste T (2009). Network Filtering for Big
#' Data: Triangulated Maximally Filtered Graph.
#' Journal of complex Networks, 5(2): 161--178.
#' <https://doi.org/10.1093/comnet/cnw015>
#'
#' Prim RC (1957). Shortest connection networks and some generalizations.
#' Bell System Technical Journal, 36(6):1389--1401.
#' <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
#'
#' @return An igraph object.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Graphs creation
#' C1 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
#'                  type = "marg",
#'                  method = "BH")
#' C2 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
#'                  type = "cond",
#'                  method = "BH")
#' C3 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
#'                  type = "mst",
#'                  method = "BH")
#' C4 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
#'                  type = "tmfg",
#'                  method = "BH")
#'
#' # Graphs plots
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar= rep(2, 4))
#' plot(C1, layout=layout.circle, main= "marg"); box(col="gray")
#' plot(C2, layout=layout.circle, main= "cond"); box(col="gray")
#' plot(C3, layout=layout.circle, main= "mst"); box(col="gray")
#' plot(C4, layout=layout.circle, main= "tmfg"); box(col="gray")
#' par(old.par)
#' 
corr2graph <- function(R, n, type = "marg", method = "none",
                       alpha = 0.05, ...)
{
	# Set correlation matrix
	p <- nrow(R)
	q <- 0
	if (type == "cond") {
		q <- p - 2
		if (corpcor::is.positive.definite(R)) {
			K <- corpcor::cor2pcor(R)
			rownames(K) <- colnames(K) <- rownames(R)
		} else {
			K <- corpcor::pcor.shrink(R, verbose = TRUE)[1:p, 1:p]
		}
	} else {
		K <- R
	}

	if (type == "marg" | type == "cond") {
		# select the correlation threshold
		z <- abs(atanh(K[lower.tri(K)]))/sqrt(1/(n - 3 - q))
		p.adj <- p.adjust(2*(1 - pnorm(z)), method = method)
		Z <- min(z[which(p.adj < alpha)])
		#Z <- qnorm(alpha/2, lower.tail = FALSE)
		thr <- (exp(2*Z/sqrt(n - 3 - q)) - 1)/(exp(2*Z/sqrt(n - 3 - q)) + 1)
		A0 <- ifelse(abs(K) > thr, 1, 0)
		diag(A0) <- 0
		del <- which(colSums(A0) == 0)
		if (length(del) > 0) A <- A0[-del, -del] else A <- A0
		ug <- graph_from_adjacency_matrix(A, mode = "undirected")
	}
    if (type == "mst") {
		D <- diag(p) - K^2
		gA <- graph_from_adjacency_matrix(D, mode = "undirected",
		                                  weighted = TRUE)
		ug <- igraph::mst(gA, algorithm = "prim")
	}
	if (type == "tmfg") ug <- TMFG(K)$graph

	return(graph = ug)
}

#' @title lavaan model to graph
#'
#' @description Convert a model, specified using lavaan syntax,
#' to a graph object in either igraph or dagitty format.
#' @param model Model specified using lavaan syntax.
#' @param directed Logical value. If TRUE (default), edge directions from
#' the model will be preserved. If FALSE, the resulting graph will
#' be undirected.
#' @param psi Logical value. If TRUE (default) covariances will be converted
#' into bidirected graph edges. If FALSE, covariances will be excluded from
#' the output graph.
#' @param format Output graph format. It can be either "igraph" (default)
#' or "dagitty".
#' @param verbose Logical value. If TRUE, a plot of the output graph will
#' be generated. For large graphs, this could significantly increase
#' computation time. If FALSE (default), graph plotting will be disabled.
#' @param ... Currently ignored.
#'
#' @import lavaan
#' @import igraph
#' @importFrom dagitty graphLayout canonicalize
#' @export
#'
#' @return An igraph object.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Writing path diagram in lavaan syntax
#'
#' model<-'
#' #path model
#' Jnk ~ PKA + PKC
#' P38 ~ PKA + PKC
#' Akt ~ PKA + PIP3
#' Erk ~ PKA + Mek
#' Mek ~ PKA + PKC + Raf
#' Raf ~ PKA + PKC
#' PKC ~ PIP2 + Plcg
#' PIP2 ~ PIP3 + Plcg
#' Plcg ~ PIP3
#' #PKA ~ 1
#' #PIP3 ~ 1
#'
#' # (co)variances
#' # PIP2 ~~ PIP3
#' '
#'
#' # Graph with covariances
#' G0 <- lavaan2graph(model, psi = TRUE)
#' plot(G0, layout = layout.circle)
#'
#' # Graph without covariances
#' G1 <- lavaan2graph(model, psi = FALSE)
#' plot(G1, layout = layout.circle)
#'
lavaan2graph <- function(model, directed = TRUE, psi = TRUE,
                         format = "igraph", verbose = FALSE, ...)
{
	lav <- lavParTable(model, fixed.x = FALSE)
	lavb <- subset(lav, lav$op == "~")
	lavc <- subset(lav, lav$op == "~~" & (lav$rhs != lav$lhs))
	ftm <- data.frame(cbind(from = lavb$rhs, to = lavb$lhs, label = lavb$label),
	                  color="blue")
	if (nrow(lavc) != 0 & psi == TRUE) {
		ftmc1 <- data.frame(cbind(from = lavc$rhs, to = lavc$lhs, label = "",
		                    color = "gray60"))
		ftmc2 <- data.frame(cbind(from = lavc$lhs, to = lavc$rhs, label = "",
		                    color = "gray60"))
		ftm <- rbind(ftm, ftmc1, ftmc2)
	}
	graph <- graph_from_data_frame(ftm, directed = directed)
	if (format == "dagitty") graph <- graph2dagitty(graph, verbose)
	if (format == "igraph" & verbose) plot(graph)
	return(graph)
}

#' @title Graph to lavaan model
#'
#' @description Convert an igraph object to a model (lavaan syntax).
#' @param graph A graph as an igraph object.
#' @param nodes Subset of nodes to be included in the model. By default,
#' all the input graph nodes will be included in the output model.
#' @param ... Currently ignored.
#'
#' @import lavaan
#' @import igraph
#' @importFrom dagitty graphLayout canonicalize
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Graph (igraph object) to structural model in lavaan syntax
#' model <- graph2lavaan(sachs$graph)
#' cat(model, "\n")
#'
#' @return A model in lavaan syntax.
#'
graph2lavaan <- function(graph, nodes = V(graph)$name, ...)
{
	# Set from-to-matrix representation of edge links
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	ftm <- as_data_frame(ig)

	if (is.directed(ig) & sum(which_mutual(ig)) > 0) {
		sel <- as.numeric(c(E(ig)[which_mutual(ig)]))
		ftm <- as_data_frame(ig)[-sel,]
		ubg <- as.undirected(graph_from_data_frame(as_data_frame(ig)[sel,]))
		ftb <- as_data_frame(ubg)
	} else {
		ftb <- NULL
	}

	modelY <- modelV <- vector()
	if (is.directed(ig)) {
		for(j in 1:nrow(ftm)) {
			modelY[j] <- paste0(ftm[j, 2], "~", ftm[j, 1])
		}
		if (length(ftb) > 0) {
			for(k in 1:nrow(ftb)) {
				modelV[k] <- paste0(ftb[k, 2], "~~", ftb[k, 1])
			}
		}
	} else {
		for(j in 1:nrow(ftm)) modelY[j] <- paste0(ftm[j, 2], "~~", ftm[j, 1])
	}
	model <- paste(c(sort(modelY), modelV), collapse = "\n")
	return(model)
}

#' @title Graph conversion from igraph to dagitty
#'
#' @description Convert an igraph object to a dagitty object.
#' @param graph A graph as an igraph or dagitty object.
#' @param canonical A logical value. If TRUE, DAG conversion is enforced
#' (for \code{graph2dagitty} only). This argument is FALSE by default.
#' @param verbose A logical value. If TRUE, the output graph is shown
#' (for \code{graph2dagitty} only). This argument is FALSE by default.
#' @param ... Currently ignored.
#'
#' @import lavaan
#' @import igraph
#' @importFrom dagitty graphLayout canonicalize
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Graph as an igraph object to dagitty object
#' G <- graph2dagitty(sachs$graph)
#' plot(dagitty::graphLayout(G))
#'
#' @return A dagitty object.
#'
graph2dagitty <- function (graph, canonical = FALSE, verbose = FALSE, ...)
{
    dg <- graph - E(graph)[which_mutual(graph)]
	ug <- as.undirected(graph - E(graph)[!which_mutual(graph)])
	ed <- attr(E(dg), "vnames")
    eb <- attr(E(ug), "vnames")
    de <- paste(gsub("\\|", "->", ed), collapse = "\n")

    if (length(eb) == 0) {
		dagi <- paste0("dag {\n", de, "\n}")
    } else {
		be <- paste(gsub("\\|", "<->", eb), collapse = "\n")
		dagi <- paste0("dag {\n", de, "\n", be, "\n}")
    }

	if (verbose) plot(dagitty::graphLayout(dagi))
	if (canonical) {
		dagi <- dagitty::canonicalize(dagi)
		if (verbose) plot(dagitty::graphLayout(dagi$g))
		return(dagi$g)
	}
    return(dagi)
}

#' @title Convert directed graphs to directed acyclic graphs (DAGs)
#'
#' @description Remove cycles and bidirected edges from a directed graph.
#'
#' @param graph A directed graph as an igraph object.
#' @param data A data matrix with subjects as rows and variables as
#' columns.
#' @param bap If TRUE, a bow-free acyclic path (BAP) is returned
#' (default = FALSE).
#' @param time.limit CPU time for the computation, in seconds
#' (default = Inf).
#' @param ... Currently ignored.
#'
#' @details The conversion is performed firstly by removing bidirected
#' edges and then the data matrix is used to compute edge P-values, through
#' marginal correlation testing (see \code{\link[SEMgraph]{weightGraph}},
#' r-to-z method). When a cycle is detected, the edge with highest
#' P-value is removed, breaking the cycle. If the bap argument is TRUE,
#' a BAP is generated merging the output DAG and the bidirected edges
#' from the input graph.
#'
#' @import igraph
#' @export
#'
#' @return A DAG as an igraph object.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' dag <- graph2dag(graph = sachs$graph, data = log(sachs$pkc))
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2), mar=rep(1, 4))
#' gplot(sachs$graph, main = "Input graph")
#' gplot(dag, main = "Output DAG")
#' par(old.par)
#'
graph2dag <- function(graph, data, bap = FALSE, time.limit = Inf, ...)
{
	if (is_dag(graph)) return(dag = graph)

	# Graph weighting by edge pvalues (r2z)
	graph <- weightGraph(graph, data, group = NULL, seed = "none")
	E(graph)$weight <- 1/(-log(E(graph)$pv))
	ftm <- as_data_frame(graph)
	wE <- ftm$weight
	names(wE) <- paste0(ftm[, 1], ":", ftm[, 2])

	# Delete all mutual edges <-> , i.e. <- & ->
	ig <- graph - E(graph)[which_mutual(graph)]
	if (is_dag(ig)) {
		cat("DAG conversion: TRUE\n")
		if (bap == TRUE) return(bap = graph)
		return(dag = ig)
	}

	# Subgraph isomorphism algorithm to detect all cycles of a given length
	find.cycles <- function(graph, k, time.limit) {
		ring <- graph.ring(k, TRUE)
		subgraph_isomorphisms(ring, graph, "lad", time.limit = time.limit)
	}

	# Function that identifies the right subisomorphisms to keep
	subisomorphism_rm_permutation <- function(si) {
		is_first_min <- function(x) {
			return(x[1] == min(x))
		}
		sel <- lapply(si, is_first_min)
		return(si[unlist(sel)])
	}

	# Function that search max(weight) edge on each cycle
	max_edges <- function(x) {
		Ec <- vector()
		for (i in 1:(length(x)-1)) Ec <- c(Ec, paste0(x[i], ":", x[i + 1]))
		Ew <- wE[which(names(wE) %in% Ec)]
		if(length(Ew) == 0) Ew <- wE[1]
		return(unlist(strsplit(names(Ew)[which(Ew == max(Ew))], ":")))
	}

	for (k in 3:vcount(graph)) {
		# Find all cycles with k vertices(edges)
		l <- find.cycles(ig, k, time.limit)
		# Remove permutations
		l <- subisomorphism_rm_permutation(si = l)
		# Extract the vertices in each cycle
		if (length(l) == 0) next
		cycles <- lapply(1:length(l), function(x) names(l[[x]]))
		# Edges with max(weight)
		l <- unique(lapply(cycles, max_edges))
		E1 <- unlist(lapply(l, function(x) paste0(x[1], "|", x[2])))
		E0 <- attr(E(ig), "vnames")
		ig <- delete_edges(ig, which(E0 %in% E1))
	}
	cat("DAG conversion :", is_dag(ig), "\n")

	# Add all mutual edges <-> , i.e. <- & ->
	if (bap == TRUE) {
		e <- as_edgelist(graph - E(graph)[!which_mutual(graph)])
		return(bap = add_edges(ig, as.vector(t(e))))
	}

	return(graph = ig)
}

#' @title Assign edge orientation of an undirected graph
#'
#' @description Assign edge orientation of an undirected graph
#' through a given reference directed graph.
#'
#' @param ug An undirected graph as an igraph object.
#' @param dg A directed reference graph.
#' @param ... Currently ignored.
#'
#' @import igraph
#' @export
#'
#' @return A directed graph as an igraph object.
#'
#' @examples
#'
#' # Graphs definition
#' G0 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc), type = "marg")
#'
#' # Reference graph-based orientation
#' G1 <- orientEdges(ug = G0, dg = sachs$graph)
#'
#' # Graphs plotting
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2), mar=rep(2,4))
#' plot(G0, layout=layout.circle, main = "Input undirected graph")
#' plot(G1, layout=layout.circle, main = "Output directed graph")
#' par(old.par)
#'
orientEdges<- function(ug, dg, ...)
{
	if (is_directed(ug)){
	 return(message("ERROR: the input graph is a Directed graph !"))
	}
	if (!is_directed(dg)){
	 return(message("ERROR: the reference graph is an Undirected graph !"))
	}
	mg <- as.directed(ug, mode = "mutual")
	exy0 <- attr(E(mg), "vnames")
	exy1 <- attr(E(dg)[which_mutual(dg) == FALSE], "vnames")
	exy2 <- exy0[which(exy0 %in% exy1)]
	if (length(exy2) == 0) return(graph = mg)
	str2 <- strsplit(exy2,"\\|")
	ftm2 <- matrix(unlist(str2),nrow=length(str2),byrow=TRUE)
	g2 <- graph_from_edgelist(ftm2, directed=TRUE)
	ug0 <- difference(ug, as.undirected(g2))
	mg0 <- as.directed(ug0, mode = "mutual")
	g <- graph.union(mg0, g2)
	E1 <- attr(E(g), "vnames")
	E0 <- attr(E(g2), "vnames")
	E(g)$color<- ifelse(E1 %in% E0, "blue", "gray")
	V(g)$color <- colorMatch(ug, g)[[1]]
	return(graph = g)
}

#' @title Vertex and edge graph coloring on the base of fitting
#'
#' @description Add vertex and edge color atrributes to an igraph object,
#' based on a fitting results data.frame generated by
#' \code{\link[SEMgraph]{SEMrun}}.
#' @param est A data.frame of estimated parameters and p-values, derived
#' from the \code{fit} object returned by \code{\link[SEMgraph]{SEMrun}}.
#' As an alternative, the user may provide a "gest" or "dest" data.frame
#' generated by \code{\link[SEMgraph]{SEMrun}}.
#' @param graph An igraph object.
#' @param group group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects.
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}. By default, method is set
#' to "none" (i.e., no multiple test correction).
#' @param alpha Significance level for node and edge coloring
#' (by default, \code{alpha = 0.05}).
#' @param vcolor A vector of three color names. The first color is given
#' to nodes with P-value < alpha and beta < 0, the third color is given
#' to nodes with P-value < alpha and beta > 0, and the second is given
#' to nodes with P-value > alpha. By default,
#' \code{vcolor = c("lightblue", "white", "pink")}.
#' @param ecolor A vector of three color names. The first color is given
#' to edges with P-value < alpha and regression coefficient < 0, the
#' third color is given to edges with P-value < alpha and regression
#' coefficient > 0, and the second is given to edges with P-value > alpha.
#' By default, \code{vcolor = c("blue", "gray50", "red2")}.
#' @param ewidth A vector of two values. The first value refers to the
#' basic edge width (i.e., edges with P-value > alpha), while the second
#' is given to edges with P-value < alpha. By default ewidth = c(1, 2).
#' @param ... Currently ignored.
#'
#' @import igraph
#' @export
#'
#' @return An igraph object with vertex and edge color and width attributes.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' \donttest{
#'
#' # Model fitting: node perturbation
#' sem1 <- SEMrun(graph = alsData$graph, data = alsData$exprs,
#'                group = alsData$group,
#'                fit = 1)
#' est1 <- parameterEstimates(sem1$fit)
#'
#' # Model fitting: edge perturbation
#' sem2 <- SEMrun(graph = alsData$graph, data = alsData$exprs,
#'                group = alsData$group,
#'                fit = 2)
#' est20 <- subset(parameterEstimates(sem2$fit), group = 1)[, -c(4, 5)]
#' est21 <- subset(parameterEstimates(sem2$fit), group = 2)[, -c(4, 5)]
#'
#' # Graphs
#' g <- alsData$graph
#' x <- alsData$group
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(1,4))
#' gplot(colorGraph(est = est1, g, group = x, method = "BH"),
#'       main = "vertex differences")
#' gplot(colorGraph(est = sem2$dest, g, group = NULL),
#'       main = "edge differences")
#' gplot(colorGraph(est = est20, g, group = NULL),
#'       main = "edges for group = 0")
#' gplot(colorGraph(est = est21, g, group = NULL),
#'       main = "edges for group = 1")
#' par(old.par)
#'
#' }
#'
colorGraph <- function (est, graph, group, method = "none", alpha = 0.05,
                       vcolor = c("lightblue","white", "pink"),
                       ecolor = c("royalblue3", "gray50", "red2"),
                       ewidth = c(1, 2), ...)
{
	E(graph)$color <- "gray50"
	if (!is.null(group)) {

	if (colnames(est)[4] == "est") {
		B <- est[est$op == "~",]
		G <- B[B$rhs == "group",]
		B <- B[-c(1:nrow(G)),]
		vnames <- gsub("z", "", G$lhs)
		G$pvalue <- p.adjust(G$pvalue, method = method)
		Vr <- vnames[G$pvalue < alpha & G$est < 0]
		Va <- vnames[G$pvalue < alpha & G$est > 0]
		V(graph)$color <- ifelse(V(graph)$name %in% Vr, vcolor[1],
		                         ifelse(V(graph)$name %in% Va,
		                         vcolor[3], vcolor[2]))
		enames <- gsub("z", "", paste0(B$rhs, "|", B$lhs))
		B$pvalue <- p.adjust(B$pvalue, method = method)
		Er <- enames[B$pvalue < alpha & B$est < 0]
		Ea <- enames[B$pvalue < alpha & B$est > 0]
		E(graph)$color <- ifelse(attr(E(graph), "vnames") %in%
		Er, ecolor[1], ifelse(attr(E(graph), "vnames") %in%
		                      Ea, ecolor[3], ecolor[2]))
		E(graph)$width <- ifelse(E(graph)$color == ecolor[2],
		                         ewidth[1], ewidth[2])

	} else if (colnames(est)[2] == "Stat") {
			G <- est[-c(1:3),]
			vnames <- gsub("X", "", rownames(G))
			G$pvalue <- p.adjust(G$pvalue, method = method)
			Vr <- vnames[G$pvalue < alpha & G$Stat < 0]
			Va <- vnames[G$pvalue < alpha & G$Stat > 0]
			V(graph)$color <- ifelse(V(graph)$name %in% Vr, vcolor[1],
			                         ifelse(V(graph)$name %in% Va,
			                                vcolor[3], vcolor[2]))
		}
	}
	if (is.null(group)) {
		if (colnames(est)[4] == "est") {
			B <- est[est$op == "~", ]
			enames <- gsub("z", "", paste0(B$rhs, "|", B$lhs))
			B$pvalue <- p.adjust(B$pvalue, method = method)
			Er <- enames[B$pvalue < alpha & B$est < 0]
			Ea <- enames[B$pvalue < alpha & B$est > 0]
			E(graph)$color <- ifelse(attr(E(graph), "vnames") %in%
			   Er, ecolor[1], ifelse(attr(E(graph), "vnames") %in%
			                         Ea, ecolor[3], ecolor[2]))
			E(graph)$width <- ifelse(E(graph)$color == ecolor[2],
			                         ewidth[1], ewidth[2])

		} else if (colnames(est)[4] == "d_est") {
			B <- est[est$op == "~", ]
			enames <- gsub("z", "", paste0(B$rhs, "|", B$lhs))
			B$pvalue<- p.adjust(B$pvalue, method=method)
			Er <- enames[B$pvalue < alpha & B$d_est < 0]
			Ea <- enames[B$pvalue < alpha & B$d_est > 0]
			E(graph)$color <- ifelse(attr(E(graph), "vnames") %in%
			   Er, ecolor[1], ifelse(attr(E(graph), "vnames") %in%
			                         Ea, ecolor[3], ecolor[2]))
			E(graph)$width <- ifelse(E(graph)$color == ecolor[2],
			                         ewidth[1], ewidth[2])
		}
	}
	return(graph)
}

colorMatch <- function(g1, g2, ...)
{
	if (!is.null(V(g1)$color)) {
		idx <- match(V(g2)$name, V(g1)$name)
		Vcol <- V(g1)$color[idx]
	} else {
		Vcol <- rep("white", vcount(g2))
	}
	if (!is.null(E(g1)$color)) {
		idx <- match(attr(E(g2), "vnames"), attr(E(g1), "vnames"))
		Ecol <- E(g1)$color[idx]
	} else {
		Ecol <- rep("gray60", ecount(g2))
	}
	if (!is.null(E(g1)$width)) {
		idx <- match(attr(E(g2), "vnames"), attr(E(g1), "vnames"))
		Ewid <- E(g1)$width[idx]
	} else {
		Ewid <- rep(1, ecount(g2))
	}
	return(list(Vcol, Ecol, Ewid))
}

Brown.test <- function(x, p, theta = NULL, tail = "both", ...)
{
	# From two-sided to one-sided (positive or negative) tests
	if (tail == "positive") p <- ifelse(theta > 0, p/2, 1 - p/2)
	if (tail == "negative") p <- ifelse(theta > 0, 1 - p/2, p/2)

	# Fisher's (1932, 4th ed.) combined X2 test
	if (is.null(x)) return(1 - pchisq(q = -2*sum(log(p)), df = 2*length(p)))

	# Brown's (1975) combined X2 test
	tmp <- c(-2.59, -2.382, -2.17, -1.946, -1.709, -1.458, -1.194,
	         -0.916, -0.625, -0.320, 0, 0.334, 0.681, 1.044,
	         1.421, 1.812, 2.219, 2.641, 3.079, 3.531, 4)

	s2X2 <- 4*ncol(x) + 2*sum(approx(seq(-1, 1, .1), tmp,
	                      xout = cor(x)[which(as.vector(lower.tri(cor(x))))])$x)
	EX2 <- 2*ncol(x)

	fX2 <- -2*sum(log(p))
	pX2 <- 1 - pchisq(q = fX2/(s2X2/(2*EX2)), df = 2*EX2^2/s2X2)

	return(pX2)
}

#' @title Optimal model search strategies
#'
#' @description Four model search strategies are implemented combining
#' \code{SEMdag()}, \code{SEMbap()}, and \code{extendGraph()} functions.
#' All strategies estimate a DAG through the adjusted (de-correlate)
#' data matrix Z by iteratively update DAG and Z.
#' @param graph Input network as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param gnet Reference directed network used to validate and import
#' nodes and interactions.
#' @param search Search strategy. Four model search strategies are available:
#' \itemize{
#' \item "outer". The estimated DAG is extended using
#' \code{\link[SEMgraph]{extendGraph}} to find new indirect paths (i.e.,
#' inferred directed connections that may hide new mediators). New
#' interactions and mediators will be searched and imported from the
#' reference network (argument \code{gnet}, see above). Both DAG and
#' extended graph complexity can be controlled with \code{beta} > 0 and
#' \code{d} > 1 arguments, respectively (see below). The term "outer"
#' means that new model mediator variables are imported from an external
#' resource (i.e., the reference network).
#' \item "inner". This strategy is analogous to the "outer" one,
#' but disables external mediator search. In other words, new indirect
#' paths are generated by adding new interactions of the input model, so
#' that mediators will be nodes already present in the input graph. The
#' reference network is still used to validate new model paths. Also in
#' this case, \code{beta} > 0 and \code{d} > 1 are used.
#' \item "direct". The input graph structure is improved through direct
#' (i.e., adjacent) link search, followed by interaction validation and
#' import from the reference network, with no mediators
#' (i.e., \code{d = 1}).
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
#' by the argument gnet. If the edges of the reference interactome are
#' weighted by P-value, as defined by the \code{E(graph)$pv} attribute,
#' the shortest path with the smallest sum of weights will be chosen (e.g.,
#' see \code{\link[SEMgraph]{weightGraph}} for graph weighting options).
#' @param alpha Significance level for false discovery rate (FDR) used
#' for either local d-separation tests (below \code{limit}) or conditional
#' independence (CI) test (above \code{limit}). This argument is used to
#' control data de-correlation. A higher \code{alpha} level includes more
#' hidden covariances, thus considering more sources of confounding.
#' If \code{alpha = 0}, data de-correlation is disabled.
#' By default, \code{alpha = 0.05}.
#' @param pstop A logical value. With the argument \code{pstop = TRUE}
#' (default), the algorithm can be halted when the Shipley's global model
#' test P-value > 0.05. If \code{pstop = FALSE}, the model search
#' algorithm stops when no additional edges can be added to the estimated
#' DAG.
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
#' any reference and the output model structure will be completely
#' determined by data.
#' Output model complexity can be limited using arguments \code{d} and
#' \code{beta}.
#' While d is fixed to 0 or 1 in "basic" or "direct", respectively;
#' we suggest starting with \code{d = 2} (only one mediator)
#' for the other two strategies.
#' For knowledge-based strategies, we suggest to to start with
#' \code{beta = 0.1}. Then, beta can be relaxed (0 to < 0.1) to improve
#' model fitting, if needed. Since data-driven models can be complex,
#' we suggest to start from beta = 0.1 when using the "basic" strategy.
#' The beta value can be relaxed until a good model fit is obtained.
#' Argument alpha determines the extent of data adjustment: lower alpha
#' values for FDR correction correspond to a smaller number of significant
#' confounding factors, hence a weaker correction
#' (default \code{alpha = 0.05}).
#'
#' @import lavaan
#' @import igraph
#' @import GGMncv
#' @importFrom glmnet glmnet
#' @importFrom RcppEigen fastLm
#' @importFrom stats na.omit var cov qchisq pchisq p.adjust
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom flip flip plot
#' @export
#'
#' @return The output model as well as the adjusted dataset are returned
#' as a list of 3 objects:
#' \itemize{
#' \item "fit", the fitted output model (lavaan object);
#' \item "graph", the output model as an igraph object;
#' \item "data", the adjusted dataset.
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
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
#'       d = 2, search = "inner", beta = 0.05, alpha = 0.05)
#' m3 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
#'       d = 2, search = "outer", beta = 0.05, alpha = 0.05)
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
modelSearch <- function(graph, data, gnet = NULL, d = 2, search = "basic",
                        beta = 0, alpha = 0.05, pstop = TRUE, limit = 30000,
                        verbose = FALSE, ...)
{

	# Set graph and data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	dataY <- as.matrix(data[, nodes])

	# Stepwise search
	if (search == "basic") d <- 0
	if (search == "direct") d <- 1
	Gt <- ig
	Zt <- dataY
	k <- 1

	while (k > 0) {
		Zt1 <- quiet(SEMbap(Gt, Zt, method = "BH", alpha = alpha,
		                    limit = limit))

		if (is.null(Zt1)) {
			cat("k =", k, "Searching for missing covariances ... 0\n")
			Zt1$data <- Zt
			Gt1 <- quiet(SEMdag(Gt, Zt1$data, gnet = gnet, d = d,
			                    beta = beta)$dag)
			break
		}

		nE <- ecount(Zt1$guu)
		cat("k =", k, "Searching for missing covariances ...", nE, "\n")
		Gt1 <- quiet(SEMdag(Gt, Zt1$data, gnet = gnet, d = d, beta = beta)$dag)
		if (pstop) {
			ctest <- Shipley.test(Gt1, Zt1$data, limit = limit, verbose = FALSE)
			pv <- ctest$ctest[3]
		} else {
			pv <- 0
		}

		if (pv > 0.05 | ecount(Gt1) - ecount(Gt) == 0) break
		Gt <- Gt1
		Zt <- Zt1$data
		k <- k + 1
	}
	cat("Done.\n\n")

	E1 <- attr(E(Gt1), "vnames")
	E0 <- attr(E(ig), "vnames")
	E(Gt1)$color <- ifelse(E1 %in% E0, "blue", "red")
	Gt1.red <- Gt1 - E(Gt1)[which(E(Gt1)$color == "blue")]
	if(ecount(Gt1.red) == 0) {
		return(message("ERROR: no new edges inferred.
		 Try decreasing the beta threshold !"))
	}


	# Extended graph

	if (search == "basic") {
		dataZ <- Zt1$data
		ig <- Gt1

	} else if (search != "basic") {

		if (search == "direct" | search == "inner") {
			ext <- extendGraph(g = list(Gt1,Gt1.red), Zt1$data,
			                   gnet = gnet,
			                   verbose = verbose)

		} else if (search == "outer") {
			ext <- extendGraph(g = list(Gt1, Gt1.red), data, gnet = gnet,
							   verbose = verbose)
		}

		sel <- which(V(ext$Ug)$color == "yellow")
		yellow <- names(V(ext$Ug)[sel])
		new <- setdiff(yellow, V(ig)$name)
		V(ext$Ug)$color[sel] <- ifelse(yellow %in% new, "yellow", "green")
		dataZ <- cbind(Zt1$data, data[, new])
		ig <- ext$Ug
	}

	if (verbose) gplot(ig, main = "Estimated Extended Graph")
	if (!pstop) {
		C_test <- Shipley.test(ig, dataZ, limit = limit, verbose = TRUE)
	} else {
		print(data.frame(C_test = ctest$ctest[1],
		                 df = ctest$ctest[2],
		                 pvalue = round(ctest$ctest[3], 6)))
	}
	cat("\n")
	fit <- SEMrun(ig, dataZ, algo = "ricf")

	return(list(fit = fit$fit, graph = ig, data = dataZ))
}

#' @title Perturbed path search utility
#'
#' @description This function uses \code{\link[SEMgraph]{SEMace}} to find
#' significant causal effects between source-sink pairs and
#' \code{\link[SEMgraph]{SEMpath}} to fit them and test their edge
#' perturbation.
#' @param graph Input network as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. Group specification enables edge perturbation
#' testing. By default, \code{group = NULL}.
#' @param ace A data.frame generated by \code{\link[SEMgraph]{SEMace}}.
#' If NULL, \code{\link[SEMgraph]{SEMace}} will be automatically run.
#' @param path If \code{path = "directed"}, all directed paths between
#' the two nodes will be included in the fitted model.
#' If \code{path = "shortest"}, only shortest paths will be considered.
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}.
#' By default, \code{method = "none"} (i.e., no multiple test correction).
#' @param alpha Significance level for ACE selection (by default,
#' \code{alpha = 0.05}).
#' @param verbose Show the significant directed (or shortest) paths
#' inside the input graph.
#' @param ... Currently ignored.
#'
#' @import lavaan
#' @import igraph
#' @importFrom stats approx
#' @importFrom dagitty paths
#' @importFrom utils head
#' @export
#'
#' @return A list of 3 objects:
#' \itemize{
#' \item "paths", list of paths as igraph objects;
#' \item "fit", fitting results for each path as a lavaan object;
#' \item "dfp", a data.frame containing SEM global fitting statistics.
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @examples
#'
#' \donttest{
#'
#' # Find and evaluate significantly perturbed paths
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' adjData <- SEMbap(graph = alsData$graph, data = als.npn)$data
#'
#' ace <- SEMace(graph = alsData$graph, data = adjData,
#'               group = alsData$group)
#' ace <- ace[order(ace$pvalue),]
#' print(ace)
#'
#' paths <- pathFinder(graph = alsData$graph, data = adjData,
#'                     group = alsData$group,
#'                     ace = ace)
#'
#' head(parameterEstimates(paths$fit$P19))
#' gplot(paths$paths$P19)
#'
#' path19 <- SEMpath(graph = alsData$graph, data = adjData,
#'                   group = alsData$group,
#'                   from = "7133",
#'                   to = "4747",
#'                   path = "directed",
#'                   verbose = TRUE)
#'
#' }
#'
pathFinder <- function(graph, data, group = NULL, ace = NULL, path = "directed",
                       method = "none", alpha = 0.05, verbose = FALSE, ...)
{
	if (is.null(ace)) {
		ace <- SEMace(graph, data, group, method = method, alpha = alpha)
	}
	ace <- ace[ace$pvalue < alpha,]
	ace <- ace[order(ace$pvalue),]
	sources <- as.character(ace$source)
	sinks <- as.character(ace$sink)
	paths <- list()
	lav <- list()
	res <- NULL
	N <- nrow(ace)
	cat("\nFound", N, "significant ACEs.\n\n")
	if (N == 0) return(list(paths = NULL, fit = NULL, dfp = NULL))
	for (i in 1:N) {
		fit <- quiet(SEMpath(graph, data, group, from = sources[i],
		                     to = sinks[i],
		                     path = path,
		                     verbose = verbose))
		if(is.null(fit)) next
		if (!is.null(group) & vcount(fit$graph) > 100) {
			dev_df <- fit$fit$ricf$dev/fit$fit$ricf$df
			srmr <- fit$fit$fitIdx$srmr
			pv1 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$Stat,
			                  tail = "positive")
			pv2 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$Stat,
			                  tail = "negative")
		} else {
			dev_df <- fitMeasures(fit$fit, "chisq")/fitMeasures(fit$fit, "df")
			srmr <- fitMeasures(fit$fit, "srmr")
			pv1 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$est,
			                  tail = "positive")
			pv2 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$est,
			                  tail = "negative")
		}

		dfp <- data.frame(pathId = paste0("P", rownames(ace)[i]),
		                  sink = sinks[i],
		                  op = "<-",
		                  source = sources[i],
		                  N.nodes = vcount(fit$graph),
		                  N.edges = ecount(fit$graph),
		                  dev_df = round(dev_df, 3),
		                  srmr = round(srmr, 3),
		                  pv.act = round(pv1, 6),
		                  pv.inh = round(pv2, 6))

		res <- rbind(res, dfp)
		paths[[i]] <- fit$graph
		lav[[i]] <- fit$fit
	}
	rownames(res) <- NULL
	names(paths) <- res$pathId
	names(lav) <- res$pathId
	print(head(res))
	return(list(paths = paths, fit = lav, dfp = res))
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
#' @import igraph
#' @import lavaan
#' @importFrom stats cor approx
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom ggm fitAncestralGraph
#' @importFrom flip flip plot
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
extractClusters <- function(graph, data, group = NULL, membership = NULL,
                            map = FALSE, verbose = FALSE, ...)
{
	if (is.null(membership)) {
		membership <- clusterGraph(graph, type = "ebc", HM = "none",
		                           size = 5,
		                           verbose = FALSE)
	}
	clusters <- cplot(graph, membership, l=layout.auto, map, verbose)[-1]
	if ("HM9999" %in% names(clusters)) N <- length(clusters) - 1
	N <- length(clusters)
	res <- NULL
	lav <- list()
	for (i in 1:N) {
		fit <- quiet(SEMfit(clusters[[i]], data, group))
		if (is.null(fit)) next
		if (!is.null(group) & vcount(clusters[[i]]) > 100) {
			dev_df <- fit$fit$ricf$dev/fit$fit$ricf$df
			srmr <- fit$fit$fitIdx[3]
			pv1 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$Stat,
			                  tail = "positive")
			pv2 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$Stat,
			                  tail = "negative")
		} else {
			dev_df <- fitMeasures(fit$fit, "chisq")/fitMeasures(fit$fit, "df")
			srmr <- fitMeasures(fit$fit, "srmr")
			pv1 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$est,
			                  tail = "positive")
			pv2 <- Brown.test(x = fit$dataXY[, -1], p = fit$gest$pvalue,
			                  theta = fit$gest$est,
			                  tail = "negative")
		}

		dfc <- data.frame(cluster = names(clusters)[i],
		                  N.nodes = vcount(clusters[[i]]),
		                  N.edges = ecount(clusters[[i]]),
		                  dev_df = round(dev_df, 3),
		                  srmr = round(srmr, 3),
		                  pv.act = round(pv1, 6),
		                  pv.inh = round(pv2, 6))

		res <- rbind(res, dfc)
		lav[[i]] <- fit$fit
	}
	rownames(res) <- NULL
	names(lav) <- names(clusters)
	print(res)
	return(list(clusters = clusters, fit = lav, dfc = res))
}

#' @title Pairwise plotting of multivariate data
#'
#' @description Display a pairwise scatter plot of two datasets for a
#' random selection of variables. If the second dataset is not given,
#' the function displays a histogram with normal curve superposition.
#'
#' @param x A matrix or data.frame (n x p) of continuous data.
#' @param y A matrix or data.frame (n x q) of continuous data.
#' @param size number of rows to be sampled (default \code{s = nrow(z)}).
#' @param r number of rows of the plot layout (default \code{r = 4}).
#' @param c number of columns of the plot layout (default \code{r = 4}).
#' @param ... Currently ignored.
#'
#' @importFrom graphics par hist curve legend
#' @importFrom stats dnorm
#' @export
#'
#' @return No return value
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#' adjdata <- SEMbap(sachs$graph, log(sachs$pkc))$data
#' rawdata <- log(sachs$pkc)
#' pairwiseMatrix(adjdata, rawdata, size = 1000)
#'
pairwiseMatrix<- function (x, y = NULL, size = nrow(x), r = 4, c = 4, ...)
{
    if (r * c > ncol(x)) {
  	 r <- r - 1
	 c <- c - 1
	 p <- 1:(r * c)
    }else{
	 p <- sample(1:ncol(x), size = r * c)
	}
    n <- sample(1:nrow(x), size = size)
    vnames <- colnames(x)
    if (is.null(y)) {
        xx <- x[, vnames]
		old.par <- par(no.readonly = TRUE)
		par(mfrow = c(r, c), mar = rep(3, 4))
		for (j in p) {
            h <- hist(xx[n, j], breaks = 30, freq = FALSE, col = "lightblue", main = vnames[j])
            x <- seq(-4, +4, by = 0.02)
            curve(dnorm(x), add = TRUE, col = "blue", lwd = 2)
        }
		on.exit(par(old.par))
    }
    else {
        xx <- x[, vnames]
		yy <- y[, vnames]
		old.par <- par(no.readonly = TRUE)
        par(mfrow = c(r, c), mar = rep(2, 4))
        for (j in p) {
            r <- round(cor(xx[n, j], yy[n, j]), 2)
            plot(xx[n, j], yy[n, j], main = vnames[j])
            legend("topleft", paste0("r = ", r), bty = "n", cex = 1, text.col = "blue")
        }
		on.exit(par(old.par))
    }
}

#' @title Node ancestry utilities
#'
#' @description Get ancestry for a collection of nodes in a graph.
#' These functions are wrappers for the original \code{SEMID} R package.
#'
#' @param g An igraph object.
#' @param nodes the nodes in the graph of which to get the ancestry.
#'
#' @references
#' Rina Foygel Barber, Mathias Drton and Luca Weihs (2019). SEMID:
#' Identifiability of Linear Structural Equation Models. R package
#' version 0.3.2. <https://CRAN.R-project.org/package=SEMID/>
#'
#' @examples
#'
#' # Get all ancestors
#' an <- V(sachs$graph)[ancestors(sachs$graph, "Erk")]; an
#'
#' # Get parents
#' pa <- V(sachs$graph)[parents(sachs$graph, "PKC")]; pa
#'
#' # Get descendants
#' de <- V(sachs$graph)[descendants(sachs$graph, "PKA")]; de
#'
#' # Get siblings
#' sib <- V(sachs$graph)[siblings(sachs$graph, "PIP3")]; sib
#'
#' @return a sorted vector of nodes.
#' @name ancestry

#' @rdname ancestry
#' @export
ancestors <- function(g, nodes)
{
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "in", unreachable = F)$order,
                  na.last = NA))
}

#' @rdname ancestry
#' @export
descendants <- function(g, nodes)
{
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "out", unreachable = F)$order,
                  na.last = NA))
}

#' @rdname ancestry
#' @export
parents <- function(g, nodes)
{
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "in"))))
}

#' @rdname ancestry
#' @export
siblings <- function(g, nodes)
{
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "out"))))
}
