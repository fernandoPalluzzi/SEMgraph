# SEMgraph 0.3.3
#
# This is the SEMgraph package for statistical graph analysis using
# Structural Equation Models
# Typically, SEMgraph requires three input elements:
# 1) Interactome (e.g. KEGG singaling pathways or STRING PPI)
# 2) Quantitative data (e.g. GWAS, DNA-methylation arrays, RNA-seq)
# 3) Case/Control identifiers
#

arrowDirection <- function(ug, dg, ...) {
	ftm <- as_edgelist(ug)
	for (k in 1:nrow(ftm)) { #k=7
		exy <- c(ftm[k, 1], ftm[k, 2])
		eyx <- c(ftm[k, 2], ftm[k, 1])
		exy.ids <- get.edge.ids(dg, exy)
		eyx.ids <- get.edge.ids(dg, eyx)
		if(exy.ids > 0 & eyx.ids == 0) ftm[k, 1:2] <- exy
		if(exy.ids == 0 & eyx.ids > 0) ftm[k, 1:2] <- eyx
	}
	return(graph_from_edgelist(ftm, directed = TRUE))
}

#' @title Graph properties summary and graph correction
#'
#' @description Produces a summary of network properties and returns 
#' the largest network component without self-loops.
#' @param g Input network as an igraph object.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @export
#'
#' @return Network largest component without self-loops.
#'
#' @examples
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#'
properties <- function(g)
{
    if (!is_igraph(g)) ig <- graph_from_graphnel(g) else ig <- g
    ig <- simplify(ig, remove.loops = TRUE)
    gcs <- igraph::decompose.graph(ig)
	vsize <- lapply(1:length(gcs), function(x) vcount(gcs[[x]]))
	ig1 <- gcs[[which.max(vsize)]]
	#ig1
	cat("Frequency distribution of graph components")
	print(table(sapply(gcs, vcount)))
	cat("Percent of vertices in the giant component:",
	    round(100*vcount(ig1)/vcount(ig), 1), "%\n\n")
	
	print(c(is.simple = is_simple(ig1),
	        #is.connected = is_connected(ig1),
	        is.dag = is_dag(ig1),
	        is.directed = is_directed(ig1),
	        is.weighted = is_weighted(ig1)))
		 	
	#require(SEMID)
	#L <- as.matrix(as_adj(ig1))
	#O <- diag(rep(0, vcount(ig1)))
	#print(c(is.globally.indentifiable = SEMID::graphID.globalID(L, O)))
	
	print(c(which.mutual = table(which_mutual(ig1))))
	print(E(ig1)[which_mutual(ig1)])
	
	return(list(ig1 = ig1, gcs = gcs[-which.max(vsize)]))
}

#' @title Network merging utility
#'
#' @description Merge a list of directed and/or undirected networks into one.
#' @param g List of igraph objects to be merged.
#' @param gref Reference input graph, used to define color code: 
#' "yellow", common nodes between merged and reference graphs; 
#' "lightblue", nodes in the merged graph that are absent in reference one. 
#' If gref = NULL (default), no nodes colour will be added.
#' @param gnet External directed interaction network as an igraph object. 
#' When a directed interactome is provided, missing edge directions will 
#' be recovered from it, and a directed output network will be enforced. 
#' If gnet = NULL (default), a merged undirected network will be produced.
#' @param verbose A logical value. If FALSE (default), the output merged 
#' graph will not be plotted to screen.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @export
#'
#' @return A merged network as an igraph object.
#' @seealso \code{\link{igraph}}
#' @examples
#' # Data loading
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#' 
#' # Generating graphs
#' fit <- SEMfit(graph, data, group, B = NULL, perm = 10000)
#' ndf <- vcount(graph)*(vcount(graph) - 1)/2 - ecount(graph)
#' ggm <- SEMggm(fit = fit, gnet = kegg, d = 2, perm = 10000, alpha = 1/ndf)
#' 
#' # Defining graph variables
#' ig <- ggm$graph$ig        # directed graph
#' guu <- ggm$graph$guu      # undirected graph
#' 
#' # Merging without external directed intractome
#' UG <- mergeGraph(g = list(ig, guu))
#' gplot(UG)
#' 
#' # Recover edge directions from external intractome
#' DG <- mergeGraph(g = list(ig, guu), gnet = kegg)
#' gplot(DG)
#'
mergeGraph <- function(g = list(), gref = NULL, gnet = NULL, verbose = FALSE, ...)
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
				g1 <- c(g1, list(arrowDirection(ug = g[[i]], dg = gnet)))
			} else {
				g1 <- c(g1, list(g[[i]]))
			}
		}
	}
	
	Ug <- igraph::graph.union(g1)
	# Vertex color attribute with (UX,UY,LX,LY,LM) nodes
	if (is_igraph(gref)) {
		V(Ug)$color <- ifelse(V(Ug)$name %in% V(gref)$name,
		                      "lightblue",
		                      "yellow")
	} else {
		V(Ug)$color <- "white"
	}
	
	# Edge weight attribute (example con weight = pvalue)
	#W <- NULL
	#for (i in 1:length(g)) {
		#w <- get.edge.attribute(Ug, paste0("weight_", i))
		#Ug[is.na(w)] <- 1
		#W <- cbind(W, w)
		#Ug <- remove.edge.attribute(Ug, paste0("weight_", i))
	#}
	#E(Ug)$weight <- apply(W, 1, min)
	#E(Ug)$weight <- 1
	if (verbose) plot(Ug)
	#gplot(Ug)
	
	return(Ug)
}

#' @title Network nodes merging by a user-defined membership attribute
#'
#' @description Merge groups of network nodes using a custom membership 
#' attribute (e.g., cluster membership).
#' @param graph Network as an igraph object.
#' @param membership Cluster membership. A vector of cluster membership 
#' identifiers, where vector names correspond to network node names. 
#' Network clustering can be done using \code{\link[SEMgraph]{clusterGraph}}.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @importFrom graph combineNodes
#' @export
#'
#' @return A network with merged nodes as an igraph object.
#' @seealso \code{\link[SEMgraph]{clusterGraph}}
#' @examples
#' # Data loading
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#' 
#' # Extracting hidden modules with SEMfsr
#' fsr.uv <- SEMfsr(graph = graph, data = data, group = group, 
#'                  type = "ebc",
#'                  HM = "UV",
#'                  size = 15)
#' 
#' # Merging nodes into modules
#' cg <- mergeNodes(graph = graph, membership = fsr.uv$M)
#'
mergeNodes <- function(graph, membership, ...)
{ 
	# Set membership object
	if (is.numeric(membership)) {
		nodes <- names(membership)
		membership <- paste0("GM", membership)
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
	for (i in 1:length(LM)) { #i=3
		gLMi <- graph::combineNodes(LM[[i]], gLM, names(LM)[i], mean)
		gLM <- gLMi
	}
	
	ig <- graph_from_graphnel(gLM)
	if (length(V(ig)$color) == 0) V(ig)$color <- "white"
	V(ig)$color[substr(V(ig)$name, 2, 2) == "M"] <- "green"
	vcol <- V(ig)$color
	names(vcol) <- V(ig)$name
	gplot(ig)
	
	return(gLM = ig)
}

#' @title Graph plotting with renderGraph
#'
#' @description Wrapper for \code{\link[Rgraphviz]{renderGraph}}.
#' @param x An igraph object.
#' @param psize Node size (by default, psize = 80).
#' @param main plot main title.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @importFrom graph nodes edgeNames isDirected nodeRenderInfo edgeRenderInfo graphRenderInfo
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @export
#'
gplot<- function(x, psize=80, main="", ...)
{
	# set graphNEL object
	g<- as_graphnel(x)
	vcol<- V(x)$color
	vsize<- V(x)$size
	vlab<- V(x)$label
	ecol<- E(x)$color
	elwd<- E(x)$width
	elab<- E(x)$label
	
	if(length(vcol)>0) {
	 names(vcol)<- V(x)$name
	 vcol<- vcol[graph::nodes(g)]
	}
	if(length(vsize)>0) {
	 names(vsize)<- V(x)$name
	 vsize<- 10*vsize[graph::nodes(g)]
	}else{
	 vsize<- rep(psize, length=vcount(x))
	 names(vsize)<- V(x)$name
	}
	if(length(vlab)>0) {
	 names(vlab)<- V(x)$name
	 vlab<- vlab[graph::nodes(g)]
	}
	if(length(ecol)>0) {
	 names(ecol)<- gsub("\\|", "~", attr(E(x), "vnames"))
	 ecol<- ecol[graph::edgeNames(g,recipEdges="distinct")] 
	}
	if(length(elwd)>0) {
	 names(elwd)<- gsub("\\|", "~", attr(E(x), "vnames"))
	 elwd<- elwd[graph::edgeNames(g,recipEdges="distinct")]
	}
	if(length(elab)>0) {
	 names(elab)<- gsub("\\|", "~", attr(E(x), "vnames"))
	 elab<- elab[graph::edgeNames(g,recipEdges="distinct")] 
	}else{
	 elab<- rep("", ecount(x))
	 names(elab)<- gsub("\\|", "~", attr(E(x), "vnames"))
	}
	
	if (graph::isDirected(g)) {
	 g<- Rgraphviz::layoutGraph(g, layoutType="dot", edgeAttrs=list(label=elab))
	}else{
	 g<- Rgraphviz::layoutGraph(g, layoutType="fdp", edgeAttrs=list(label=elab))
	}
	graph::nodeRenderInfo(g)<- list(col="black", fill=vcol, lty=1, label=vlab,lwd=1, 
	        textCol="black", fontsize=15, shape="circle", width=vsize, height=vsize)
	graph::edgeRenderInfo(g)<- list(col=ecol, lty=1, lwd=elwd)
	graph::graphRenderInfo(g)<- list(main=main)
	Rgraphviz::renderGraph(g)
}

#' @title Correlation matrix to graph
#'
#' @description Convert a correlation matrix to an igraph object.
#' @param R Correlation matrix.
#' @param n. Sample size (i.e., the number of subjects).
#' @param alpha Significance level used to compute the correlation threshold.
#' By default, alpha = 0.05.
#' @param type. Correlation type: "marg" for marginal correlation, and 
#' "cond" for conditional correlation.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import lavaan
#' @import igraph
#' @export
#'
#' @return An igraph object.
#'
corr2graph<- function(R, n, alpha=0.05, type="marg", ...)
{
	# select correlation matrix:
	p<- nrow(R)
	if (type == "marg") {
	 q <- 0
	 K <- R
	}else if (type== "cond") {
	 q <- p-2
	 if( corpcor::is.positive.definite(R) ) {
	  K <- corpcor::cor2pcor(R)
	  rownames(K)<- colnames(K)<- rownames(R)
	 } else K <- corpcor::pcor.shrink(R, verbose=TRUE)[1:p,1:p]
	}

	# select the correlation threshold:
	z <- abs(atanh(K[lower.tri(K)]))/sqrt(1/(n-3-q))
	pBH <- p.adjust(2*(1- pnorm(z)), "BH")
	Z <- min( z[which(pBH < alpha)] )
	#Z <- qnorm(alpha/2, lower.tail=FALSE) #1.959964
	thr <- (exp(2*Z/sqrt(n-3-q))-1)/(exp(2*Z/sqrt(n-3-q))+1) #thr

	# from adjacency matrix=A -> to=graph:
	A<- ifelse(abs(K) > thr, 1, 0)
	diag(A)<- 0 #sum(A)/2
	del<- which( colSums(A) == 0 )
	if( length(del)>0 ) A<- A[-del,-del]
	graph<- graph_from_adjacency_matrix(A, mode="undirected")
	plot(graph)
	
	return( graph )
}

#' @title Path diagram to graph
#'
#' @description Convert a path diagram, specified using lavaan syntax, 
#' to an igraph object.
#' @param model Path diagram using lavaan syntax.
#' @param directed Logical value. If TRUE (default), edge directions from 
#' the path diagram will be preserved. If FALSE, the resulting graph will 
#' be undirected.
#' @param psi Logical value. If TRUE (default) covariances will be converted 
#' into bidirected graph edges. If FALSE, covariances will be excluded from 
#' the output graph.
#' @param verbose Logical value. If TRUE (default), a plot of the output 
#' graph will be generated. For large graphs, this could significantly 
#' increase computation time. If FALSE, graph plotting will be disabled.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import lavaan
#' @import igraph
#' @export
#'
#' @return An igraph object.
#'
#' @examples
#' 
#' # Writing path diagram in lavaan syntax
#' 
#' model <- '
#' # path diagram
#' y1 ~ x1
#' y2 ~ x1
#' y3 ~ y1 + y2
#' y4 ~ x1 + y1 + y3
#' # covariances
#' y1 ~~ y2'
#' 
#' # Converting path diagram into graph
#' graph <- lavaan2graph(model)
#' graph
#'
lavaan2graph<- function(model, directed=TRUE, psi=TRUE, verbose=TRUE, ...)
{
	lav<- lavParTable(model, fixed.x=FALSE)
	lavb<- subset(lav, lav$op == "~")
	lavc<- subset(lav, lav$op == "~~" & (lav$rhs != lav$lhs))
	ftm<- data.frame(cbind(from=lavb$rhs, to=lavb$lhs, label=lavb$label), color="gray50")
	if(nrow(lavc) != 0 & psi == TRUE) {
	 ftmc1<- data.frame(cbind(from=lavc$rhs, to=lavc$lhs, label="", color="blue"))
	 ftmc2<- data.frame(cbind(from=lavc$lhs, to=lavc$rhs, label="", color="blue"))
	 ftm<- rbind(ftm,ftmc1,ftmc2)
	}
	graph<- graph_from_data_frame(ftm, directed=directed)
	if (verbose) plot(graph)
	return( graph )
}

#' @title Graph to dagitty DAG format
#'
#' @description Convert an igraph object to a dagitty directed acyclic 
#' graph (DAG).
#' @param graph igraph object.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @importFrom dagitty graphLayout
#' @export
#'
#' @return An dagitty graph object.
#'
graph2dagitty<- function(graph, ...)
{
	ed <- attr(E(graph), "vnames")[which_mutual(graph) == FALSE]
	eb <- attr(E(graph), "vnames")[which_mutual(graph) == TRUE]
	de <- paste(gsub("\\|", "->", ed), collapse = "\n")
	if(length(eb) == 0) {
		dag<- paste0("dag {\n",de,"\n}") #cat(dag)
	} else {
		eb <- eb[1:(length(eb)/2)]
		be <- paste(gsub("\\|", "<->", eb), collapse = "\n")
		dag <- paste0("dag {\n", de, "\n", be, "\n}") #cat(dag)
	}
	plot(dagitty::graphLayout(dag))
	return(dag)
}

#' @title Graph to path diagram
#'
#' @description Convert an igraph object to a path diagram, specified 
#' using lavaan syntax.
#' @param graph An igraph object.
#' @param nodes Subset of nodes to be included in the path diagram.
#' By default, all the graph nodes will be included inside the path diagram.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import lavaan
#' @import igraph
#' @export
#'
#' @return A path diagram in lavaan syntax.
#'
graph2lavaan<- function(graph, nodes=V(graph)$name, ...)
{
	# Set from-to-matrix representation of edge links
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	ftm<- as_data_frame(ig) #head(ftm)
		
	if( is.directed(ig) & sum(which_mutual(ig))>0 ){
	 sel<- as.numeric(c(E(ig)[which_mutual(ig)]))
	 ftm<- as_data_frame(ig)[-sel,]
	 ubg<- as.undirected(graph_from_data_frame(as_data_frame(ig)[sel,]))
	 ftb<- as_data_frame(ubg)
	}else{ ftb<- NULL }

	modelY<- modelV<- vector()
	if ( is.directed(ig) ) {
	 for(j in 1:nrow(ftm)) modelY[j]<- paste0(ftm[j,2],"~",ftm[j,1])
	 if ( length(ftb)>0  ) for(k in 1:nrow(ftb)) modelV[k]<- paste0(ftb[k,2],"~~",ftb[k,1])
	}else{
	 for(j in 1:nrow(ftm)) modelY[j]<- paste0(ftm[j,2],"~~",ftm[j,1])
	}
	model<- paste(c(sort(modelY), modelV), collapse="\n")
	return( model )
}	

#' @title Convert directed graphs to bow-free acyclic path diagrams (BAPs)
#'
#' @description Convert cycle-containing directed graphs into BAPs.
#' When converting a graph to a SEM, it is desirable to have a recursive
#' (i.e., acyclic) model, that can be viewed as a directed acyclic graph
#' (DAG). BAPs are DAGs extension, allowing the presence of bidirected 
#' edges that can be interpreted as covariances (or latent variables
#' connecting pairs of nodes). Although cycles are permitted by SEM
#' theory, it is possible to convert them into BAPs. The graph2bap 
#' function enables topological sorting of the input graph until a cycle
#' is found. The first node of the cycle is defined by the hierarchical 
#' order, and thus the last edge of the cycle is determined (namely, 
#' the terminal edge). The terminal directed edge of each cycle is changed 
#' into a bidirected one, converting the input directed graph to a BAP.
#' @param graph A directed graph as an igraph object.
#' @param ... arguments to be passed to or from other methods.
#'
#' @import igraph
#' @export
#'
#' @return A BAP as an igraph object.
#'
graph2bap <- function(graph, ...)
{
	# DAG partial topological sort with a warning issue
	idx <- as.numeric(topo_sort(graph, mode = "out"))
	sg1 <- induced_subgraph(graph, vids = V(graph)$name[-idx])
	# Find cycles and reverse terminal edges
	cycles <- FindCycles(sg1)
	add <- function(x) {
		p <- length(x)
		edge <- c(x[p], x[p-1])
	}
	sg2 <- add_edges(sg1, unlist(unique(lapply(cycles, add))))
	return(bap = graph.union(g = list(graph, sg2)))
}

FindCycles <- function(g)
{
	Cycles <- NULL
	for(v1 in V(g)) {
		if(degree(g, v1, mode = "in") == 0) next
		GoodNeighbors <- neighbors(g, v1, mode = "out")
		GoodNeighbors <- GoodNeighbors[GoodNeighbors > v1]
		for(v2 in GoodNeighbors) {
			TempCyc <- lapply(all_simple_paths(g, v2,v1, mode = "out"),
			                  function(p) c(v1, p))
			TempCyc <- TempCyc[which(sapply(TempCyc, length) > 3)]
			TempCyc <- TempCyc[sapply(TempCyc, min) == sapply(TempCyc, "[", 1)]
			Cycles <- c(Cycles, TempCyc)
		}
	}
	return(Cycles)
}
