#  SEMgraph library
#  Copyright (C) 2019-2023 Mario Grassi; Fernando Palluzzi; Barbara Tarantino
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

#' @title SEM-based gene set analysis (GSA)
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
#' fitting has been combined with up- or down-regulation of genes from KEGG
#' to obtain overall pathway perturbation as follows: 
#' \itemize{
#' \item The node perturbation is defined as activated when the minimum among
#' the p-values is positive; if negative, the status is inhibited. 
#' \item Up- or down- regulation of genes (derived from KEGG database) has
#' been obtained from the weighted adjacency matrix of each pathway as column
#' sum of weights over each source node. If the overall sum of node weights
#' is below 1, the pathway is flagged as down-regulated otherwise as up-regulated. 
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
#' Grassi, M., Tarantino, B. SEMgsa: topology-based pathway enrichment analysis with
#' structural equation models. BMC Bioinformatics 23, 344 (2022).
#' https://doi.org/10.1186/s12859-022-04884-8
#' 
#' @examples
#'
#' \dontrun{
#'
#' # Nonparanormal(npn) transformation
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' # Selection of FTD-ALS pathways from kegg.pathways.Rdata
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
	 cpval <- fit$pval
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
	   PVAL = 2 * min(cpval[1], cpval[2]))
	 res.tbl <- rbind(res.tbl,df)
	}
 
 	ADJP <- p.adjust(res.tbl$PVAL, method="bonferroni")
	res.tbl<- cbind(res.tbl, ADJP)
	rownames(res.tbl) <- names(DEG) <- names(g)
	gsa <- res.tbl[order(res.tbl$ADJP),]
	DEG <- DEG[rownames(gsa)]
	
	return( list(gsa=gsa, DEG=DEG) )
}

#' @title Gene set analysis extraction utility
#'
#' @description Extract differentially regulated nodes (DRNs) and interaction 
#' network(s) from the \code{\link[SEMgraph]{SEMgsa}} output.
#' The interaction network is obtained by the union of all the perturbed 
#' pathways.
#'
#' @param gsa Gene Set Analysis object generated by 
#' \code{\link[SEMgraph]{SEMgsa}}.
#' @param reference Reference list of biological pathways. Every element of the 
#' list must be an igraph object. SEMgraph provides both Reacrome and KEGG 
#' pathways list. By default, reactome.pathways is used.
#' @param lcc A logical value. If TRUE (default) only the largest connected 
#' component of the resulting network is considered.
#' @param min.size Minimum required pathway size (default = 5).
#' @param max.size Maximum allowed pathway size (default = 300).
#' @param min.degs Minimum required number of DRNs (default = 2).
#' @param act.alpha Pathway activation significance level (default = 0.05).
#' @param inh.alpha Pathway inhibition significance level (default = 0.05).
#' @param adj.alpha Global cutoff over the adjusted perturbation p-value 
#' (default = 0.05).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of two objects:
#' \enumerate{
#' \item "model", The output model as an igraph object;
#' \item "clusters", A list containing model clustering results.
#' It includes tree more objects:
#' \itemize{
#' \item "membership", A vector reporting cluster membership for each node;
#' \item "fit", hidden model fitting as a lavaan object;
#' \item "dataHM", data matrix containing cluster scores (hidden variables) 
#' alongside the original data variables.
#' }
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#' 
#' @references
#'
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
#'
#' @seealso \code{\link[SEMgraph]{SEMgsa}}
#'
#' @examples
#'
#' \donttest{
#'
#' # Nonparanormal(npn) transformation
#' 
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' # Selection of FTD-ALS pathways from kegg.pathways.Rdata
#' 
#'  paths.name <- c("MAPK signaling pathway",
#'                  "Protein processing in endoplasmic reticulum",
#'                  "Endocytosis",
#'                  "Wnt signaling pathway",
#'                  "Neurotrophin signaling pathway",
#'                  "Amyotrophic lateral sclerosis")
#' 
#' j <- which(names(kegg.pathways) %in% paths.name)
#' 
#' GSA <- SEMgsa(kegg.pathways[j], als.npn, alsData$group,
#'                method = "bonferroni", alpha = 0.05,
#'                n_rep = 1000)
#' 
#' gsa.pert <- gsa.extract(GSA, kegg.pathways)
#' 
#' }
#'
gsa.extract <- function(gsa, reference = reactome.pathways, lcc = TRUE,
                        min.size = 5, max.size = 300, min.drns = 2,
                        act.alpha = 0.05, inh.alpha = 0.05, adj.alpha = 0.05,
                        ...)
{  
  message(paste0("# Evaluated pathways: ", nrow(gsa$gsa), "\n"))
  
  pathways <- gsa$gsa[gsa$gsa$No.nodes > min.size - 1
                      & gsa$gsa$No.nodes < max.size + 1
                      & gsa$gsa$No.DEGs > min.drns - 1
                      & gsa$gsa$pNa < act.alpha
                      & gsa$gsa$pNi < inh.alpha
                      & gsa$gsa$ADJP < adj.alpha,]
  
  pathways <- rownames(pathways)
  message(paste0("# Perturbed pathways: ", length(pathways), "\n"))
  
  j <- which(names(reference) %in% pathways)
  gsa.net <- graph.union(reference[j])
  if (lcc) gsa.net <- properties(gsa.net)[[1]]
  
  seeds <- unique(unlist(gsa$DEG[names(gsa$DEG) %in% pathways]))
  names(seeds) <- NULL
  
  message(paste0("\n# Differentially regulated nodes (seeds): ",
                 length(seeds), "\n"))
  
  return(list(gsa.seeds = seeds, gsa.net = gsa.net))
}

#' @title SEM-based differential causal inference (DCI)
#'
#' @description Creates a network with perturbed edges obtained from the
#' output of \code{\link[SEMgraph]{SEMace}}, comparable to the procedure in
#' Jablonski et al (2022), or of \code{\link[SEMgraph]{SEMrun}} with two-group
#' and CGGM solver, comparable to the algorithm 2 in Belyaeva et al (2021). 
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
#' @param type  Average Causal Effect (ACE) with two-group, "parents" (back-door)
#' adjustement set, and "direct" effects (\code{type = "ace"}, default), or
#' CGGM solver with two-group using a clustering method.
#' If \code{type = "tahc"}, network modules are generated using the tree
#' agglomerative hierarchical clustering method. Other non-tree clustering
#' methods from igraph package include: "wtc" (walktrap community structure
#' with short random walks), "ebc" (edge betweeness clustering), "fgc"
#' (fast greedy method), "lbc" (label propagation method), "lec" (leading
#' eigenvector method), "loc" (multi-level optimization), "opc" (optimal
#' community structure), "sgc" (spinglass statistical mechanics), "none"
#' (no breaking network structure into clusters).
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
#' library(huge)
#' data.npn<- huge.npn(alsData$exprs)
#' dim(data.npn) #160 17695
#'
#' # Extract KEGG interactome (max component)
#' KEGG<- properties(kegg)[[1]]
#' summary(KEGG)
#'
#' # KEGG modules with ALS perturbed edges using fast gready clustering
#' gD<- SEMdci(KEGG, data.npn, alsData$group, type="fgc")
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
        gD <- graph.union(gL)
    }
    else if (type == "none") {
        dest <- quiet(SEMggm2(g, data, group)$dest)
        dsub <- subset(dest, p.adjust(dest$pvalue, method = method) < alpha)
        ftm <- data.frame(from = dsub$rhs, to = dsub$lhs)
        gD <- graph_from_data_frame(ftm)
    }
    return(gD)
}

#' @title Graph properties summary and graph decomposition
#'
#' @description Produces a summary of network properties and returns
#' graph components (ordered by decreasing size), without self-loops.
#' @param graph Input network as an igraph object.
#' @param data An optional data matrix (default data = NULL) whith rows
#' corresponding to subjects, and columns to graph nodes (variables).
#' Nodes will be mapped onto variable names. 
#' @param ... Currently ignored.
#'
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
#' summary(g)
#' properties(g)
#'
properties<- function (graph, data = NULL, ...) 
{
    if (!is_igraph(graph)){
     ig <- graph_from_graphnel(graph)
    }else{ ig <- graph }
    if (!is.null(data)){
     nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
     ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
    }
    
	ig <- simplify(ig, remove.loops = TRUE)
    gcs <- igraph::decompose.graph(ig, min.vertices = 2)
    vsize <- sapply(1:length(gcs), function(x) vcount(gcs[[x]]))
	names(vsize) <- 1:length(vsize)
	gcs <- gcs[as.numeric(names(sort(vsize, decreasing = TRUE)))]
    
	cat("Frequency distribution of graph components\n\n")
    tab <- table(sapply(gcs, vcount))
    tab <- data.frame(n.nodes = as.numeric(names(tab)), n.graphs = as.numeric(tab))
    print(tab)
    cat("\nPercent of vertices in the giant component:", 
        round(100 * vcount(gcs[[1]])/vcount(ig), 1), "%\n\n")
    print(c(is.simple = is_simple(gcs[[1]]), is.dag = is_dag(gcs[[1]]), 
        is.directed = is_directed(gcs[[1]]), is.weighted = is_weighted(gcs[[1]])))
    cat("\n")
    print(c(which.mutual = table(which_mutual(gcs[[1]]))))
    
	return(gcs)
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
	if(is_igraph(graph)){
	 g <- as_graphnel(graph)
	} else {
	 g <- graph
	 graph <- graph_from_graphnel(graph)
	}
	
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
	
	g<- Rgraphviz::layoutGraph(g, layoutType=l, edgeAttrs=list(label=elab))
	graph::nodeRenderInfo(g)<- list(col = color, fill = vcol, label = vlab,
	                                lwd = lwd, lty = lty,
	                                textCol = color.txt,
	                                fontsize = fontsize, cex = cex,
	                                shape = vshape, width = w, height = h)
	graph::edgeRenderInfo(g) <- list(col = ecol, lty = lty, lwd = elwd)
	graph::graphRenderInfo(g) <- list(main = main, cex.main = cex.main,
	                                  font.main = font.main)
	Rgraphviz::renderGraph(g)
	
	return(invisible(g))
}

#' @title lavaan model to graph
#'
#' @description Convert a model, specified using lavaan syntax,
#' to an igraph object.
#' @param model Model specified using lavaan syntax.
#' @param directed Logical value. If TRUE (default), edge directions from
#' the model will be preserved. If FALSE, the resulting graph will
#' be undirected.
#' @param psi Logical value. If TRUE (default) covariances will be converted
#' into bidirected graph edges. If FALSE, covariances will be excluded from
#' the output graph.
#' @param verbose Logical value. If TRUE, a plot of the output graph will
#' be generated. For large graphs, this could significantly increase
#' computation time. If FALSE (default), graph plotting will be disabled.
#' @param ... Currently ignored.
#'
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
#' model<-"
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
#'
#' #(co)variances
#' PKA ~~ PIP3
#' "
#'
#' # Graph with covariances
#' G0 <- lavaan2graph(model, psi = TRUE)
#' plot(G0, layout = layout.circle)
#'
#' # Graph without covariances
#' G1 <- lavaan2graph(model, psi = FALSE)
#' plot(G1, layout = layout.circle)
#'
lavaan2graph<- function (model, directed = TRUE, psi = TRUE, verbose = FALSE, ...)			
{
    lav <- lavParTable(model, fixed.x = FALSE)
    lavb <- subset(lav, lav$op == "~")
    lavc <- subset(lav, lav$op == "~~" & (lav$rhs != lav$lhs))
	lavc<- lavc[lavc$user == 1,]
    ftm <- data.frame(cbind(from = lavb$rhs, to = lavb$lhs, label = lavb$label), 
        color = "black")
    if (nrow(lavc) != 0 & psi == TRUE) {
		ftmc1 <- data.frame(cbind(from = lavc$rhs, to = lavc$lhs, 
            label = "", color = "red"))
        ftmc2 <- data.frame(cbind(from = lavc$lhs, to = lavc$rhs, 
            label = "", color = "red"))
        ftm <- rbind(ftm, ftmc1, ftmc2)
    }
    graph <- graph_from_data_frame(ftm, directed = directed)
	if (verbose) gplot(graph)
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
	ftm <- igraph::as_data_frame(ig)

	if (is.directed(ig) & sum(which_mutual(ig)) > 0) {
		sel <- as.numeric(c(E(ig)[which_mutual(ig)]))
		ftm <- igraph::as_data_frame(ig)[-sel,]
		uft <- igraph::as_data_frame(ig)[sel,]
		ubg <- as.undirected(graph_from_data_frame(uft))
		ftb <- igraph::as_data_frame(ubg)
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
#' @param graph A graph as an igraph or as an adjacency matrix.
#' @param graphType character, is one of "dag" (default)' or "pdag".
#' DAG can contain the directed (->) and bi-directed (<->) edges,
#' while PDAG can contain the edges: ->, <->, and the undirected edges
#' (--) that represent edges whose direction is not known.
#' @param verbose A logical value. If TRUE, the output graph is shown.
#' This argument is FALSE by default.
#' @param ... Currently ignored.
#'
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
graph2dagitty<- function (graph, graphType = "dag", verbose = FALSE, ...) 
{
    if (!is.igraph(graph)) graph <- graph_from_adjacency_matrix(graph)
	dg <- graph - E(graph)[which_mutual(graph)]
    ug <- as.undirected(graph - E(graph)[!which_mutual(graph)])
    ed <- attr(E(dg), "vnames")
    eb <- attr(E(ug), "vnames")
    de <- paste(gsub("\\|", "->", ed), collapse = "\n")
    if (length(eb) == 0) {
        dagi <- paste0("dag {\n", de, "\n}")
    }
    else {
		be <- paste(gsub("\\|", "<->", eb), collapse = "\n")
        dagi <- paste0("dag {\n", de, "\n", be, "\n}")
    }
    if (graphType == "pdag"){
		dagi <- gsub("<->", "--", dagi)
		dagi <- gsub("dag", "pdag", dagi) #cat(dagy)
	}
	if (verbose) plot(dagitty::graphLayout(dagi))
  	return(dagi)
}

#' @title Graph conversion from dagitty to igraph
#'
#' @description Convert a dagitty object to a igraph object.
#' @param dagi A graph as a dagitty object ("dag" or "pdag").
#' @param verbose A logical value. If TRUE, the output graph is shown.
#' This argument is FALSE by default.
#' @param ... Currently ignored.
#' 
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Conversion from igraph to dagitty  (and viceversa)
#' dagi <- graph2dagitty(sachs$graph, verbose = TRUE)
#' graph <- dagitty2graph(dagi, verbose = TRUE)
#'
#' @return An igraph object.
#'
dagitty2graph<- function(dagi, verbose = FALSE, ...) 
{
	# edges to ftm
	edges<- dagitty::edges(dagi)
	dsel<- which(edges$e == "->")
	d1<- data.frame(from=edges$v[dsel], to=edges$w[dsel])
	bsel<- which(edges$e == "<->" | edges$e == "--")
	b1<- data.frame(from=edges$v[bsel], to=edges$w[bsel])
	b2<- data.frame(from=edges$w[bsel], to=edges$v[bsel])
	# ftm to graph
	ftm<- rbind(d1,b1,b2)
	graph<- graph_from_data_frame(ftm)
	if (verbose) gplot(graph)
	return(graph)
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
graph2dag<- function(graph, data, bap = FALSE, time.limit = Inf, ...)
{
	if (is_dag(graph)) return(dag=graph)
	# graph weighting by edge pvalues (r2z)
	graph <- weightGraph(graph, data)
	E(graph)$weight <- 1/(-log(E(graph )$pv))
	ftm <- igraph::as_data_frame(graph)
	wE <- ftm$weight
	names(wE) <- paste0(ftm[,1],":",ftm[,2])
	# delete all mutual edges <-> , i.e. <- & -> 
	ig <- graph - E(graph)[which_mutual(graph )]
	if (is_dag(ig) & bap == FALSE) {
	 cat("DAG conversion : TRUE\n")
	 return(dag = ig)
	}
	if (is_dag(ig) & bap == TRUE) {
	 cat("BAP conversion : TRUE\n")
	 return(bap = graph)
	}
	
	# subgraph isomorphism algorithm to detect all cycles of a given length
	# time limit: CPU time for the computation, in seconds (defaults=Inf)
	find.cycles <- function(graph, k, time.limit) {
	 ring <- graph.ring(k, TRUE)
	 subgraph_isomorphisms(ring, graph, "lad", time.limit=time.limit)
	}
	# function that identifies the right subisomorphisms to keep
	subisomorphism_rm_permutation <- function(si) {
	 is_first_min <- function(x) { return(x[1] == min(x)) }
	 sel <- lapply(si, is_first_min)
	 return(si[unlist(sel)])
	}
	# function that search max(weight) edge on each cycle
	max_edges <- function(x){
	 Ec <- vector()
	 for (i in 1:(length(x)-1)) Ec<- c(Ec, paste0(x[i],":",x[i+1]))
	 Ew <- wE[which(names(wE) %in% Ec)]
	 return(unlist(strsplit(names(Ew)[which(Ew == max(Ew))],":")))
	}

	for (k in 3:vcount(ig)){ #k=6
	 # find all cycles with k vertices(edges)
	 l <- find.cycles(ig, k, time.limit)
	 # remove permutations
	 l <- subisomorphism_rm_permutation(si=l)
	 # extract the vertices in each cycle
	 if (length(l) == 0) next
	 cycles <- lapply(1:length(l), function(x) names(l[[x]]))
	 # edges with max(weight)
	 l <- unique(lapply(cycles, max_edges))
	 E1 <- unlist(lapply(l, function(x) paste0(x[1],"|",x[2])))
	 E0 <- attr(E(ig), "vnames")
	 ig <- delete_edges(ig, which(E0 %in% E1))
	}
	if (bap == FALSE) {
	 cat("DAG conversion :", is_dag(ig),"\n")
	 return(dag = ig)
	}
    if (bap == TRUE) {
	 cat("BAP conversion :", is_dag(ig),"\n")
	 # add all mutual edges <-> , i.e. <- & ->
	 e <- as_edgelist(graph-E(graph)[!which_mutual(graph)])
	 return(bap = add_edges(ig, as.vector(t(e))))
	}
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

#' @title Assign edge orientation of an undirected graph
#'
#' @description Assign edge orientation of an undirected graph
#' through a given reference directed graph. The vertex (color)
#' and edge (color, width and weight) attributes of the input
#' undirected graph are preserved in the output directed graph. 
#'
#' @param ug An undirected graph as an igraph object.
#' @param dg A directed reference graph.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A directed graph as an igraph object.
#'
#' @examples
#'
#' # Graphs definition
#' G0 <- as.undirected(sachs$graph)
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
	 return(message(" ERROR: the input graph is a Directed graph !"))
	}
	if (!is_directed(dg)){
	 return(message(" ERROR: the reference graph is an UNdirected graph !"))
	}
	
	# graph edge comparison:
	mg <- as.directed(ug, mode = "mutual")
	exy0 <- attr(E(mg), "vnames")
	exy1 <- attr(E(dg)[which_mutual(dg) == FALSE], "vnames")
	exy2 <- exy0[which(exy0 %in% exy1)]
	if (length(exy2) == 0) return(graph = mg)
	str2 <- strsplit(exy2,"\\|")
	ftm1 <- matrix(unlist(str2),nrow=length(str2),byrow=TRUE)
	g1 <- graph_from_edgelist(ftm1, directed=TRUE)
	# plot(g1, layout=layout.circle)
	ug0 <- difference(ug, as.undirected(g1))
	mg0 <- as.directed(ug0, mode = "mutual")
	
	#output graph with attr matching:
	g <- graph.union(mg0, g1)
	gattr<- attrMatch(mg, g)
	V(g)$color <- gattr[[1]]
	E(g)$color <- gattr[[2]]
	E(g)$width <- gattr[[3]]
	E(g)$weight <- gattr[[4]]
	# plot(g, layout=layout.circle)
	
	return(graph = g)
}

attrMatch<- function(g1, g2, ...)
{
	if (!is.null(V(g1)$color)) {
	 idx<- match(V(g2)$name, V(g1)$name)
	 Vcol<- V(g1)$color[idx]
	}else{
	 Vcol<- rep(NA, vcount(g2))
	}
	if (!is.null(E(g1)$color)) {
	 idx<- match(attr(E(g2), "vnames"), attr(E(g1), "vnames"))
	 Ecol<- E(g1)$color[idx]
	}else{
	 Ecol<- rep("gray", ecount(g2))
	}
	if (!is.null(E(g1)$width)) {
	 idx<- match(attr(E(g2), "vnames"), attr(E(g1), "vnames"))
	 Ewid<- E(g1)$width[idx]
	}else{
	 Ewid<- rep(1, ecount(g2))
	}
	if (!is.null(E(g1)$weight)) {
	 idx<- match(attr(E(g2), "vnames"), attr(E(g1), "vnames"))
	 Ewei<- E(g1)$weight[idx]
	}else{
	 Ewei<- rep(1, ecount(g2))
	}
	return(list(Vcol, Ecol, Ewid, Ewei))
}

#' @title Vertex and edge graph coloring on the base of fitting
#'
#' @description Add vertex and edge color attributes to an igraph object,
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
#' est20 <- subset(parameterEstimates(sem2$fit), group == 1)[, -c(4, 5)]
#' est21 <- subset(parameterEstimates(sem2$fit), group == 2)[, -c(4, 5)]
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
			G <- est
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
    if (is.null(y)) {
		vnames <- colnames(x)
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
		vnames <- intersect(colnames(x), colnames(y))
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

#' @title SEM-based mediation analysis
#'
#' @description Perform mediation analysis through structural equation 
#' modeling (SEM). The SEM is represented by two equations: 
#' x2 ~ b*x1 and y ~ a*x1 + c*x2, where x1 is the exposure factor, 
#' y is the outcome variable, and x2 is the mediator.
#' This model is iteratively applied for each column in the input 
#' data matrix, given x1 and y. A new parameter IE := b*c is defined, 
#' corresponding to the mediation effect.
#' 
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to potential mediatiors.
#' @param x1 A numeric vector corresponding to the exposure factor.
#' @param y A numeric vector corresponding to the outcome variable.
#' @param model The user might specify a different model. To use the 
#' given mediation model, leave this argument to NULL (default).
#' @param group A vector to group subjects for multigroup analysis. 
#' By default, \code{group = NULL}.
#' @param fixed.x Logical value. If TRUE (default), the exposure factor 
#' (i.e., the exogenous variable) estimates are fixed to its sampling 
#' values. If FALSE, their estimates are considered as free parameters.
#' @param se Standard error estimation method. Use "standard" (default) 
#' for conventional standard errors, based on the inverse information 
#' matrix. Use "robust.huber.white" for pseudo ML Huber-White estimation.
#' @param std.x Logical value. If TRUE, the exposure factor is standardized 
#' (default = FALSE).
#' @param std.y Logical value. If TRUE, the outcome variable is standardized 
#' (default = FALSE).
#' @param std.m Logical value. If TRUE, the mediator is standardized 
#' (default = FALSE).
#' @param adjust.method One of the multiple testing correction methods 
#' from the \code{p.adjust} function. By default, the Benjamini-Hochberg 
#' ("BH") adjuctment method is used.
#' @param gene.symbols Logical value. If variables are genes and data 
#' matrix column names are entrez gene IDs, this argument can be set to 
#' TRUE to obtain a new column with official gene symbols 
#' (default = FALSE).
#' @param ... Currently ignored.
#' 
#' @export
#'
#' @return A list of two objects:
#' \enumerate{
#' \item "model", The output model as an igraph object;
#' \item "clusters", A list containing model clustering results.
#' It includes tree more objects:
#' \itemize{
#' \item "membership", A vector reporting cluster membership for each node;
#' \item "fit", hidden model fitting as a lavaan object;
#' \item "dataHM", data matrix containing cluster scores (hidden variables) 
#' alongside the original data variables.
#' }
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' 
#' library(huge)
#' library(org.Hs.eg.db)
#' 
#' als.npn <- huge.npn(alsData$exprs)
#'
#' # Random sampling 10 variables (possible mediators)
#' j <- sample.int(ncol(als.npn), 10)
#' 
#' # Generate a random binary exposure vector
#' exposure <- rbinom(n = nrow(als.npn), size = 1, prob = 0.3)
#' 
#' # Generate a random gaussian outcome vector
#' outcome <- rnorm(n = nrow(als.npn))
#' 
#' mediators <- sem.mediation(data = als.npn[, j],
#'                            x1 = exposure,
#'                            y = outcome,
#'                            gene.symbols = TRUE)
#' 
#' }
#'
sem.mediation <- function(data, x1, y, model = NULL, group = NULL,
                          fixed.x = TRUE, se = "robust.huber.white",
                          std.x = FALSE, std.y = FALSE, std.m = FALSE,
                          adjust.method = "BH", gene.symbols = FALSE,
                          ...) {
  if (is.null(model)) {
    model <- '
    x2 ~ b*x1
    y ~ a*x1 + c*x2
    IE := b*c
    TE := a + b*c'
  }
  
  n <- nrow(data)
  m <- ncol(data)
  mediation <- data.frame(ID = character(), r2 = numeric(), est = numeric(),
                          se = numeric(), z = numeric(),
                          pvalue = numeric())
  
  if (std.x) x1 <- scale(x1)
  if (std.y) y <- scale(y)
  
  for (i in 1:m) {
    x2 <- data[, i]
    if (std.m) x2 <- scale(x2)
    semdata <- data.frame(x1, x2, y)
    mest <- fit.mediation(model = model, data = semdata, group = group,
                          fixed.x = fixed.x, se = se)
    mest$ID <- colnames(data)[i]
    mediation <- rbind(mediation, mest[c(6, 1:5)])
  }
  
  mediation$FDR <- p.adjust(mediation$pvalue, method = adjust.method)
  if (gene.symbols) {
    mediation$symbol <- mapIds(org.Hs.eg.db, mediation$ID,
                               column = 'SYMBOL',
                               keytype = 'ENTREZID')
    mediation <- mediation[, c(1, 8, 2:7)]
  }
  mediation <- mediation[order(mediation$FDR),]
  return(mediation)
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
#'
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
  as.numeric(sort(graph.bfs(g, nodes, mode = "in", unreachable = F)$order,
                  na.last = NA))
}

#' @rdname ancestry
#' @export
descendants <- function(g, nodes)
{
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, mode = "out", unreachable = F)$order,
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

quiet <- function(x) {
	sink(tempfile())
	on.exit(sink())
	invisible(force(x))
}
