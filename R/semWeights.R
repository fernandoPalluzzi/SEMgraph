#  SEMgraph library
#  Copyright (C) 2019 Fernando Palluzzi; Mario Grassi
#  e-mail: <fernando.palluzzi@gmail.com>
#  University of Pavia, Department of Brain and Behavioral Sciences
#  Via Bassi 21, Pavia, 27100 Italy

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


#' @title Graph weighting via moderated mediation trivariate model
#'
#' @description Weight connections of a directed graph by using two-group
#' SEM fittinng. For each directed edge j -> k, the SEM tests group effect 
#' on the source node j, the target node k, and their directed 
#' interaction. The weight parameter w is thus defined as the degree 
#' centrality-scaled linear combination of the three effects. The actual 
#' edge weight will be the negative natural logarithm of the z-test = w/SE(w) 
#' P-value. Node (source and target) degree centrality scaling contrasts 
#' the effect of weights inflation towards highly connected genes (i.e., hubs). 
#' This comes from the evidence that perturbed systems often involve essential 
#' nodes with low degree (e.g., community-bridging "bottleneck" connectors) or 
#' peripheral nodes (e.g., effectors).
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#'
#' @return A weighted graph, as an igraph object, with two edge attributes:
#' zsign and pvalue. Attribute zsign may have values [-1, 0, +1], where 0
#' is assigned if Pr(|z|) > 0.1, while the sign of the z-test (either
#' positive or negative) is assigned when Pr(|z|) < 0.1. Attribute pvalue
#' correspond to the z-test P-value.
#'
#' @import igraph
#' @import lavaan
#' @export
#'
#' @references
#'
#' Palluzzi F, Ferrari R, Graziano F, Novelli V, Rossi G, Galimberti D,
#' Rainero I, Benussi L, Nacmias B, Bruni AC, Cusi D, Salvi E, Borroni B,
#' Grassi M. (2017). A novel network analysis approach reveals DNA damage,
#' oxidative stress and calcium/cAMP homeostasis-associated biomarkers
#' in frontotemporal dementia. PLoS ONE 12(10): e0185797.
#' https://doi.org/10.1371/journal.pone.0185797
#'
#' @seealso \code{\link{edgeweight.cfa}} and \code{\link{edgeweight.r2z}}
#' for weighting bidirected (or undirected) graphs
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' G <- edgeweight.sem(graph = graph, data = data, group = group)
#' E(G)$pv
#' E(G)$zsign
#'
edgeweight.cov <- function(graph, data, group)
{
	# Set genes and from-to-matrix representation of gene-gene links
	genes <- colnames(data)
	graph <- induced_subgraph(graph, vids = which(V(graph)$name %in% genes))
	degree <- igraph::degree(graph, v = V(graph), mode = "all")
	ftm <- as_data_frame(graph)
	Y <- scale(data)
	est <- NULL
	
	for (i in 1:nrow(ftm)) {
		cat("\r", "edge weight", i, "of", nrow(ftm))
		flush.console()
		YC <- data.frame(cbind(Y[, c(ftm[i, 1], ftm[i, 2])], group))
		colnames(YC)[1:2] <- c("x", "y")
		
		dx<- degree[ftm[i,1]]
		dy<- degree[ftm[i,2]]
		model<- paste0(
			'x ~ c(a1,a2)*1
			y ~ c(b1,b2)*1
			x ~~ c(c1,c2)*y 
			#w:= (a1-a2)/2+(b1-b2)/2+(c1-c2)
			w:= (a1-a2)/', dx, '+(b1-b2)/', dy, '+(c1-c2)')
		try(fit<- sem(model, data=YC, group="group", fixed.x=FALSE))
		try(est<- c(est, list(parameterEstimates(fit))))
	}
	cat("\n done.\n")
	
	zeta <- function(x) est[[x]]$z[est[[x]]$label == "w"]
	pvalue <- function(x) est[[x]]$pvalue[est[[x]]$label == "w"]
	B <- sapply(1:length(est), zeta)
	zsign <- ifelse(abs(B) < 1.64, 0, sign(B))
	pv <- sapply(1:length(est), pvalue)
	pv[is.na(pv)] <- 0.5
	pv[pv == 0] <- 1*10^-9
	pv[pv == 1] <- 1-1*10^-9
	ftm <- cbind(ftm, zsign, pv)
	gdf <- graph_from_data_frame(ftm, directed = is.directed(graph))
	
	return(graph = gdf)
}

#' @title Graph weighting via single-factor Confirmatory Factor Analysis (CFA)
#'
#' @description Weight connections of a bidirected (covariance or
#' correlation) graph using a single-factor CFA model in which the group
#' variable has a direct effect on a latent variable acting over pairs of
#' nodes. This function can be used to model correlation between nodes.
#' Covariance between pairs of nodes is absorbed by a latent variable
#' representing unknown common cause acting on both observed variables
#' (nodes). Edge weights correspond to the z-test (= estimate/standardError)
#' sign and P-value of the group effect over the latent variable.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#'
#' @return A weighted graph, as an igraph object, with two edge attributes:
#' "zsign" and "pv". Attribute "zsign" may have values [-1, 0, +1],
#' where 0 is assigned if Pr(|z|) > 0.1, while the sign of the z-test
#' (either positive or negative) is assigned when Pr(|z|) < 0.1. Attribute
#' "pv" correspond to the z-test P-value.
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cov
#' @export
#'
#' @references
#'
#' Palluzzi F, Ferrari R, Graziano F, Novelli V, Rossi G, Galimberti D,
#' Rainero I, Benussi L, Nacmias B, Bruni AC, Cusi D, Salvi E, Borroni B,
#' Grassi M. (2017). A novel network analysis approach reveals DNA damage,
#' oxidative stress and calcium/cAMP homeostasis-associated biomarkers
#' in frontotemporal dementia. PLoS ONE 12(10): e0185797.
#' https://doi.org/10.1371/journal.pone.0185797
#'
#' @seealso \code{\link{edgeweight.sem}} for weighting directed graphs,
#' and \code{\link[lavaan]{cfa}} for more information about CFA
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' G <- edgeweight.cfa(graph = graph, data = data, group = group)
#' E(G)$pv
#' E(G)$zsign
#'
edgeweight.cfa <- function(graph, data, group)
{
	# Set genes and from-to-matrix representation of gene-gene links
	genes <- colnames(data)
	graph <- induced_subgraph(graph, vids = which(V(graph)$name %in% genes))
	ftm <- as_data_frame(graph)
	Y <- scale(data)
	est <- NULL 
	
	for (i in 1:nrow(ftm)) {
		cat("\r", "edge weight", i, "of", nrow(ftm))
		flush.console()
		YC <- data.frame(cbind(Y[, c(ftm[i, 1], ftm[i, 2])], group))
		colnames(YC)[1:2] <- paste0("y", 1:2)
		if(cov(YC$y1, YC$y2) < 0) YC$y1 <- -1*YC$y1
		a <- sqrt(cov(YC$y1, YC$y2))
		
		model <- paste0(
			'f =~ ',a,'*y1+',a,'*y2
			f ~ group
			y1~~',1-a^2,'*y1
			y2~~',1-a^2,'*y2')
		
		try(fit <- cfa(model, data = YC, fixed.x = TRUE))
		try(est <- c(est, list(parameterEstimates(fit))))
		# summary(fit); inspect(fit,"est") 
	}
	cat("\n done.\n")
	names(est) <- paste0(ftm[, 1], "<->", ftm[, 2])
	
	zeta <- function(x) est[[x]]$z[est[[x]]$op == "~"]
	pvalue <- function(x) est[[x]]$pvalue[est[[x]]$op == "~"]
	B <- sapply(1:length(est), zeta)
	zsign <- ifelse(abs(B) < 1.64, 0, sign(B))
	var <- sapply(1:length(est), function(x) est[[x]]$est[6])
	pv <- sapply(1:length(est), pvalue)
	pv[is.na(pv) | var < 0] <- 0.5
	pv[pv == 0] <- 1*10^-9
	pv[pv == 1] <- 1-1*10^-9
	ftm <- cbind(ftm, zsign, pv)
	gdf <- graph_from_data_frame(ftm, directed = is.directed(graph))
	
	return(graph = gdf)
}

#' @title Graph weighting via Fisher's r-to-z transformation of correlation
#' coefficient
#'
#' @description Weight connections of a bidirected (covariance or correlation)
#' graph using two-group Fisher's transform correlation test, using standard
#' error correction according to Yuan et al. 2013. Edge weights correspond
#' to the sign and P-value of the correlation group difference test.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#'
#' @return A weighted graph, as an igraph object, with two edge attributes:
#' "zsign" and "pv". Attribute "zsign" may have values [-1, 0, +1], where
#' 0 is assigned if Pr(|z|) > 0.1, while the sign of the z-test (either
#' positive or negative) is assigned when Pr(|z|) < 0.1. Attribute "pv"
#' correspond to the z-test P-value.
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cor pnorm
#' @export
#'
#' @references
#' Yuan Z, Liu H, Zhang X, Li F, Zhao J, Zhang F, Xue F (2013).
#' From Interaction to Co-Association - A Fisher r-To-z Transformation-Based
#' Simple Statistic for Real World Genome-Wide Association Study.
#' PLoS One, 8(7): e70774. https://doi.org/10.1371/journal.pone.0070774
#'
#' @seealso \code{\link{edgeweight.sem}} for weighting directed graphs
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' G <- edgeweight.r2z(graph = graph, data = data, group = group)
#' E(G)$pv
#' E(G)$zsign
#'
edgeweight.r2z <- function(graph, data, group)
{
	genes <- colnames(data)
	#head(genes)
	graph <- induced_subgraph(graph, vids = which(V(graph)$name %in% genes))
	ftm <- as_data_frame(graph)

	n1 <- length(group[group == 1])
	n0 <- length(group[group == 0])
	n <- n1 + n0
	pv <- zsign <- vector()

	for(i in 1:nrow(ftm)){
		x <- data[, ftm[i, 1]]
		y <- data[, ftm[i, 2]]
		x1 <- data[group == 1, ftm[i, 1]]
		y1 <- data[group == 1, ftm[i, 2]]
		x0 <- data[group == 0, ftm[i, 1]]
		y0 <- data[group == 0, ftm[i, 2]]

		z <- 0.5*log((1 + cor(x, y))/(1 - cor(x, y)))
		z1 <- 0.5*log((1 + cor(x1, y1))/(1 - cor(x1, y1)))
		z0 <- 0.5*log((1 + cor(x0, y0))/(1 - cor(x0, y0)))
		fz <- exp(-0.176 + 0.09*z + 0.315*z^2)

		# u-test & p-value
		try(if(abs(z) <= 0.5) u <- (z1 - z0)/sqrt(1/(n1 - 3) + 1/(n0 - 3)))
		try(if(abs(z) > 0.5) u <- (z1 - z0)/(fz*sqrt(1/(n1 - 3) + 1/(n0 - 3))))
		try(pv[i] <- 2*pnorm(-abs(u)))
		#2*(1 - pnorm(abs(t)))
		try(zsign[i] <- ifelse(abs(u) < 1.64, 0, sign(u)))
	}

	pv[is.na(pv)] <- 0.5
	pv[pv == 0] <- 1*10^-10
	pv[pv == 1] <- 1 - 1*10^-10
	ftm <- cbind(ftm, zsign, pv)
	#head(ftm)
	gdf <- graph_from_data_frame(ftm, directed = is.directed(graph))

	return(graph = gdf)
}

#' @title Seed extraction methods
#'
#' @description Extract important nodes (i.e., seeds) from a network using
#' three methods: group influence p-value, Minimum Spanning Tree (MST),
#' and prototype node detection.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#' @param h A decimal numbed, corresponding to the prototype clustering
#' distance measure (= 1 - abs(correlation)) cutoff. By default, h = 0.2
#' (corresponding to an absolute correlation coefficient of 0.8).
#' @param alpha Significance level of the group effect over graph nodes.
#' By default, alpha = 0.05.
#'
#' @return The input graph endowed with three new binary (1: seed,
#' 0: not-seed) vertex attibutes:
#' \enumerate{
#' \item "spvlm", seed based on P-value of the simple linear regression 
#' model node ~ group
#' \item "sprot", prototype seeds derived from R method
#' \code{\link[protoclust]{protoclust}}
#' \item "smst", non-leaf nodes of the MST derived from the input network
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cor lm as.dist
#' @importFrom graphics abline
#' @importFrom protoclust protoclust protocut
#' @export
#'
#' @seealso \code{\link[stats]{lm}} for direct effect group estimation
#' on graph nodes; \code{\link[igraph]{mst}} using Prim's algorithm for
#' weighted networs; \code{\link[protoclust]{protoclust}} for prototype
#' nodes calculation
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' G <- seedweight(graph = graph, data = data, group = group)
#' table(V(G)$spvlm)
#' table(V(G)$sprot)
#' table(V(G)$smst)
#'
seedweight <- function(graph, data, group, h = 0.2, alpha = 0.05)
{
	# Set genes, graph and data objects
	genes <- colnames(data)
	#head(genes)
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% genes))
	Y <- data[, which(colnames(data) %in% V(ig)$name)]
	D <- as.dist(1 - abs(cor(Y)))

	# Seed nodes by p_values
	pv.lm <- function(x) { summary(lm(x~group))$coefficients[2, 4] }
	pvlm <- apply(Y, 2, pv.lm)
	#pvlm
	seed1 <- V(ig)$name[pvlm < alpha]
	V(ig)$spvlm <- ifelse(V(ig)$name %in% seed1, 1, 0)

	# Seed nodes by prototypes
	plot(hc <- protoclust(D))
	abline(h = 0.20, lty = 1, col = "red")
	# Cut distance threshold fixed to 0.2 (i.e., 0.8 correlation)
	cutd <- protocut(hc, h = h)
	seed2 <- hc$labels[cutd$protos]
	V(ig)$sprot <- ifelse(V(ig)$name %in% seed2, 1, 0)

	# Seed nodes by MST
	eweight <- rep(1, ecount(ig))
	mst <- minimum.spanning.tree(ig, weights = eweight, algorithm = "prim")
	seed3 <- V(mst)$name[which(igraph::degree(mst) > 1)]
	V(ig)$smst <- ifelse(V(ig)$name %in% seed3, 1, 0)

	return(graph = ig)
}
