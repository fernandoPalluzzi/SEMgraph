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

#' @title Graph weighting methods
#'
#' @description Add data-driven edge and node weights to the input graph.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group Binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects. By default, \code{group = NULL}. If group is not NULL, also
#' node weighting is actived, and node weights correspond to the sign
#' (-1 if z <= 2, +1 if z > 2, 0 otherwise) and P-value of the
#' z-test = b/SE(b) from simple linear regression y ~ x
#' (i.e., glm(node ~ group)).
#' @param method Edge weighting method. It can be one of the following:
#' \enumerate{
#' \item "r2z", weight edges of a graph using Fisher's r-to-z transform 
#' (Fisher, 1915) to test the correlation coefficient of pairs 
#' of interacting nodes, if \code{group == NULL}. Otherwise, the correlation
#' difference between group will be tested and edge weights correspond to
#' the sign (-1 if z <= 2, +1 if z > 2, 0 otherwise) and P-value of the group
#' correlation difference.
#' \item "sem", edge weights are defined by a SEM model that implies 
#' testing the group effect simultaneously on the source node and the sink mode.
#' A new parameter w is defined as the weighted sum of the total effect 
#' of the group on source and sink nodes, adjusted by node degree centrality, 
#' and edge weights correspond to the sign (-1 if z <= 2, +1 if z > 2, 
#' 0 otherwise) and P-value of the z-test = w/SE(w). 
#' Not available if \code{group == NULL}.
#' \item "cov", edge weights are defined by a new parameter w combining 
#' the group effect on the source node (mean group difference, adjusted 
#' by source degree centrality), the sink node (mean group difference, 
#' adjusted by sink degree centrality), and the source-sink interaction 
#' (correlation difference). Edge weights correspond to the sign
#' (-1 if z <= 2, +1 if z > 2, 0 otherwise) and P-value of the 
#' z-test = w/SE(w) of the combined difference of the group over 
#' source node, sink node, and their connection. 
#' Not available if \code{group == NULL}.
#' \item "cfa", edge weights are defined by a CFA1 model that implies 
#' testing the group effect, w on a latent variable (LV) with observed 
#' indicators two interacting nodes, fixing loading coefficients and residual
#' variances for model identification. Edge weights correspond to the sign
#' (-1 if z <= 2, +1 if z > 2, 0 otherwise) and P-value of the z-test = w/SE(w)
#' of the group effect on the LV. Not available if \code{group == NULL}.
#' }
#' @param limit An integer value corresponding to the number of graph 
#' edges. Beyond this limit, multicore computation is enabled to reduce 
#' the computational burden. 
#' By default, \code{limit = 10000}.
#' @param ... Currently ignored.
#'
#' @return A weighted graph, as an igraph object.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 2022;, btac567, https://doi.org/10.1093/bioinformatics/btac567
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
#' # New nodes attributes
#' head(V(G)$pv); summary(V(G)$pv)
#' head(V(G)$zsign); table(V(G)$zsign)
#' 
weightGraph<- function(graph, data, group = NULL, method = "r2z", limit = 10000, ...) 
{
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	ig <- quiet(properties(ig)[[1]])
	degree <- igraph::degree(ig, v = V(ig), mode = "all")
	ftm <- as_data_frame(ig)
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
	gdf <- graph_from_data_frame(ftm, directed = is.directed(graph))
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

#' @title Active module identification
#'
#' @description Uses different information flow and tree-based strategies 
#' for identifying active modules (e.g., disease modules), including a 
#' perturbed subset of nodes and edges. 
#' Function scalability enables graph reduction at both pathway and 
#' entire interactome scales.
#' @param graph An igraph object.
#' @param type Module identification method. If \code{type = "kou"}, 
#' the Steiner tree algorithm will be applied. 
#' If \code{type = "usp"}, the resulting graph will be the union of all 
#' significant shortest paths. If \code{type = "rwr"}, the random walk 
#' with restart algorithm will be enabled. Finally, if \code{type = "hdi"}, 
#' the heat diffusion algorithm is used.
#' @param seed Either a user-defined vector containing seed node names 
#' or one among: "pvlm", "proto", or "qi", corresponding to the seed 
#' name attribute yielded by \code{\link[SEMgraph]{weightGraph}}.
#' @param eweight Edge weight type derived from
#' \code{\link[SEMgraph]{weightGraph}} or from user-defined distances. 
#' This option determines the weight-to-distance transform. If set to 
#' "none" (default), edge weights will be set to 1. 
#' If \code{eweight = "kegg"}, repressing interactions (-1) will be set 
#' to 1 (maximum distance), neutral interactions (0) will be set to 0.5, 
#' and activating interactions (+1) will be set to 0 (minimum distance).
#' If \code{eweight = "zsign"}, all significant interactions will be set 
#' to 0 (minimum distance), while non-significant ones will be set to 1.
#' If \code{eweight = "pvalue"}, weights (p-values) will be transformed 
#' to the inverse of negative base-10 logarithm. 
#' If \code{eweight = "custom"}, the algorithm will use the distance 
#' measure specified by the user as "weight" edge attribute.
#' @param alpha Significance level to assess shortest paths significance, 
#' when type is "usp". By default, \code{alpha = 0.05}.
#' @param top Number of top nodes for the "rwr" and "hdi" algorithms. The 
#' output subgraph is induced by the top-n ranking nodes.
#' By default, \code{top = 100} (i.e., the top-100 of nodes are selected).
#' @param limit An integer value corresponding to the number of graph 
#' edges. If \code{type = "usp"}, beyond this limit, multicore computation 
#' is enabled to reduce the computational burden. 
#' By default, \code{limit = 10000}.
#' @param ... Currently ignored.
#'
#' @details Graph filtering algorithms include:
#' \enumerate{
#' \item "kou", the Steiner tree connecting a set of seed nodes, using 
#' the algorithm suggested by Kou et al. (1981);
#' \item "usp", generates a subnetwork as the union of the significant 
#' (P-value < alpha) shortest paths between the seeds set;
#' \item "rwr", Random Walk with Restart, a wrapper for \code{random.walk}
#' function of the R package \code{diffusr};
#' \item "hdi", Heat Diffusion algorithm, a wrapper for \code{heat.diffusion}
#' function of the R package \code{diffusr}.
#' }
#'
#' @return An active module, an igraph object with colored nodes
#' (seed = "green", and connector = "white").
#'
#' @import igraph
#' @importFrom diffusr random.walk heat.diffusion
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' 
#' Palluzzi F, Grassi M (2021). SEMgraph: An R Package for Causal Network 
#' Analysis of High-Throughput Data with Structural Equation Models. 
#' <arXiv:2103.08332>
#' 
#' Kou L, Markowsky G, Berman L (1981). A fast algorithm for Steiner trees.
#' Acta Informatica, 15(2): 141-145. <https://doi.org/10.1007/BF00288961>
#'
#' Simon Dirmeier (2018). diffusr: Network Diffusion Algorithms. R
#' package version 0.1.4.
#' <https://CRAN.R-project.org/package=diffusr/>
#'
#' @examples
#' 
#' # Graph weighting
#' G <- weightGraph(graph = sachs$graph, data = sachs$pkc, group = sachs$group,
#'                  method = "r2z",
#'                  seed = c(0.05, 0.5, 0.5))
#' 
#' # RWR algorithm, seeds and edge P-values as weights
#' R1 <- activeModule(graph = G, type = "kou", seed = "pvlm", eweight = "pvalue")
#' R2 <- activeModule(graph = G, type = "kou", seed = "proto", eweight = "pvalue")
#' R3 <- activeModule(graph = G, type = "kou", seed = "qi", eweight = "pvalue")
#' 
#' # Graphs
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mar=rep(2, 4))
#' plot(G, layout = layout.circle, main = "input graph")
#' box(col = "gray")
#' plot(R1, layout = layout.circle, main = "lm P-value (alpha = 0.05)")
#' box(col = "gray")
#' plot(R2, layout = layout.circle, main = "prototype (h = 0.5)")
#' box(col = "gray")
#' plot(R3, layout = layout.circle, main = "closeness (q = 0.5)")
#' box(col = "gray")
#' par(old.par)
#' 
activeModule <- function(graph, type, seed, eweight = "none", alpha = 0.05,
                         top = 100, limit = 10000, ...)
{
	if (length(eweight) == 1) {
	 if (eweight == "kegg") eweight <- (1 - E(graph)$weight)/2
	 else if (eweight == "zsign") eweight <- 1 - abs(E(graph)$zsign)
	 else if (eweight == "pvalue") eweight <- 1/(-log10(E(graph)$pv))
	 else if (eweight == "custom") eweight <- E(graph)$weight
	 else if (eweight == "none") eweight <- rep(1, ecount(graph))
	}
	
	if (length(seed) == 1) {
	 if (seed == "pvlm") seed <- V(graph)$name[V(graph)$pvlm == 1]
	 else if (seed == "proto") seed <- V(graph)$name[V(graph)$proto == 1]
	 else if (seed == "qi") seed <- V(graph)$name[V(graph)$qi == 1]
	} else {
	 seed <- seed[seed %in% V(graph)$name]
	}
	
	if (type == "kou" & length(seed) != 0) {
	 R <- SteinerTree(graph, seed = seed, eweight = eweight)
	
	} else if (type == "usp" & length(seed) != 0) {
	 R <- USPG(graph, seed = seed, eweight = eweight, alpha = alpha, limit = limit)
	
	} else if (type == "rwr" & length(seed) != 0) {
	 R <- RWR(graph, seed = seed, eweight = eweight, algo = "rwr", top = top)
	
	} else if (type == "hdi" & length(seed) != 0) {
	 R <- RWR(graph, seed = seed, eweight = eweight, algo = "hdi", top = top)
	}
	V(R)$color <- ifelse(V(R)$name %in% seed, "green", "white")
	
	return(R)
}

seedweight <- function(ig, data, group, alpha, h, q, ...)
{
	# Set data object
	Y <- data[, which(colnames(data) %in% V(ig)$name)]
	D <- as.dist(1 - abs(cor(Y)))
	
	# Seed nodes by p_values
	pv.lm <- function(x) { summary(lm(x~group))$coefficients[2, 4] }
	if (!is.null(group)) {
		pvlm <- apply(Y, 2, pv.lm)
		p.adj <- p.adjust(pvlm, method="BH")
		seed1 <- V(ig)$name[p.adj < alpha]
		V(ig)$pvlm <- ifelse(V(ig)$name %in% seed1, 1, 0)
	} else {
		V(ig)$pvlm <- 0
	}
	
	# Seed nodes by prototypes
	#plot(hc <- protoclust::protoclust(D))
	hc <- protoclust::protoclust(D)
	#abline(h = h, lty = 1, col = "red")
	
	# Cut distance threshold fixed to 0.2 (i.e., 0.8 correlation)
	cutd <- protoclust::protocut(hc, h = h)
	seed2 <- hc$labels[cutd$protos]
	V(ig)$proto <- ifelse(V(ig)$name %in% seed2, 1, 0)

	# Seed nodes by closeness
	suppressWarnings(qi <- igraph::closeness(ig, mode = "all", weights = NA))
	seed3 <- V(ig)$name[which(qi > quantile(qi, probs = q))]
	V(ig)$qi <- ifelse(V(ig)$name %in% seed3, 1, 0)
	
	return (ig)
}

RWR <- function(graph, seed, eweight, algo, top, ...)
{
	E(graph)$weight <- eweight
	W <- as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)
	p0 <- ifelse(V(graph)$name %in% seed, 1, 0)
	q <- 1-top/vcount(graph)
	
	if (algo == "rwr") {
		pt <- diffusr::random.walk(p0 = p0, W, r = 0.5)
		score <- pt$p.inf
		#score <- rowSums(pt$transition.matrix)
		top <- rownames(W)[score > quantile(score, q)]
	} else {
		#ht <- heat.diffusion(h0 = pval, W, t = 0.5)
		ht <- diffusr::heat.diffusion(h0 = p0, W, t = 0.5)
		top <- rownames(W)[ht > quantile(ht, q)]
	}
	
	return(graph = induced_subgraph(graph, top))
}

SteinerTree <- function(graph, seed, eweight, ...)
{
	# Define graph, edge weights, and distance matrix
	E(graph)$weight <- eweight
	D <- igraph::distances(graph, v = seed, to = seed,
	                       mode = "all",
	                       weights = eweight)

	# Step 1: complete undirected distance graph Gd for terminal nodes
	Gd <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE)
	Gd <- Gd - igraph::edges(E(Gd)[which(E(Gd)$weight == Inf)])

	# Step 2: MST T1 of the complete distance graph Gd
	T1 <- minimum.spanning.tree(Gd, weights = NULL, algorithm = "prim")

	# Step 3: for each edge in T1, replace it with the shortest path in ig
	edge_list <- as_edgelist(T1)
	N <- nrow(edge_list)
	subgraph <- vector()

	for (n in 1:N) {
		i <- edge_list[n, 1]
		j <- edge_list[n, 2]

		# Extract from ig all nodes of the shortest paths between edges of T1
		path <- shortest_paths(graph, from = V(graph)[i], to = V(graph)[j],
		                       mode = "all",
		                       weights = eweight,
		                       output = "both")
		vpath <- V(graph)$name[path$vpath[[1]]]
		subgraph <- igraph::union(subgraph, vpath)
	}

	# Step 4: MST Ts of the extracted (induced) sub-graph Gs of ig
	Gs <- induced_subgraph(graph, unique(subgraph))
	Ts <- minimum.spanning.tree(Gs, weights = NULL, algorithm = "prim")

	# Step 5: Pruning non-seed genes with degree=1 (one at time) from Ts
	St <- Ts
	idx <- ifelse(V(St)$name %in% seed == TRUE, FALSE, TRUE)
	i <- 1
	I <- length(V(St)[idx]) + 1
	while(i < I) {
		K <- igraph::degree(St, v = V(St), mode = "all")
		todel <- names(which(K == 1))
		todel <- todel[which(!todel %in% seed)]
		if(length(todel) > 0) {
			St <- igraph::delete.vertices(St, todel)
		}
		i <- i + 1
	}
	
	return(St)
}

USPG <- function(graph, seed, eweight, alpha, limit, ...)
{
	# Define graph, edge weights, and distance matrix
	E(graph)$weight <- eweight
	if (is.null(E(graph)$pv)) E(graph)$pv <- rep(0, ecount(graph))
	D <- igraph::distances(graph, v = seed, to = seed, mode = "all",
	                       weights = eweight)
	
	# Complete directed distance graph for terminal nodes
	Gd <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE)
	Gd <- Gd - igraph::edges(E(Gd)[which(E(Gd)$weight == Inf)])
	#plot(as_graphnel(Gd)); Gd; E(Gd)$weight
	
	# For each edge in Gd, replace it with the shortest path
	ftm <- as_edgelist(Gd)
	N <- nrow(ftm)
	if (alpha == 1) alpha <- N
	
	local <- function(x) {
		i <- x[[1]]
		j <- x[[2]]
		
		# Extract nodes from the shortest paths between edges of Gd
		path <- shortest_paths(graph, from = V(graph)[i], to = V(graph)[j],
		                       mode = "all", weights = E(graph)$weight,
		                       output = "both")
		vpath <- V(graph)$name[path$vpath[[1]]]
		r <- length(vpath) - 1
		pvalue <- E(graph)$pv[path$epath[[1]]]
		
		# Fisher's combined significance test of the shortest path
		ppath <- 1 - pchisq(-2*sum(log(pvalue)), df = 2*length(pvalue))
		ppath[is.na(ppath)] <- 0.5
		if (ppath < alpha/N) {
			ftm <- lapply(1:r, function(x) {
						data.frame(from = vpath[x], to = vpath[x + 1])
					})
		} else {
			ftm <- NULL
		}
		do.call(rbind, lapply(ftm, as.data.frame))
	}
	
	x <- split(ftm, f = seq(nrow(ftm)))
	message("Edge weigthing of ", length(x), " edges ...")
	op <- pbapply::pboptions(type = "timer", style = 2)
	
	if (length(x) > limit) {
		n_cores <- parallel::detectCores()/2
		cl <- parallel::makeCluster(n_cores)
		parallel::clusterExport(cl, c("local", "graph", "alpha", "N"),
		                        envir = environment())
		ftm <- pbapply::pblapply(x, local, cl = cl)
		parallel::stopCluster(cl)
	} else {
		est <- pbapply::pblapply(x, local, cl = NULL)
	}
	ftm <- do.call(rbind, lapply(est, as.data.frame))
	
	# Merging the shortest paths
	if( !is.null(ftm) ) {
		del <- which(duplicated(ftm) == TRUE)
		if(length(del) > 0) ftm <- ftm[-del,]
		Gs <- simplify(graph_from_data_frame(ftm, directed = FALSE))
		if(is.directed(graph)) Gs <- orientEdges(ug = Gs, dg = graph)
	} else {
		Gs <- make_empty_graph(0)
	}
	
	return(Gs)
}

TMFG <- function(cormat, ...)
{
    n <- ncol(cormat)
    cormat <- abs(cormat)
    in_v <- matrix(nrow = nrow(cormat), ncol = 1)
    ou_v <- matrix(nrow = nrow(cormat), ncol = 1)
    tri <- matrix(nrow = ((2 * n) - 4), ncol = 3)
    separators <- matrix(nrow = n - 4, ncol = 3)
    s <- rowSums(cormat*(cormat > mean(matrix(unlist(cormat),nrow = 1)))*1)
    in_v[1:4] <- order(s, decreasing = TRUE)[1:4]
    ou_v <- setdiff(1:nrow(in_v), in_v)
    tri[1,] <- in_v[1:3,]
    tri[2,] <- in_v[2:4,]
    tri[3,] <- in_v[c(1, 2, 4),]
    tri[4,] <- in_v[c(1, 3, 4),]
	
    S <- matrix(nrow = (3 * nrow(cormat) - 6), ncol = 3)
    if (cormat[in_v[1], in_v[2]] > cormat[in_v[2], in_v[1]]) {
		S[1,] <- c(in_v[1], in_v[2], 1)
    } else {
		S[1, ] <- c(in_v[2], in_v[1], 1)
    }
    
    if (cormat[in_v[1], in_v[3]] > cormat[in_v[3], in_v[1]]) {
		S[2,] <- c(in_v[1], in_v[3], 1)
    } else {
		S[2,] <- c(in_v[3], in_v[1], 1)
    }
    
    if (cormat[in_v[1], in_v[4]] > cormat[in_v[4], in_v[1]]) {
		S[3,] <- c(in_v[1], in_v[4], 1)
    } else {
		S[3,] <- c(in_v[4], in_v[1], 1)
    }
    
    if (cormat[in_v[2], in_v[3]] > cormat[in_v[3], in_v[2]]) {
		S[4,] <- c(in_v[2], in_v[3], 1)
    } else {
		S[4,] <- c(in_v[3], in_v[2], 1)
    }
    
    if (cormat[in_v[2], in_v[4]] > cormat[in_v[4], in_v[2]]) {
		S[5,] <- c(in_v[2], in_v[4], 1)
    } else {
		S[5,] <- c(in_v[4], in_v[2], 1)
    }
    
    if (cormat[in_v[3], in_v[4]] > cormat[in_v[4], in_v[3]]) {
            S[6,] <- c(in_v[3], in_v[4], 1)
    } else {
            S[6,] <- c(in_v[4], in_v[3], 1)
    }
    
	 gain <- matrix(-Inf, nrow = n, ncol = (2 * (n - 2)))
    gain[ou_v, 1] <- rowSums(cormat[ou_v, (tri[1,])])
    gain[ou_v, 2] <- rowSums(cormat[ou_v, (tri[2,])])
    gain[ou_v, 3] <- rowSums(cormat[ou_v, (tri[3,])])
    gain[ou_v, 4] <- rowSums(cormat[ou_v, (tri[4,])])
    ntri <- 4
    gij <- matrix(nrow = 1, ncol = ncol(gain))
    v <- matrix(nrow = 1, ncol = ncol(gain))
    ve <- array()
    tr <- 0
    for (e in 5:n) {
		if (length(ou_v) == 1) {
			ve <- ou_v
			v <- 1
			w <- 1
			tr <- which.max(gain[ou_v,])
        } else {
			for (q in 1:ncol(gain)) {
				gij[, q] <- max(gain[ou_v, q])
				v[, q] <- which.max(gain[ou_v, q])
				tr <- which.max(gij)
            }
            ve <- ou_v[v[tr]]
            w <- v[tr]
        }
        ou_v <- ou_v[-w]
        in_v[e] <- ve
        for (u in 1:length(tri[tr, ])) {
            cou <- 6 + ((3 * (e - 5)) + u)
            S[cou, ] <- cbind(ve, tri[tr, u], 1)
        }
        separators[e - 4, ] <- tri[tr, ]
        tri[ntri + 1, ] <- cbind(rbind(tri[tr, c(1, 3)]), ve)
        tri[ntri + 2, ] <- cbind(rbind(tri[tr, c(2, 3)]), ve)
        tri[tr, ] <- cbind(rbind(tri[tr, c(1, 2)]), ve)
        gain[ve, ] <- 0
        gain[ou_v, tr] <- rowSums(cormat[ou_v, tri[tr, ], drop = FALSE])
        gain[ou_v, ntri + 1] <- rowSums(cormat[ou_v, tri[ntri + 1, ],
                                        drop = FALSE])
        gain[ou_v, ntri + 2] <- rowSums(cormat[ou_v, tri[ntri + 2, ],
                                        drop = FALSE])
        ntri <- ntri + 2
    }
   	cliques <- rbind(in_v[1:4], (cbind(separators, in_v[5:ncol(cormat)])))
    L <- S
    L[, 1] <- S[, 2]
    L[, 2] <- S[, 1]
    K <- rbind(S, L)
    x <- as.matrix(Matrix::sparseMatrix(i = K[, 1], j = K[, 2], x = K[, 3]))
    diag(x) <- 1
    for (r in 1:nrow(x)) for (z in 1:ncol(x)) {
        if (x[r, z] == 1) {
            x[r, z] <- cormat[r, z]
        }
    }
    colnames(x) <- colnames(cormat)
    x <- as.data.frame(x)
    row.names(x) <- colnames(x)
    x <- as.matrix(x)
	x <- x - diag(nrow(x))
	gtmf <- graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)
    
	return(list(graph = gtmf, separators = separators, cliques = cliques))
}
