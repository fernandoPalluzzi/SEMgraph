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
#' (-1 if z<-2, +1 if z>2, 0 otherwise) and P-value of the z-test = b/SE(b)
#' from simple linear regression y ~ x (i.e., glm(node ~ group)).
#' @param method Edge weighting method. It can be one of the following:
#' \enumerate{
#' \item "r2z", weight edges of a graph using Fisher's r-to-z transform 
#' (Fisher, 1915) to test the correlation coefficient of pairs 
#' of interacting nodes, if \code{group == NULL}. Otherwise, the correlation
#' difference between group will be tested and edge weights correspond to
#' the sign (-1 if z<-2, +1 if z>2, 0 otherwise) and P-value of the group
#' correlation difference.
#' \item "sem", edge weights are defined by a SEM model that implies 
#' testing the group effect simultaneously on the source node and the sink mode.
#' A new parameter w is defined as the weighted sum of the total effect 
#' of the group on source and sink nodes, adjusted by node degree centrality, 
#' and edge weights correspond to the sign (-1 if z<-2, +1 if z>2, 0 otherwise)
#' and P-value of the z-test = w/SE(w). Not available if \code{group == NULL}.
#' \item "cov", edge weights are defined by a new parameter w combining 
#' the group effect on the source node (mean group difference, adjusted 
#' by source degree centrality), the sink node (mean group difference, 
#' adjusted by sink degree centrality), and the source-sink interaction 
#' (correlation difference). Edge weights correspond to the sign
#' (-1 if z<-2, +1 if z>2, 0 otherwise) and P-value of the z-test = w/SE(w)
#' of the combined difference of the group over source node, sink node, and
#' their connection. Not available if \code{group == NULL}.
#' \item "cfa", edge weights are defined by a CFA1 model that implies 
#' testing the group effect, w on a latent variable (LV) with observed indicators
#' two interacting nodes, fixing loading coefficients and residual
#' variances for model identification. Edge weights correspond to the sign
#' (-1 if z<-2, +1 if z>2, 0 otherwise) and P-value of the z-test = w/SE(w)
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
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
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
	ftm <- igraph::as_data_frame(ig)
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
