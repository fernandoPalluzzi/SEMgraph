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

#' @title Compute the Average Causal Effect (ACE) for a given source-sink pair
#'
#' @description Compute total effects as ACEs of variables X
#' on variables Y in a directed acyclic graph (DAG). The ACE will be estimated
#' as the path coefficient of X (i.e., theta) in the linear equation
#' Y ~ X + Z. The set Z is defined as the adjustment (or conditioning) set of
#' Y over X, applying various adjustement sets. Standard errors (SE),
#' for each ACE, are computed following the \code{lm} standard procedure
#' or a bootstrap-based procedure (see \code{\link[boot]{boot}} for details).
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If \code{group = NULL} (default), group influence
#' will not be considered.
#' @param type character Conditioning set Z. If "parents" (default) the
#' Pearl's back-door set (Pearl, 1998), "minimal" the dagitty minimal set
#' (Perkovic et al, 2018), or "optimal" the O-set with the smallest
#' asymptotic variance (Witte et al, 2020) are computed.
#' @param effect character X to Y effect. If "all" (default) all effects from
#' X to Y, "source2sink" only effects from source X to sink Y, or "direct"
#' only direct effects from X to Y are computed.
#' @param method Multiple testing correction method. One of the values
#' available in \code{\link[stats]{p.adjust}}.
#' By default, \code{method = "BH"} (i.e., FDR multiple test correction).
#' @param alpha Significance level for ACE selection (by default,
#' \code{alpha = 0.05}).
#' @param boot The number of bootstrap samplings enabling bootstrap
#' computation of ACE standard errors. If \code{NULL} (default), bootstrap
#' is disabled.
#' @param ... Currently ignored.
#'
#' @return A data.frame of ACE estimates between network sources and sinks.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Pearl J (1998). Graphs, Causality, and Structural Equation Models.
#' Sociological Methods & Research, 27(2):226-284.
#' <https://doi.org/10.1177/0049124198027002004>
#'
#' Perkovic E, Textor J, Kalisch M, Maathuis MH (2018). Complete graphical
#' characterization and construction of adjustment sets in Markov equivalence
#' classes of ancestral graphs. Journal of Machine Learning Research, 18:1-62.
#' <http://jmlr.org/papers/v18/16-319.html>
#'
#' Witte J, Henckel L, Maathuis MH, Didelez V (2020). On efficient
#' adjustment in causal graphs. Journal of Machine Learning Research, 21:1-45.
#' <http://jmlr.org/papers/v21/20-175.htm>
#'
#' @examples
#'
#' # ACE without group, O-set, all effects:
#' ace1 <- SEMace(graph = sachs$graph, data = log(sachs$pkc),
#'                group = NULL, type = "optimal", effect = "all",
#'                method = "BH", alpha = 0.05, boot = NULL)
#' print(ace1)
#'
#' # ACE with group perturbation, Pa-set, direct effects:
#' ace2 <- SEMace(graph = sachs$graph, data = log(sachs$pkc),
#'                group = sachs$group, type = "parents", effect = "direct",
#'                method = "none", alpha = 0.05, boot = NULL)
#' print(ace2)
#'
SEMace<- function(graph, data, group=NULL, type="parents", effect="all", method="BH", alpha=0.05, boot=NULL, ...)
{
	# Set igraph and dagitty graph objects
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	
	if( !is_dag(ig) & type != "parents" ){
	 cat("\nWARNING: input graph is not acyclic!\n")
     cat(" Applying graph -> DAG conversion.\n")
	 dag<- graph2dag(ig, data) #del cycles & all <->
	 if ( !is_dag(dag) ) return (theta=NULL)
	}else{ dag<- ig }
	dagy<- graph2dagitty(dag, verbose=FALSE)
	#plot(dagitty::graphLayout(dagy))
	
	# Set distance matrix, distance graph from source to target nodes
	D<- distances(dag, mode="out", weights=NA)
	D<- ifelse(D == Inf, 0, D) #sum(D>0)
	if (effect == "all") {
	 gD<- igraph::simplify(graph_from_adjacency_matrix(D, mode="directed", weighted=TRUE))
	}
	if (effect == "source2sink") {
	 din<- igraph::degree(dag, mode= "in")
	 Vx<- V(dag)$name[din == 0]
	 dout<- igraph::degree(dag, mode= "out")
	 Vy<- V(dag)$name[dout == 0]
	 Dxy<- D[c(Vx,Vy),c(Vx,Vy)]
	 gD<- igraph::simplify(graph_from_adjacency_matrix(Dxy, mode="directed", weighted=TRUE))
	}
	if (effect == "direct") {
	 Dd<- ifelse(D == 1, 1, 0)
	 gD<- igraph::simplify(graph_from_adjacency_matrix(Dd, mode="directed", weighted=TRUE))
	}
	cat("\n","Frequency distribution of path length from X to Y :")
	print(table(E(gD)$weight)); cat("\n")

	# Compute total effect (theta=ACE)
	theta<- NULL
	res<- NULL
	ftm<- igraph::as_data_frame(gD)
	for (i in 1:nrow(ftm)) {
	 cat("\r","ACE=", i, "of", nrow(ftm)) 
     flush.console()
	 x<- ftm[i,1]# x
	 y<- ftm[i,2]# y
	 
	 # Adjustement Z SET:
	  if (type == "parents") { #backdoor (Pearl, 1988)
	  z<- setdiff(V(dag)$name[SEMgraph::parents(dag, x)], x)
	 }
	 if (type == "minimal") { #Perkovic et al (2018)
	  SET<- dagitty::adjustmentSets(dagy, x, y, type="minimal", effect="total")
	  z<- unlist(SET[[1]])
	 }
	 if (type == "optimal") { #Witte et al (2020)
	  #paths<- all_simple_paths(dag, from=x, to=y, mode="out")
	  #cn<- setdiff(unique(names(unlist(paths))),x)
	  paths<- dagitty::paths(dagy, from=x, to=y, directed=TRUE)$paths
	  cn<- setdiff(unique(unlist(strsplit(gsub("->","",paths),"  "))),x)
	  pa_cn<- V(dag)$name[SEMgraph::parents(dag, cn)]
	  forb<- c(V(dag)$name[SEMgraph::descendants(dag, cn)],x)
	  z<- setdiff(pa_cn, forb)
	 }
 
	 # LM fitting:
	 Z<- scale(data[,c(y,x,z)])
	 #if( nrow(Z) < ncol(Z) ) boot <- 1000
	 if( is.null(group) ) {
	  if(is.null(boot)) {
	   try(est<- lmest(x, y, Z))
	  }else{
	   try(est<- boot.lmest(x, y, Z, R=boot))
	  }
	 }else{
	  try(est<- lmest2(x, y, Z, group, boot))
	 }
	 theta<- rbind(theta,cbind(pathL=ftm[i,3],est))
 	}
    cat("\n")
	theta<- subset(theta, p.adjust(theta$pvalue, method=method) < alpha)
	class(theta)<- c("lavaan.data.frame" ,"data.frame")
	return( theta )
}

boot.lmest<- function(x, y, Z, R,...)
{
	# LM fitting y ~ x + Z :
	est<- function(Z, i) { 
	 stats::lm.fit(as.matrix(Z[i,-1]),Z[i,1])$coefficients[1]
	}
	#ncpus<- parallel::detectCores(logical = FALSE)
	#xboot<- boot::boot(Z, est, R=1000, parallel="snow", ncpus=ncpus)
	xboot<- boot::boot(Z, est, R=R)
	t0<- xboot$t0
	se<- sd(xboot$t)
	z<- t0/se
	res<- data.frame(sink=y, op="<-", source=x, est=t0, se=se, z=z,
	 pvalue= 2*(1-pnorm(abs(z))),
	 ci.lower= (t0-1.96*se),
	 ci.upper= (t0+1.96*se))
	return( res )
}

lmest<- function(x, y, Z, ...)
{
	 # LM fitting y ~ x + Z :
	 fit<- stats::lm.fit(as.matrix(Z[,-1]),Z[,1])
	 est<- as.numeric(fit$coefficients)[1]
	 sigma<- sum(fit$residuals^2)/fit$df.residual #sqrt(sigma)
	 #if( is.na(sigma) ) return(NULL)
	 X<- as.matrix(Z[,-1])
	 r<- fit$rank
	 p<- ncol(X)
	 if( r == p ) {
	  se<- sqrt(sigma*diag(solve(t(X)%*%X)))[1]
	 }else{
	  C<- t(X)%*%X # cross-product matrix
	  E<- eigen(C) # Eigenvalues-eigenvectors of C
	  W<- E$vectors[1:p,1:r]%*%diag(1/E$values[1:r])%*%t(E$vectors[1:p,1:r])
	  se<- sqrt(sigma*diag(W))[2]
	 }
	 z<- est/se
	 res<- data.frame(sink= y, op= "<-", source= x, est= est, se= se, z= z,
	  pvalue= 2*(1-pnorm(abs(z))),
	  ci.lower= (est-1.96*se),
	  ci.upper= (est+1.96*se))
	 return( res )
}

lmest2<- function(x, y, Z, group, boot,...)
{
	# LM fitting y ~ x + Z  w/n groups
	Z1<- as.matrix(Z[group==1,])
	Z0<- as.matrix(Z[group==0,])
	#if(nrow(Z1) < ncol(Z1) | nrow(Z0) < ncol(Z0)) boot <- 1000
	if (is.null(boot)) {
	 est1<- lmest(x, y, Z1)
	 est0<- lmest(x, y, Z0)
	}else{
	 est1<- boot.lmest(x, y, Z1, R=boot)
	 est0<- boot.lmest(x, y, Z0, R=boot)
	}
	d_est<- est1$est - est0$est
	d_se<- sqrt(est1$se^2 + est0$se^2)
	d_z<- d_est/d_se
	pvalue<- 2*(1-pnorm(abs(d_z)))
	d_lower<- d_est - 1.96*d_se
	d_upper<- d_est + 1.96*d_se
	res<- cbind(est1[,1:3], d_est, d_se, d_z, pvalue, d_lower, d_upper)
	return( res )
} 

#' @title Search for directed or shortest paths between pairs of source-sink nodes
#'
#' @description Find and fit all directed or shortest paths between two
#' source-sink nodes of a graph.
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If \code{NULL} (default), group influence will
#' not be considered.
#' @param from Starting node name (i.e., source node).
#' @param to Ending node name (i.e., sink node).
#' @param path If \code{path = "directed"}, all directed paths between
#' the two nodes will be included in the fitted model.
#' If \code{path = "shortest"}, only shortest paths will be returned.
#' @param verbose Show the directed (or shortest) path between the
#' given source-sink pair inside the input graph.
#' @param ... Currently ignored.
#'
#' @return A list of four objects: a fitted model object of class
#' \code{\link[lavaan]{lavaan}} ("fit"), aggregated and node-specific
#' group effect estimates and P-values ("gest"), the extracted subnetwork
#' as an igraph object ("graph"), and the input graph with a color
#' attribute mapping the chosen path ("map").
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#'
#' # Directed path fitting
#' path <- SEMpath(graph = sachs$graph, data = log(sachs$pkc),
#'                 group = sachs$group,
#'                 from = "PIP3",
#'                 to = "Erk",
#'                 path = "directed")
#'
#' # Summaries
#' summary(path$fit)
#' print(path$gest)
#'
#' # Graphs
#' gplot(path$map, main="path from PiP2 to Erk")
#' plot(path$map, layout=layout.circle, main="path from PiP2 to Erk")
#'
SEMpath <- function(graph, data, group, from, to, path, verbose = FALSE, ...)
{
	# Set igraph object
	nodes <- colnames(data)
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	if (distances(ig, from, to, mode = "out", weights = NA) == Inf) {
		message("ERROR: infinite distance from", from, "to", to, ".")
		return(NULL)
	}

	if (path == "shortest") {
		# Set shortest path nodes
		paths<- all_shortest_paths(ig, from, to, mode = "out", weights = NA)
		nodes<- unique(names(unlist(paths$res)))
	} else if (path == "directed") {
		# Set directed path nodes
		#paths <- all_simple_paths(ig, from, to, mode = "out")
		#nodes <- unique(names(unlist(paths)))
		dagi <- graph2dagitty(graph2dag(ig, data))
		paths <- dagitty::paths(dagi, from, to, directed = TRUE)$paths
		nodes <- unique(unlist(strsplit(gsub("->", "", paths), "  ")))
	}

	# Plot selected paths
	ig1 <- induced_subgraph(ig, vids = which(V(ig)$name %in% nodes))
	V(ig)$color <- "white"
	V(ig)$color[V(ig)$name %in% V(ig1)$name] <- "gold"
	E(ig)$color <- "gray40"
	E(ig)$color[attr(E(ig), "vnames") %in% attr(E(ig1),
	            "vnames")] <- "darkorange3"
	E(ig)$width <- 1
	E(ig)$width[attr(E(ig), "vnames") %in% attr(E(ig1), "vnames")] <- 2.5
	if (verbose) {
		gplot(ig)
	}

	# SEM fitting
	cat("Path:", from, "->", to, "size-", c(vcount(ig1), ecount(ig1)),
	    "--\n\n")
	if (is.null(group)) {
		sem <- SEMfit(graph = ig1, data = data, group = NULL)
	} else {
		sem <- SEMfit(graph = ig1, data = data, group = group)
	}

	return(list(fit = sem$fit, gest = sem$gest, graph = ig1, map = ig))
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
#' By default, \code{method = "BH"} (i.e., FDR multiple test correction).
#' @param alpha Significance level for ACE selection (by default,
#' \code{alpha = 0.05}).
#' @param ... Currently ignored.
#'
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
#' print(paths$dfp)
#' head(parameterEstimates(paths$fit[[1]]))
#' gplot(paths$paths[[1]])
#'
#' }
#'
pathFinder <- function(graph, data, group = NULL, ace = NULL, path = "directed",
						method = "BH", alpha = 0.05, ...)
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
	if (N == 0) {
	 return(message("STOP: found 0 significant ACEs !"))
	}
	dag <- graph2dag(graph, data)

	for (i in 1:N) { #i=1
	 cat("\r","ACE=", i, "of", N) 
	 flush.console()
	 from <- sources[i]
	 to <- sinks[i]
	 if (distances(dag, from, to, mode = "out", weights = NA) == Inf) next
	 if (distances(dag, from, to, mode = "out", weights = NA) < 2) next
	 fit <- quiet(SEMpath(dag, data, group, from, to, path, verbose = FALSE))
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
					   n.nodes = vcount(fit$graph),
					   n.edges = ecount(fit$graph),
					   dev_df = round(dev_df, 3),
					   srmr = round(srmr, 3),
					   V.pv.act = round(pv1, 6),
					   V.pv.inh = round(pv2, 6))

	 res <- rbind(res, dfp)
	 paths <- c(paths, list(fit$graph))
	 lav <- c(lav, list(fit$fit))
	}
	
	if (is.null(res)){
	 return(message("\nFound 0 significant ACEs with > 2 nodes !"))
	}else{
	 cat("\nFound", nrow(res), "significant ACEs with > 2 nodes\n\n")
	 rownames(res) <- NULL
	 names(paths) <- res$pathId
	 names(lav) <- res$pathId
	}
	
	return(list(paths = paths, fit = lav, dfp = res))
}
