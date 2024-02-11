#  SEMgraph library
#  Copyright (C) 2019-2024 Mario Grassi; Fernando Palluzzi; Barbara Tarantino 
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
#' node weighting is actived, and node weights correspond to the attribute:
#' V(graph)$pv (P-value of the z-test = b/SE(b) from simple linear regression
#' y ~ x, i.e., lm(node ~ group)) and V(graph)$sign (-1 if z<-2, +1 if z>2,
#' 0 otherwise).
#' @param method Edge weighting method. It can be one of the following:
#' \enumerate{
#' \item "r2z", weight edges are defined using Fisher's r-to-z transform
#' (Fisher, 1915) to test the correlation coefficient of pairs of interacting
#' nodes, if \code{group=NULL}. Otherwise, the difference between group of
#' the r-to-z trasform will be tested. Edge weights correspond to the attribute:
#' E(graph)$pv (P-value of the z-test) and E(graph)$sign (-1 if z<-2, +1 if z>2,
#' 0 otherwise).
#' \item "sem", edge weights are defined by a SEM model that implies 
#' testing the group effect simultaneously on source and sink nodes.
#' A new parameter w is defined as the weighted sum of the total effect 
#' of the group on source and sink nodes, adjusted by node degree centrality. 
#' Edge weights correspond to the attribute: E(graph)$pv (P-value of the
#' z-test = w/SE(w)) and E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise). 
#' Not available if \code{group=NULL}.
#' \item "cov", edge weights are defined by a new parameter w combining 
#' the group effect on the source node (mean group difference, adjusted 
#' by source degree centrality), the sink node (mean group difference, 
#' adjusted by sink degree centrality), and the source--sink interaction 
#' (correlation difference). Edge weights correspond to the attribute:
#' E(graph)$pv (P-value of the z-test = w/SE(w) of the combined difference
#' of the group over source node, sink node, and their connection) and
#' E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise).
#' Not available if \code{group=NULL}.
#' \item "cfa", edge weights are defined by a CFA1 model that implies 
#' testing the group effect, w on a latent variable (LV) with observed
#' indicators two interacting nodes, fixing loading coefficients and residual
#' variances for model identification. Edge weights correspond to the
#' attribute: E(graph)$pv (P-value of the z-test = w/SE(w) of the group
#' effect on the LV) and E(graph)$sign (-1 if z<-2, +1 if z>2, 0 otherwise).
#' Not available if \code{group=NULL}.
#' }
#' @param limit An integer value corresponding to the number of graph 
#' edges. Beyond this limit, multicore computation is enabled to reduce 
#' the computational burden. By default, \code{limit = 10000}.
#' @param ... Currently ignored.
#'
#' @return A weighted graph, as an igraph object.
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Grassi M, Tarantino B (2023). [Supplementary material of] SEMtree: tree-based structure
#' learning methods with structural equation models. 
#' Bioinformatics, 39 (6), 4829–4830 <https://doi.org/10.1093/bioinformatics/btad377>
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
#' # New node attributes
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

#' @title Transform data methods
#'
#' @description Implements various data trasformation methods with
#' optimal scaling for ordinal or nominal data, and to help relax
#' the assumption of normality (gaussianity) for continuous data.
#' 
#' @param x A matrix or data.frame (n x p). Rows correspond to subjects, and
#' columns to graph nodes.
#' @param method Trasform data method. It can be one of the following:
#' \enumerate{
#' \item "npn" (default), performs nonparanormal(npn) or semiparametric
#' Gaussian copula model (Liu et al, 2009), estimating the Gaussian copula
#' by marginally transforming the variables using smooth ECDF functions.
#' The npn distribution corresponds to the latent underlying multivariate
#' normal distribution, preserving the conditional independence structure
#' of the original variables.
#' \item "spearman", computes a trigonometric trasformation of Spearman
#' rho correlation for estimation of latent Gaussian correlations
#' parameter of a nonparanormal distribution (Harris & Dorton (2013),
#' and generates the data matrix with the exact same sample covariance
#' matrix as the estimated one.
#' \item "kendall", computes a trigonometric trasformation of Kendall
#' tau correlation for estimation of latent Gaussian correlations
#' parameter of a nonparanormal distribution (Harris & Dorton (2013),
#' and generates the data matrix with the exact same sample covariance
#' matrix as the estimated one.
#' \item "polichoric", computes the polychoric correlation matrix and
#' generates the data matrix with the exact same sample covariance matrix
#' as the estimated one. The polychoric correlation (Olsson, 1974) is a
#' measure of association between two ordinal variables. It is based on the
#' assumption that two latent bivariate normally distributed random variables
#' generate couples of ordinal scores. Tetrachoric (two binary variables) and
#' biserial (an ordinal and a numeric variables) correlations are special cases.
#' \item "lineals", performs optimal scaling in order to achieve linearizing
#' transformations for each bivariate regression between pairwise variables for
#' subsequent structural equation models using the resulting correlation
#' matrix computed on the transformed data (de Leeuw, 1988).
#' \item "mca", performs optimal scaling of categorical data by Multiple
#' Correspondence Analysis (MCA, a.k.a homogeneity analysis) maximizing
#' the first eigenvalues of the trasformed correlation matrix. The estimates
#' of the corresponding structural parameters are consistent if the underlying
#' latent space of the observed variables is unidimensional.
#' }
#' @param ... Currently ignored.
#'
#' @details Nonparanormal trasformation is computationally very efficient
#' and only requires one ECDF pass of the data matrix. Polychoric correlation
#' matrix is computed with the \code{lavCor()} function of the \code{lavaan}
#' package. Optimal scaling (lineals and mca) is performed with the
#' \code{lineals()} and \code{corAspect()} functions of the \code{aspect}
#' package (Mair and De Leeuw, 2008). To note, SEM fitting of the generate data
#' (fake data) must be done with a covariance-based method and bootstrap SE,
#' i.e., with \code{SEMrun(..., algo="ricf", n_rep=1000)}.
#' 
#' @return A list of 2 objects is returned:
#' \enumerate{
#' \item "data", the matrix (n x p) of n observations and p transformed
#' variables or the matrix (n x p) of simulate observations based on the
#' selected correlation matrix.  
#' \item "catscores", the category weights for "lineals" or "mca"
#' methods or NULL otherwise.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' 
#' Liu H, Lafferty J, and Wasserman L (2009). The Nonparanormal: Semiparametric Estimation of
#' High Dimensional Undirected Graphs. Journal of Machine Learning Research 10(80): 2295-2328
#' 
#' Harris N, and Drton M (2013). PC Algorithm for Nonparanormal Graphical Models.
#' Journal of Machine Learning Research 14 (69): 3365-3383
#' 
#' Olsson U (1979). Maximum likelihood estimation of the polychoric correlation coefficient.
#' Psychometrika, 44(4), 443-460.
#' 
#' Mair P, and De Leeuw J (2008). Scaling variables by optimizing correlational and
#' non-correlational aspects in R. Journal of Statistical Software, 32(9), 1-23.
#' 
#' de Leeuw J (1988). Multivariate analysis with linearizable regressions. Psychometrika,
#' 53, 437-454.
#'
#' @examples
#' 
#' #... with continuous ALS data
#' graph<- alsData$graph
#' data<- alsData$exprs; dim(data)
#' X<- data[, colnames(data) %in% V(graph)$name]; dim(X)
#'
#' npn.data<- transformData(X, method="npn")
#' sem0.npn<- SEMrun(graph, npn.data$data)
#'
#' mvnS.data<- transformData(X, method="spearman")
#' sem0.mvnS<- SEMrun(graph, mvnS.data$data)
#'
#' mvnK.data<- transformData(X, method="kendall")
#' sem0.mvnK<- SEMrun(graph, mvnK.data$data)
#' 
#' #...with ordinal (K=4 categories) ALS data
#' Xord <- data.frame(X)
#' Xord <- as.data.frame(lapply(Xord, cut, 4, labels = FALSE))
#' colnames(Xord) <- sub("X", "", colnames(Xord))
#' 
#' \dontrun{
#'
#' mvnP.data<- transformData(Xord, method="polychoric")
#' sem0.mvnP<- SEMrun(graph, mvnP.data$data, algo="ricf", n_rep=1000)
#'
#' }
#'
#' lin.data<- transformData(Xord, method="lineals")
#' sem0.lin<- SEMrun(graph, lin.data$data)
#' lin.data$catscores; head(lin.data$data)
#'
#' #...with nominal (K=4 categories) ALS data
#' mca.data<- transformData(Xord, method="mca")
#' sem0.mca<- SEMrun(graph, mca.data$data)
#' mca.data$catscores
#' 
#' # plot colored graphs
#' #par(mfrow=c(3,2), mar=rep(1,4))
#' #gplot(sem0.npn$graph, l="fdp", main="ALS npm")
#' #gplot(sem0.mvnS$graph, l="fdp", main="ALS mvnS")
#' #gplot(sem0.mvnK$graph, l="fdp", main="ALS mvnK")
#' #gplot(sem0.mvnP$graph, l="fdp", main="ALS mvnP")
#' #gplot(sem0.lin$graph, l="fdp", main="ALS lin")
#' #gplot(sem0.mca$graph, l="fdp", main="ALS mca")
#'
transformData <- function (x, method = "npn", ...)
{
	n <- nrow(x)
	p <- ncol(x)
	x.col <- colnames(x)
	x.row <- rownames(x)
	catscores <- NULL
	
	if (method == "npn") {
		cat("Conducting the nonparanormal transformation via shrunkun ECDF...")
		z <- qnorm(apply(x, 2, rank)/(n + 1))
		z <- z/sd(z[, 1])	
	}
	if (method == "spearman") { 
		cat("Simulating gaussian data via Spearman correlations...")
		x <- 2 * sin(pi/6 * cor(x, method = "spearman"))
		z <- generateData(Sest = x, n = n, p = p)
    }
	if (method == "kendall") { 
		cat("Simulating gaussian data via Kendall correlations...")
		x <- sin(pi/2 * cor(x, method = "kendall"))
		z <- generateData(Sest = x, n = n, p = p)
    }
	if (method == "polychoric") { 
		cat("Simulating gaussian data via polychoric correlations...")
		#x <- suppressWarnings(lavCor(x, ordered = names(x))[1:p,1:p])
		x <- lavCor(x, ordered=names(x), cor.smooth=TRUE)[1:p,1:p]
		z <- generateData(Sest = x, n = n, p = p)
    }
	if (method == "lineals") {
		cat("Conducting the optimal (ordinal) linearizing transformation...")
		z <- aspect::lineals(x, level = "ordinal")
		catscores<- z$catscores
		z <- sqrt(n-1)*z$scoremat
    }
	if (method == "mca") {
		cat("Conducting the first solution of Multiple Correspondence Analysis...")
		z <- aspect::corAspect(x, aspect = "aspectEigen")
		catscores<- z$catscores
		z <- sqrt((n-1)/n)*z$scoremat
	}

	cat("done.\n")
	colnames(z) <- x.col
    rownames(z) <- x.row
	v0 <- which(apply(z, 2, var) == 0)
	if (length(v0) > 0) z <- cbind(z[,-v0], x[,v0])

	return(list(data = z, catscores = catscores))
}

generateData <- function(Sest, n, p, ...)
{
	if (!corpcor::is.positive.definite(Sest)){
	 Sest <- corpcor::cor.shrink(Sest, verbose = FALSE)[1:p,1:p]
	}
	e <- eigen(Sest)
	sqrt.true.cov.mat <- e$vectors%*%sqrt(diag(e$values))

	fake.data <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sest)
	samp.cov.mat <- cov(fake.data)
	e <- eigen(samp.cov.mat)
	sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))

	fake.data <- t(sqrt.true.cov.mat%*%solve(sqrt.samp.cov.mat,t(fake.data)))
	fake.data <- as.data.frame(fake.data)

	return(fake.data)
}

#' @title SEM-based out-of-sample predictions
#'
#' @description Given the values of (observed) x-variables in a structural equation
#' model, this function may be used to predict the values of (observed) y-variables.
#' Response variables (y) represent sink nodes, and predictor variables (x)
#' might consist of either (i) just source nodes or (ii) source and mediators from
#' the fitted graph structure.
#'
#' @param object An object, as that created by the function \code{SEMrun()} with the
#' argument \code{fit} set to \code{fit = 0} or \code{fit = 1}.
#' @param newdata An optional matrix with rows corresponding to subjects, and
#' columns to graph nodes (variables). If \code{object$fit} is a model with the
#' group variable (\code{fit = 1}), the first column of newdata must be the
#' new group binary vector (0=control, 1=case). As a default \code{newdata = NULL}, 
#' meaning that the K-fold cross validation is applied on the \code{object$data}.
#' Conversely, if the argument \code{newdata} is specified, this matrix will be
#' used for testing (out-of-sample predictions) and \code{object$data} will be
#' used for training.
#' @param K_fold The number of subsets (folds) into which the data will be 
#' partitioned for performing K-fold cross-validation. The model is refit K times, 
#' each time leaving out one of the K folds (default, K_fold=5). If the argument 
#' \code{newdata} is specified, the K-fold cross validation will not be done.
#' @param source A logical value. If FALSE (default), the predictor variables (x) 
#' include source and mediators. If TRUE, x includes only the source nodes. 
#' @param verbose A logical value. If FALSE (default), the processed graph 
#' will not be plotted to screen.
#' @param ... Currently ignored.
#'
#' @details The function uses a SEM-based predictive approach (Rooij et al., 2022)
#' to produce predictions while accounting for the given graph structure. Predictions
#' (for y given x) are based on the (joint y and x) model-implied variance-covariance
#' (Sigma) matrix and mean vector (Mu) of the fitted SEM, and the standard expression
#' for the conditional mean of a multivariate normal distribution. Thus, the structure
#' described in the SEM is taken into consideration, which differs from ordinary
#' least squares (OLS) regression. Note that if the model is saturated (and hence
#' df = 0), or when \code{source = TRUE}, i.e., the set of predictors will include only
#' the source nodes, the SEM-based predictions are identical or similar to OLS
#' predictions. 
#'
#' @return A list of 3 objects:
#' \enumerate{
#' \item "yobs", the matrix of observed continuous values of sink nodes based on
#' out-of-bag samples. 
#' \item "yhat", the matrix of continuous predicted values of sink nodes ased on
#' out-of-bag samples.
#' \item "PE", vector of the prediction error equal to the Root Mean Squared Error
#' (RMSE) for each out-of-bag sink prediction. The first value of PE is the total
#' RMSE, where we sum over all sink nodes.
#' }
#'
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references 
#' 
#' de Rooij M, Karch JD, Fokkema M, Bakk Z, Pratiwi BC, and Kelderman H
#' (2023). SEM-Based Out-of-Sample Predictions, Structural Equation Modeling:
#' A Multidisciplinary Journal, 30:1, 132-148
#' <https://doi.org/10.1080/10705511.2022.2061494>
#' 
#' @examples
#'
#' # load ALS data
#' ig<- alsData$graph
#' X<- alsData$exprs
#' X<- transformData(X)$data
#' group<- alsData$group
#' 
#' #...with train-test (0.8-0.2) samples
#' set.seed(1)
#' train<- sample(1:nrow(X), 0.8*nrow(X))
#'
#' # SEM fitting
#' #sem0<- SEMrun(ig, X[train,], algo="lavaan", SE="none")
#' #sem0<- SEMrun(ig, X[train,], algo="ricf", n_rep=0)
#' sem0<- SEMrun(ig, X[train,], algo="cggm")
#' 
#' # predictors, source+mediator variables
#' res1<- predictSink(sem0, newdata=X[-train,]) 
#' print(res1$PE)
#' 
#' # predictors, source variables
#' res2<- predictSink(sem0, newdata=X[-train,], source=TRUE) 
#' print(res2$PE)
#' 
#' #...with 5-fold cross-validation samples
#' set.seed(2)
#'
#' # SEM fitting
#' #sem0<- SEMrun(ig, X, algo="lavaan", SE="none")
#' #sem0<- SEMrun(ig, X, algo="ricf", n_rep=0)
#' sem0<- SEMrun(ig, X, algo="cggm")
#' 
#' # predictors, source+mediator variables	
#' res3<- predictSink(sem0, K_fold = 5, verbose=TRUE)
#' print(res3$PE)
#' 
#' # predictors, source variables
#' res4<- predictSink(sem0, K_fold = 5, source=TRUE, verbose=TRUE) 
#' print(res4$PE)
#'
#' \dontrun{
#'
#' #...with 10-fold cross-validation samples and 10-iterations
#' 
#' # SEM fitting
#' #sem1<- SEMrun(ig, X, group, algo="lavaan", SE="none")
#' #sem1<- SEMrun(ig, X, group, algo="ricf", n_rep=0)
#' sem1<- SEMrun(ig, X, group, algo="cggm")
#' 
#' # predictors, source+mediator+group variables
#' res<- NULL
#' for (r in 1:10) {
#' 	set.seed(r)
#' 	cat("rep = ", r, "\n")
#' 	resr<- predictSink(sem1, K_fold = 10)
#' 	res<- rbind(res, resr$PE)
#' }
#' res
#' apply(res, 2, mean)
#'
#' }
#'
predictSink <- function(object, newdata = NULL, K_fold = 5, source = FALSE,
						verbose = FALSE, ...)
{
	# set graph, predictors and outcomes
	stopifnot(inherits(object$fit, c("lavaan", "RICF", "GGM")))
	#stop("ERROR: in SEMrun(..., fit = .) cannot be fit=2 (for now).")
	graph<- object$graph
	data<- object$data
	if (!is.na(data[1,1])) graph<- map_group(graph)
	if (!is.null(newdata)) {
		idx<- colnames(data)[colnames(data) %in% colnames(newdata)]
		train <- data[,idx]
		test <- newdata[,idx]
		data<- rbind(train, test)
	}
	idv<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids = V(graph)$name %in% idv)
	din<- igraph::degree(graph, mode = "in")
	dout<- igraph::degree(graph, mode = "out")
	V(graph)$color[din == 0]<- "cyan"
	V(graph)$color[dout == 0]<- "orange"
	yn<- V(graph)$name[dout == 0] #sink
	xn<- V(graph)$name[dout != 0] #source+mediator
	if (source == TRUE) {
		xn<- V(graph)$name[din == 0] #source
		graph<- map_source(xn, yn)
	}
	if (verbose) gplot(graph)

	# K-fold cross-validation train indices
	if (is.null(newdata)) {
		Y <- rowMeans(data[,yn])
		idx <- createFolds(y=Y, k=K_fold, list=TRUE, returnTrain=TRUE)
	}else{
		K_fold <- 1
		idx<- list(1:nrow(train))
	}

	# SEM fitting on train data and predition on test data
	yobs<- NULL
	yhat<- NULL

	for (k in 1:K_fold) { #k=1
	  if (K_fold != 1) {
		message("Fold: ", k) 
		fit<- quiet(SEMrun(graph, data[idx[[k]],], algo="ricf", n_rep=0))
	  }else{
		fit<- object
	  }
	  if (inherits(fit$fit, "lavaan")) {
		xnames<- paste0("z", xn)
		ynames<- paste0("z", yn)
		Sxx<- fitted(fit$fit)$cov[xnames, xnames]
		Sxy<- fitted(fit$fit)$cov[xnames, ynames]
		#mx<- fitted(fit$fit)$mean[xnames]
		#my<- fitted(fit$fit)$mean[ynames]
	  }else{
		Sxx<- fit$fit$Sigma[xn, xn]
		Sxy<- fit$fit$Sigma[xn, yn]
	  }
	  mx<- rep(0, length(xn))
	  my<- rep(0, length(yn))
	  xtest<- as.matrix(data[-idx[[k]], xn])
	  xtest<- scale(xtest, center = mx, scale = TRUE)
	  n<- nrow(xtest)
	  py<- length(yn)
	  My<- matrix(my, n, py, byrow = TRUE)
	  if (corpcor::is.positive.definite(Sxx)) {
		yhatk<- My + xtest %*% solve(Sxx) %*% Sxy
	  }else{
		yhatk<- My + xtest %*% Sxy
	  }
	  yobs<- rbind(yobs, scale(data[-idx[[k]], yn]))
	  yhat<- rbind(yhat, yhatk)
	}

	PE<- sqrt(colMeans((yobs - yhat)^2))
	PE<- c(RMSEp=sqrt(mean(PE^2)), PE)
	if (verbose) print(PE)

	return(list(yobs=yobs, yhat=yhat, PE=PE))
}

createFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE, ...) 
{
    #createFolds() function from "caret" package (author: Max Kuhn)
	#All rights reserved. See the file COPYING for license terms.

	if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
        y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) { #i=1
            min_reps <- numInClass[i]%/%k
            if (min_reps > 0) {
                spares <- numInClass[i]%%k
                seqVector <- rep(1:k, min_reps)
                if (spares > 0) 
                  seqVector <- c(seqVector, sample(1:k, spares))
                foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
            }
            else {
                foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                  size = numInClass[i])
            }
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
            sep = "")
        if (returnTrain) 
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
	
    return(out)
}

map_source <- function(xn, yn, verbose=FALSE, ...)
{
	gout <- make_empty_graph(length(c(xn,yn)))
	V(gout)$name <- c(xn,yn)
	E <- NULL
	 for(k in 1:length(xn)){
	  for(j in 1:length(yn)){ #i=1
		E <- c(E, xn[k], yn[j])
	  }
	 }
	gout <- gout + igraph::edges(E)  
	if (verbose) gplot(gout)
	
	return(gout)	
}

map_group <- function(graph, verbose=FALSE, ...)
{
	gout <- graph + igraph::vertices("group")
	nodes<- V(graph)$name
	E <- NULL
	 for(v in 1:length(nodes)){
	  	E <- c(E, "group", nodes[v])
	  }
	gout <- gout + igraph::edges(E)
	#V(gout)$color[V(gout)$name == "outcome"] <- "green"
	if (verbose) gplot(gout)
	
	return(gout)	
}

map_outcome <- function(graph, verbose=FALSE, ...)
{
	gout <- graph + igraph::vertices("outcome")
	dout<- igraph::degree(graph, mode = "out")
	leaf<- V(graph)$name[dout == 0]
	E <- NULL
	 for(v in 1:length(leaf)){
	  	E <- c(E, leaf[v], "outcome")
	  }
	gout <- gout + igraph::edges(E)
	V(gout)$color[V(gout)$name == "outcome"] <- "green"
	if (verbose) gplot(gout)
	
	return(gout)	
}

map_LV <- function(graph, data, cg=NULL, verbose=FALSE, ...)
{
	VH <- colnames(data)[grepl("LV",colnames(data))]
	gH <- graph + igraph::vertices(VH)
	E <- NULL
	 for(v in 1:length(VH)){
	  if (!is.null(cg)) graph<- cg[[v]]
	  for(i in 1:vcount(graph)){ #i=1
		E <- c(E, VH[v], V(graph)$name[i])
	  }
	 }
	gH <- gH + igraph::edges(E)
	if (verbose) gplot(gH, l="fdp")
	
	return(gH)	
}
