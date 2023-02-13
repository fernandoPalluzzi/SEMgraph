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
#' for graph clustering.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
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
		#LX <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "L"]

		# Latent Variables(LV) model
		K <- as.numeric(names(table(membership)))
		LV <- NULL
		for(k in 1:length(K)) {
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[, which(colnames(dataY) %in% Xk)])
			fa1 <- factor.analysis(Y = Y, r = 1, method = "ml")$Z
			LV <- cbind(LV, fa1)
		}
		colnames(LV) <- paste0("LV", K)
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
		#LY <- V(gLM)$name[substr(V(gLM)$name, 1, 1) == "C"]

		# Composite Variables(CV) model
		K <- as.numeric(names(table(membership)))
		CV <- NULL
		for(k in 1:length(K)) {
			Xk <- subset(names(membership), membership == K[k])
			Y <- as.matrix(dataY[,which(colnames(dataY) %in% Xk)])
			pc1 <- factor.analysis(Y = Y, r = 1, method = "pc")$Z
			CV <- cbind(CV, pc1)
	}
	colnames(CV) <- paste0("CV", K)
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
			spc1 <- factor.analysis(Y = as.matrix(Uk), r = 1,
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

#' @title Factor analysis for high dimensional data
#'
#' @description Wrapper for Factor Analysis with potentially high dimensional variables 
#' implement in the "cate" R package (Author: Jingshu Wang [aut], Qingyuan Zhao [aut, cre]
#' Maintainer: Qingyuan Zhao <qz280@cam.ac.uk>) that is optimized for the high dimensional
#' problem where the number of samples n is less than the number of variables p.
#'
#' @param Y data matrix, a n*p matrix
#' @param r number of factors (default, r =1)
#' @param method algorithm to be used, "pc" (default) or "ml"
#'
#' @details The two methods extracted from "cate" are quasi-maximum likelihood (ml), and
#' principal component analysis (pc). The ml is iteratively solved the EM algorithm
#' using the PCA solution as the initial value. See Bai and Li (2012) for more details.
#'
#' @return a list of objects
#' \describe{
#' \item{Gamma}{estimated factor loadings}
#' \item{Z}{estimated latent factors}
#' \item{Sigma}{estimated noise variance matrix}
#' }
#'
#' @references
#'
#' Jushan Bai and Kunpeng Li (2012). Statistical Analysis of Factor Models of High
#' Dimension. The Annals of Statistics, 40 (1), 436-465
#' <https://doi.org/10.1214/11-AOS966> 
#' 
#' Jingshu Wang and Qingyuan Zhao (2020). cate: High Dimensional Factor Analysis
#' and Confounder Adjusted Testing and Estimation. R package version 1.1.1.
#' <https://CRAN.R-project.org/package=cate>
#'
#' @examples
#' 
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' ## pc
#' pc<- factor.analysis(Y = als.npn, r = 2, method = "pc")
#' head(pc$Gamma)
#' head(pc$Z)
#' head(pc$Sigma)
#' 
#' ## ml
#' ml <- factor.analysis(Y = als.npn, r = 2, method = "ml")
#' head(ml$Gamma)
#' head(ml$Z)
#' head(ml$Sigma)
#'
#' @export
#'
factor.analysis <- function(Y, r = 1, method = "pc") {

	if (r == 0) {
		return(list(Gamma = NULL,
					Z = NULL,
					Sigma = apply(Y, 2, function(v) mean(v^2))))
	}
	if (method == "pc") {
		fa.pc(Y, r)
	} else if (method == "ml") {
		fa.em(Y, r)
	}
}

fa.pc <- function(Y, r) {

	svd.Y <- svd(Y)

	Gamma <- svd.Y$v[, 1:r] %*% diag(svd.Y$d[1:r], r, r) / sqrt(nrow(Y))
	Z <- sqrt(nrow(Y)) * svd.Y$u[, 1:r]

	Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) mean(x^2))

	return(list(Gamma = Gamma,
				Z = Z,
				Sigma = Sigma))
}

fa.em <- function(Y, r, tol = 1e-6, maxiter = 1000) {

    ## A Matlab version of this EM algorithm was in
    ## http://www.mathworks.com/matlabcentral/fileexchange/28906-factor-analysis/content/fa.m

    ## The EM algorithm:
    ##
    ## Y = Z Gamma' + E Sigma^{1/2}                        (n * p)
    ## mle to estimate Gamma (p * r) and Sigma             (p * p)

    ## E step:
    ## EZ = Y (Gamma Gamma' + Sigma)^{-1} Gamma            (n * r)
    ## VarZ = I - Gamma' (Gamma Gamma' + Sigma)^{-1} Gamma (r * r)
    ## EZ'Z = n VarZ + EZ' * EZ                            (r * r)

    ## M step:
    ## update Gamma: Gamma = (Y' EZ)(EZ'Z)^{-1}
    ## update Sigma:
    ## Sigma = 1/n diag(Y'Y - Gamma EZ' Y - Y' EZ' Gamma' + Gamma EZ'Z Gamma')

    ## The log-likelihood (simplified)
    ## llh <- -log det(Gamma Gamma' + Sigma) - tr((Gamma Gamma' + Sigma)^{-1} S)
    ##

    ## For details, see http://cs229.stanford.edu/notes/cs229-notes9.pdf

    p <- ncol(Y)
    n <- nrow(Y)

    ## initialize parameters
	init <- fa.pc(Y, r)
    Gamma <- init$Gamma
	#Gamma <- matrix(runif(p * r), nrow = p)
    #invSigma <- 1/colMeans(Y^2)
    invSigma <- 1/init$Sigma
	#invSigma <- 1/colMeans((Y - sqrt(n) * start.svd$u %*% t(Gamma))^2)
    llh <- -Inf # log-likelihood

    ## precompute quantitites
    I <- diag(rep(1, r))
    ## diagonal of the sample Covariance S
    #sample.var <- apply(Y, 2, var)
	sample.var <- colMeans(Y^2)

    ## compute quantities needed
    tilde.Gamma <- sqrt(invSigma) * Gamma
    M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, symmetric = TRUE)
    YSG <- Y %*% (invSigma * Gamma)

	logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
    B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
    logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
    llh <- -logdetY - logtrY

	converged <- FALSE
    for (iter in 1:maxiter) {

        ## E step:
        ## Using Woodbury matrix identity:
        ## tilde.Gamma = Sigma^{-1/2} Gamma
        ## VarZ = (I + tilde.Gamma' tilde.Gamma)^{-1}
        ## EZ = Y Sigma^{-1} Gamma VarZ
        varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
        EZ <- YSG %*% varZ
        EZZ <- n * varZ + t(EZ) %*% EZ

        ## M step:
        eigenEZZ <- eigen(EZZ, symmetric = TRUE)
        YEZ <- t(Y) %*% EZ
        ## EZ'Z = G'G
        G <- sqrt(eigenEZZ$values) * t(eigenEZZ$vectors)
        ## updating invSigma
        invSigma <- 1/(sample.var - 2/n * rowSums(YEZ * Gamma) +
                           1/n * rowSums((Gamma %*% t(G))^2))
        ## updating Gamma
        Gamma <- YEZ %*% eigenEZZ$vectors %*%
            (1/eigenEZZ$values * t(eigenEZZ$vectors))

        ## compute quantities needed
        tilde.Gamma <- sqrt(invSigma) * Gamma
        M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
        eigenM <- eigen(M, T)
        YSG <- Y %*% (invSigma * Gamma)

        ## compute likelihood and check for convergence
        old.llh <- llh
        ## Ussing Woodbury matrix identity
        ## log det(Gamma Gamma' + Sigma) =
        ## log [det(invSigma^{-1})det(M)]
        ## tr((Gamma Gamma' + Sigma)^{-1} S) =
        ## tr(invSigma S) - tr(B'B)/n
        ## where B = (M')^{-1/2} YSG'
        logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
        B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
        logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
        llh <- -logdetY - logtrY

        if (abs(llh - old.llh) < tol * abs(llh)) {
            converged <- TRUE
            break
        }
    }

    # GLS to estimate factor loadings
    svd.H <- svd(t(Gamma) %*% (invSigma * Gamma))
    Z <- Y %*% (invSigma * Gamma) %*% (svd.H$u %*% (1/svd.H$d * t(svd.H$v)))

    return(list(Gamma = Gamma,
                Sigma = 1/invSigma,
                Z = Z,
                niter = iter,
                converged = converged))

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
#' @return (i) clusters: list clusters as igraph objects;
#' (ii) fit: list of fitting results for each cluster as a lavaan object;
#' (iii) dfc: data.frame of summary results
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
#' print(clusters$dfc)
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
		 pv1<- 1
		 pv2<- 1
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
	if (map) {
		verbose <- FALSE
	    #plot(graph, layout = l)
		gplot(graph, l = "fdp")
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
				gH <- HM[[x]]
				E(gH)$weight <- 1
				plot(gH, vertex.color = vcol[x], vertex.size = d[V(HM[[x]])$name], 
				 layout = l, main = paste0("Hidden Module ", M[x]))
			Sys.sleep(1)
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
#' columns to graph nodes. If \code{membership} is not NULL, is currently
#' ignored, \code{data = NULL}.
#' @param h Cutting the minimax clustering at height, h = 1 - abs(cor(j,k)),
#' yielding a merged node (and a reduced data set) in which every node in the
#' cluster has correlation of at least cor(j,k) with the prototype node.
#' By default, \code{h = 0.5}, i.e. cor(j,k) = 0.5.
#' @param membership Cluster membership. A vector of cluster membership
#' identifiers as numeric values, where vector names correspond to graph
#' node names. By default, \code{membership = NULL}.
#' @param HM Hidden cluster label. If membership is derived from clusterGraph:
#' HM = "LV", a latent variable (LV) will be defined as common unknown cause
#' acting on cluster nodes. If HM = "CV", cluster nodes will be considered as
#' regressors of a latent composite variable (CV). Finally, if HM = "UV", an
#' unmeasured variable (UV) is defined, where source nodes of the module (i.e.,
#' in-degree = 0) act as common regressors influencing the other nodes
#' via an unmeasured variable. By default, \code{HM = NULL}
#' @param verbose A logical value. If FALSE (default), the merged graphs will
#' not be plotted to screen.
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
#' M <- mergeNodes(G, data = alsData$exprs, h = 0.5, verbose=TRUE)
#'
#' # Gene memberships with EBC method and size=10
#' m <- clusterGraph(G, type = "ebc", size = 10)
#' M <- mergeNodes(G, membership = m, HM = "LV", verbose=TRUE)
#'
#' # Gene memberships defined by user
#' c1 <- c("5894", "5576", "5567", "572", "598")
#' c2 <- c("6788", "84152", "2915", "836", "5530")
#' c3 <- c("5603", "6300", "1432", "5600")
#' m <- c(rep(1,5), rep(2,5), rep(3,4))
#' names(m) <- c(c1, c2, c3)
#' M <- mergeNodes(G, membership = m, HM = "CV", verbose=TRUE)
#'
mergeNodes<- function(graph, data, h=0.5, membership=NULL, HM=NULL, verbose=FALSE, ...)
{
	# Set membership object :
	if (is.numeric(membership)){
	 nodes <- names(membership)
	 if (is.null(HM)) HM <- "HM"
	 membership <- paste0(HM, membership)
	 names(membership) <- nodes
	}else{
	 membership <- prototype(graph, data, h=h, size=3)
	}
	
	LM<- NULL
	for ( i in 1:length(table(membership)) ) {
	 m<- names(table(membership))[i]
	 LMi<- V(graph)$name[which(V(graph)$name %in% names(membership)[membership==m])]
	 LM<- c(LM, list(LMi))
	}
	names(LM)<- names(table(membership))

	# visualize graph object :
	gLM<- as_graphnel(graph_from_edgelist(as_edgelist(graph)))
	for ( i in 1:length(LM) ) {
	 gLMi<- graph::combineNodes(LM[[i]], gLM, names(LM)[i], mean)
	 gLM<- gLMi
	}

	ig<- graph_from_graphnel(gLM)
	if( length(V(ig)$color) == 0 ) V(ig)$color<- "white"
	V(ig)$color[substr(V(ig)$name,1,2) == HM] <- "pink"
	V(ig)$color[substr(V(ig)$name,1,1) == "p"] <- "pink"
	V(ig)$name <- gsub("p", "", V(ig)$name)
	if (verbose) gplot(ig)

	return(list(gLM=ig, membership=membership))
}

prototype<- function(graph, data, h, size, ...)
{
	# Set graph & data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	Y <- data[,nodes] #colnames(Y); head(Y)
	ig <- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	D <- as.dist(1-abs(cor(Y)))

	# hierarchical clustering with prototypes
	hc <- protoclust::protoclust(D)
	plot(hc);abline(h=h, lty=1, col="red")
	cutd <- protoclust::protocut(hc, k=NULL, h=h)
	protos <- hc$labels[cutd$protos]
	cln <- sort(cutd$cl) # as.numeric(cln); names(cln)

	# cluster membership with nodes > size
	nrep <- as.numeric(table(cln))
	cl <- unlist(lapply(1:length(protos), function(x) rep(paste0("p",protos[x]), nrep[x])))
	names(cl) <- names(cln)
	csize <- names(table(cl))[table(cl) >= size]
	cl <- cl[cl %in% csize] # table(cl)

	return(membership = cl)
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
