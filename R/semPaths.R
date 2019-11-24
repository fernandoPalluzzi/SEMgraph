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


SEMmodel <- function(graph, nodes, group, B, ...)
{
	# Set nodes and from-to-matrix representation of gene-gene links
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	ig <- simplify(ig, remove.loops = TRUE)
	ftm <- as_data_frame(ig)
	#head(ftm)
	if (nrow(ftm) == 0) return(list(model = NULL, graph = NULL))

	if(is.directed(ig) & sum(which_mutual(ig)) > 0) {
		sel <- as.numeric(c(E(ig)[which_mutual(ig)]))
		ftm <- as_data_frame(ig)[-sel,]
		ubg <- as.undirected(graph_from_data_frame(as_data_frame(ig)[sel,]))
		ftb <- as_data_frame(ubg)
	} else {
		ftb <- NULL
	}

	# Estimated "common" beta (or psi) coefficients with MLE method:
	if (length(B) == 0) {
		modelY <- modelV <- vector()
		if (is.directed(ig)) {
			for(j in 1:nrow(ftm)) {
				modelY[j] <- paste0("z", ftm[j, 2], "~z", ftm[j, 1])
			}
			if (length(ftb) > 0 ) {
				for(k in 1:nrow(ftb)) {
					modelV[k] <- paste0("z", ftb[k, 2], "~~z", ftb[k, 1])
				}
			}
		} else {
			for(j in 1:nrow(ftm)) {
				modelY[j] <- paste0("z", ftm[j, 2], "~~z", ftm[j, 1])
			}
		}
	}

	# Fixed "common" beta (or psi) coefficients with KEGG attribute:
	if (length(B) > 0) {
		modelY <- modelV <- vector()
		w <- ftm$weight
		if (is.directed(ig)) {
			for(j in 1:nrow(ftm)) {
				modelY[j] <- paste0("z", ftm[j, 2], "~", w[j]*B,
				                    "*z", ftm[j, 1])
			}
			if (length(ftb) > 0) {
				for(k in 1:nrow(ftb)) {
					modelV[k] <- paste0("z", ftb[k, 2], "~~z", ftb[k, 1])
				}
			}
		} else {
			for(j in 1:nrow(ftm)) {
				modelY[j] <- paste0("z", ftm[j, 2], "~~", w[j]*B,
				                    "*z", ftm[j, 1])
			}
		}
	}

	# group mean differences effect
	modelC <- sort(paste0("z", nodes, "~", "group"))
	# equal residual variance
	#modelV <- sort(paste0("z", nodes, "~~v*z", nodes))
	# unequal residual variance
	#modelV <- sort(paste0("z", nodes, "~~z", nodes))
	if( length(group) == 0 ) {
		model <- paste0(c(sort(modelY), modelV))
	} else {
		model <- paste0(c(modelC, sort(modelY), modelV))
	}

	return(list(model = model, graph = ig))
}

#' @title Fit a network as a Structural Equation Model (SEM)
#'
#' @description SEMfit converts a network with directed (x -> y),
#' undirected (x -- y), or mixed (directed: x -> y and bidirected:
#' x <-> y) edges to a SEM, and fits it through the
#' \code{\link[lavaan]{lavaan}} function, via Maximum Likelihood
#' Estimation (estimator = "ML", default estimator in
#' \code{\link[lavaan]{lavOptions}}). If a binary group
#' (i.e., case/control) variable is present, a group effect is added to
#' the model as an exogenous variable. SEMfit can handle loops,
#' although double links are not allowed, such as: self-loops, direct
#' feedbacks, and bow (i.e., double directed/bidirected) links. If a
#' directed or mixed network is given as input, SEMfit outputs three main
#' sets of parameter estimates: (i) group effect estimations on network
#' nodes (gamma), and aggregted group effects: D = sum of group
#' effects, adjusted by the residual variance; A = sum of the tagret
#' nodes perturbation (i.e., group effect) accumulation from source
#' nodes; and E = sum of the source nodes perturbation (i.e., group
#' effect) emission towards target nodes; (ii) beta coefficients (network
#' interaction effects), and (iii) residual (co)variances (psi).
#' When the group variable is present, beta and psi coefficients, common
#' to both groups, will be estimated. If the input is an undirected
#' network, only psi coefficients will be estimated. In the latter case,
#' two sets of parameters will be estimated: (i) group effect estimations
#' on network nodes (gamma), and aggregted group effects: D = sum of group
#' effects, adjusted by residual variances and covariances; V = sum of
#' group effects, adjusted by the residual variances; and C = sum of
#' the between-node covariance perturbation (i.e., group effect); and
#' (ii) variance and covariance common to both groups.
#' P-values for parameter sets are computed by conventional z-test
#' (= estimate/SE), through \code{\link[lavaan]{lavaan}}. However, for
#' very large networks, permutation test P-values can be computed
#' (see \code{\link[flip]{flip}} for details).
#'
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param B Node-node interaction fixed weight. If B is NULL (default),
#' beta coefficients will be estimated by MLE. If B is numeric, it will
#' be used as a scaling factor for the edge weights in the graph object
#' (graph attribute E(graph)$weight). Since SEMgraph scales data before
#' model fitting, we suggest a grid search for the optimal B value in
#' the interval [0, 0.3]. As a rule of thumb, to our experience,
#' B = 0.1 performs well on any network.
#' @param perm Number of permutations. By default, perm is set to 0 and
#' conventional standard errors will be computed. If perm > 1, P-values
#' will be computed from a moment-based chi-squared approximation
#' derived from the empirical distribution of permuted data (Larson and
#' Owen, 2015). To reduce computational time costs per permutation, we
#' suggest perm = 500 (this will leave P-values precision almost
#' unaltered). In perm = 1, no P-values are calculated.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "fit", SEM fitted object of class lavaan;
#' \item "gest", group effect estimates and P-values on subgraph nodes;
#' \item "model", SEM model as a string;
#' \item "graph", the induced subgraph of the input network mapped on
#' data variables;
#' \item "dataXY", input data subset mapping graph nodes, plus group at
#' the first column (if no group is specified, this column will take NA
#' values).
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cor
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom flip flip plot
#' @export
#'
#' @references
#'
#' Yves Rosseel (2012). lavaan: An R Package for Structural Equation
#' Modeling. Journal of Statistical Software, 48(2): 1-36.
#' URL http://www.jstatsoft.org/v48/i02/
#'
#' Pepe D, Grassi M (2014). Investigating perturbed pathway modules
#' from gene expression data via Structural Equation Models. BMC
#' Bioinformatics, 15: 132.
#' URL https://doi.org/10.1186/1471-2105-15-132
#'
#' Larson JL and Owen AB (2015). Moment based gene set tests. BMC
#' Bioinformatics, 16: 132. https://doi.org/10.1186/s12859-015-0571-7
#'
#' @seealso See \code{\link[lavaan]{lavaan}} for SEM fitting, and
#' \code{\link[flip]{flip}} for permutation.
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#'
#' # Conventional standard errors computation
#' start <- Sys.time()
#' fit1 <- SEMfit(graph, data, group, B = NULL, perm = 0)
#' end <- Sys.time()
#' time1 <- end - start
#'
#' # Moment-based chi-squared approximation (Larson and Owen, 2015)
#' start <- Sys.time()
#' fit2 <- SEMfit(graph, data, group, B = NULL, perm = 10000)
#' end <- Sys.time()
#' time2 <- end - start
#'
#' # Execution time comparison
#' time1
#' time2
#'
#' # if perm = 0 (i.e., canonical SE computation - no permutations):
#' # SEM goodness of fit indices
#' fitMeasures(fit1$fit, c("rmsea", "rmsea.pvalue", "srmr"))
#' R2 <- inspect(fit1$fit, "rsquare"); R2; mean(R2)
#' # SEM parameter estimates
#' summary(fit1$fit)
#' parameterEstimates(fit1$fit)[1:10,]
#' # Group effect estimates and aggregated group effects (i.e., D, A, E)
#' fit1$gest
#' pvalue <- fit1$gest$pvalue[-c(1:3)]
#' # Multiple testing correction
#' length(which(p.adjust(pvalue, method = "BH") < 0.05))
#'
#' # if perm > 1 (i.e., moment-based chi-squared approximation):
#' # SEM goodness of fit indices
#' fitMeasures(fit2$fit, c("rmsea", "rmsea.pvalue", "srmr"))
#' R2 <- inspect(fit2$fit, "rsquare"); R2; mean(R2)
#' # SEM parameter estimates
#' summary(fit2$fit)
#' parameterEstimates(fit2$fit)[1:10,]
#' # Group effect estimates and aggregated group effects (i.e., D, A, E)
#' fit2$gest
#' pvalue <- fit2$gest@res$pchisq[-c(1:3)]    # chi-squared P-value
#' # Multiple testing correction
#' length(which(p.adjust(pvalue, method = "BH") < 0.05))
#' # Plotting flip permutation space
#' par(mfrow = c(2, 2))
#' flip::plot(fit2$gest[-c(1:3)])
#' flip::plot(fit2$gest[1, 2])
#' flip::plot(fit2$gest[1, 3])
#' flip::plot(fit2$gest[2, 3])
#'
SEMfit <- function(graph, data, group = NULL, B = NULL, perm = 0, ...)
{
	# Set data from graph:
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- as.matrix(data[, nodes])
	colnames(dataY) <- paste0("z", nodes)
	dataXY <- cbind(group, dataY)
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	if (is.positive.definite(cor(dataXY)[1:p, 1:p])) {
		covXY <- cor(dataXY)[1:p, 1:p]
	} else {
		covXY <- cor.shrink(dataXY, verbose = TRUE)[1:p, 1:p]
	}

	# SEMmodel via from-to-matrix representation of node-node links:
	model <- SEMmodel(graph, nodes, group, B)
	#cat(model)
	if (length(model[[1]]) == 0) {
		cat("NULL SEM fitting: No.edges = 0 !\n\n")
		return(list(fit = NULL, model = NULL, graph = NULL))
	}

	# SEM fitting
	#fit <- sem(model[[1]], sample.cov = covXY, sample.nobs = n, se = "standard", fixed.x = TRUE)
	#int.ov.free: If FALSE, the intercepts of the observed variables are fixed to zero
	#auto.var: If TRUE, the residual variances are set free
	#auto.cov.y: If TRUE, the covariances of dependent variables are set free
	if (perm == 0) SE <- "standard" else SE <- "none"
	fit <- lavaan(model[[1]], sample.cov = covXY,
	              sample.nobs = n,
	              se = SE,
	              fixed.x = FALSE,
	              int.ov.free = TRUE,
	              auto.var = TRUE,
	              auto.cov.y = FALSE)
	if (fit@Fit@converged == TRUE) {
		srmr <- as.numeric(fitMeasures(fit, c("srmr")))
		cat("Model converged:", fit@Fit@converged, "srmr:", srmr, "\n\n")
	} else {
		cat("Model converged:", fit@Fit@converged, "srmr:", NA, "\n\n")
		return(list(fit = NULL, gest = NULL))
	}

	gest <- SEMflip(fit = fit, data = dataXY, group = group, perm = perm)
	if (!is.na(gest[[2]])) {
		cat("P-value of overall node group differences:", gest[[2]], "\n\n")
	}
	if (length(group) == 0) dataXY <- cbind(group = rep(NA, n), dataXY)
	colnames(dataXY) <- gsub("z", "", colnames(dataXY))

	return(list(fit = fit, gest = gest[[1]], model = model[[1]],
	            graph = model[[2]], dataXY = dataXY))
}

#' @title Fit a network as a two-group SEM and test group influence on
#' network edges
#'
#' @description SEMfit2 converts an input network to a SEM and fits it
#' within each group (i.e., two-group SEM). SEMfit2 outputs two sets of
#' parameter estimates for each group: beta coefficients (network
#' interaction effects) and residual (co)variances if the input is a
#' directed or mixed graph (see \code{\link[SEMgraph]{SEMfit}} for
#' details), or variances and covariances if the input is an undirected
#' graph.
#'
#' Instead of modeling group effects as an exogenous variable, SEMfit2
#' estimates the differences of the beta and/or psi coefficients
#' (network edges) between groups.
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects.
#' @param fit Optional argument. In place of the arguments graph, data
#' and group, a fitted model object returned by \code{\link{SEMfit}} or
#' \code{\link{SEMbap}} functions, can be specified.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of two objects: a fitted model object of class lavaan
#' ("fit"), and the estimate group differences ("dest") of the beta
#' coefficients, covariances, or both if the input model is a directed,
#' undirected, or mixed graph, respectively.
#'
#' @importFrom corpcor is.positive.definite cor.shrink
#' @export
#' @references
#' Yves Rosseel (2012). lavaan: An R Package for Structural Equation
#' Modeling. Journal of Statistical Software, 48(2), 1-36.
#' URL http://www.jstatsoft.org/v48/i02/
#'
#' Pepe D, Grassi M (2014). Investigating perturbed pathway modules
#' from gene expression data via Structural Equation Models.
#' BMC Bioinformatics, 15: 132.
#' URL https://doi.org/10.1186/1471-2105-15-132
#'
#' @seealso \code{\link[lavaan]{lavaan}}
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cor pnorm pchisq
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#'
#' # Fitting from network data
#' fit2 <- SEMfit2(graph, data, group)
#' summary(fit2$fit)
#' dest <- fit2$dest
#' # Multiple testing correction
#' length(which(p.adjust(dest[, 7], method = "BH") < 0.05))
#'
#' # Fitting from SEMfit output
#' fit1 <- SEMfit(graph, data, group, B = NULL, perm = 10000)
#' fit2 <- SEMfit2(fit = fit1)
#' summary(fit1$fit)
#' dest <- fit1$dest
#' # Multiple testing correction
#' length(which(p.adjust(dest[, 7], method = "BH") < 0.05))
#'
SEMfit2 <- function(graph, data, group, fit = list(), ...)
{
	# SEMmodel via from-to-matrix representation of gene-gene links
	if (length(fit) == 0) {
		nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
		dataY <- data[, nodes]
		colnames(dataY) <- paste0("z", nodes)
		p <- ncol(dataY)
		#head(dataY)
		model <- SEMmodel(graph, nodes, group = NULL, B = NULL)[[1]]
	} else {
		p <- ncol(fit$dataXY) - 1
		dataY <- fit$dataXY[, -1]
		colnames(dataY) <- paste0("z", colnames(dataY))
		group <- fit$dataXY[, 1]
		model <- fit$model[-c(1:p)]
	}

	# covariances for cases (GROUP 1)
	data1 <- dataY[group == 1,]
	n1 <- nrow(data1)
	if (is.positive.definite(cor(data1)[1:p, 1:p])) {
		cov1 <- cor(data1)[1:p, 1:p]
	} else {
		cov1 <- cor.shrink(data1,verbose = TRUE)[1:p, 1:p]
	}

	# covariances for controls (GROUP 0)
	data0 <- dataY[group == 0,]
	n0 <- nrow(data0)
	if (is.positive.definite(cor(data0)[1:p, 1:p])) {
		cov0 <- cor(data0)[1:p, 1:p]
	} else {
		cov0 <- cor.shrink(data0, verbose = TRUE)[1:p, 1:p]
	}

	# SEM fitting:
	fit <- lavaan(model, sample.cov = list(cov0, cov1),
	              sample.nobs = list(n0, n1),
	              se = "standard",
	              fixed.x = FALSE,
	              int.ov.free = TRUE,
	              auto.var = TRUE,
	              auto.cov.y = FALSE)
	#group.equal=c("residuals","residual.covariances"))
	if (fit@Fit@converged == TRUE) {
		srmr <- fitMeasures(fit, c("srmr"))
		cat("Model converged:", fit@Fit@converged, "srmr:", srmr, "\n\n")
	} else {
		cat("Model converged:", fit@Fit@converged, "srmr:", NA, "\n\n")
		return(list(fit = NULL, dest = NULL))
	}
	#summary(fit)

	# z-test of the beta's (psi's) differences between groups
	est <- parameterEstimates(fit)
	est1 <- est[est$group == 1,][1:length(model),]
	est2 <- est[est$group == 2,][1:length(model),]
	d_est <- est1$est - est2$est
	d_se <- sqrt(est1$se^2 + est2$se^2)
	pvalue <- 2*(1 - pnorm(abs(d_est/d_se)))
	d_lower <- d_est - 1.96*d_se
	d_upper <- d_est + 1.96*d_se

	dest <- cbind(est1[, 1:3], d_est, d_se, d_z = d_est/d_se, pvalue,
	              d_lower, d_upper)
	dest <- data.frame(lapply(dest,
	                   function(y) if(is.numeric(y)) round(y, 3) else y))

	pX2 <- 1 - pchisq(-2*sum(log(dest[, 7])), df = 2*length(dest[, 7]))
	cat("Fisher's combined pvalue of edge group differences:", pX2, "\n\n")

	return(list(fit = fit, dest = dest))
}

#' @title Bow-free graphical lasso data-driven and knowledge-based
#' interaction search
#'
#' @description Fast estimation of sparse inverse covariance using Gaussian
#' Graphical Models (GGMs), with de-sparsified graphical lasso method
#' (D-S_GL; see \code{\link[SILGGM]{SILGGM}}). Significant bow-free covariaces
#' will be selected. Only covariances connecting nodes without direct links
#' in the input graph (i.e., bow-free covariances) and present in the
#' reference interactome (i.e., biologically validated) are added to the
#' input graph. Once new covariances has been added, the resulting mixed
#' network is fitted (see \code{\link[SEMgraph]{SEMfit}} for further details).
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param gnet External interaction network as an igraph object. Interaction
#' data from this network will be used to integrate additional interaction
#' information inside the model. Two preset databases are available:
#' (i) kegg, for KEGG signaling pathways (directed), and (ii) string, for
#' STRING protein interactions (undirected). Direct and indirect interactions
#' will be added on the base of the selected interactome. Note that, if
#' the input graph is undirected, an undirected reference interactome
#' is expected.
#' @param d An integer value indicating the maximum length of indirect
#' interactions between pairs of bow-free nodes. If d = 1, direct
#' interactions between bow-free nodes will be searched in the reference
#' interactome. If d > 1, indirect interactions of length d or shorter
#' (i.e., with at most d - 1 connectors) between bow-free nodes will be
#' searched.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param B Node-node interaction fixed weight. If B is NULL (default),
#' beta coefficients will be estimated by MLE. If B is numeric, it will
#' be used as a scaling factor for the edge weights in the graph object
#' (graph attribute E(graph)$weight). Since SEMgraph scales data before
#' model fitting, we suggest a grid search for the optimal B value in
#' the interval [0, 0.3]. As a rule of thumb, to our experience, B = 0.1
#' performs well on any network.
#' @param perm Number of permutations. By default, perm is set to 0 and
#' conventional standard errors will be computed. If perm > 1, P-values
#' will be computed from a moment-based chi-squared approximation derived
#' from the empirical distribution of permuted data (Larson and Owen, 2015).
#' To reduce computational time costs per permutation, we suggest perm = 500
#' (this will leave P-values precision almost unaltered). In perm = 1,
#' no P-values are calculated.
#' @param alpha Significance level used for GGM search. Alpha is set to
#' 0.05 by default.
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen, saving execution time.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "fit", SEM fitted object of class lavaan;
#' \item "gest", group effect estimates and P-values on subgraph nodes;
#' \item "model", SEM model as a string;
#' \item "graph", list of four graphs:
#' \itemize{
#' \item "ig", the induced subgraph mapped on input data,
#' \item "guu", the undirected GGM graph,
#' \item "gLV", the directed latent GGM graph,
#' \item "Ug", union of graphs "ig" and "guu";
#' }
#' \item "dataXY", input data subset mapping graph nodes, plus group at
#' the first column (if no group is specified, this column will take NA
#' values).
#' }
#'
#' @import igraph
#' @import lavaan
#' @import SILGGM
#' @importFrom stats na.omit var
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom flip flip plot
#' @export
#'
#' @references
#'
#' Zhang R, Ren Z, Chen W (2018). SILGGM: An extensive R package for
#' efficient statistical inference in large-scale gene networks.
#' PLoS Comput. Biol., 14(8): e1006369.
#' https://doi.org/10.1371/journal.pcbi.1006369
#'
#' Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse covariance
#' estimation with the graphical lasso. Biostatistics, 9(3): 432-441.
#' https://doi.org/10.1093/biostatistics/kxm045
#'
#' @seealso \code{\link[lavaan]{lavaan}}, \code{\link{igraph}},
#' \code{\link[SILGGM]{SILGGM}}
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#' # Network degrees of freedom
#' ndf <- vcount(graph)*(vcount(graph) - 1)/2 - ecount(graph)
#' ggm <- SEMggm(graph = graph, data = data, gnet = kegg, group = group,
#'               d = 2,
#'               perm = 10000,
#'               alpha = 1/ndf,
#'               verbose = TRUE)
#'
#' # Results summary
#' summary(ggm$gest)
#' pval <- ggm$gest@res$pchisq[-c(1:3)]
#' length(which(p.adjust(pval, method = "BH") < 0.1))
#'
#' # Plot output graphs
#' par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
#' ig <- ggm$graph$ig; E(ig)$weight <- 1
#' guu <- ggm$graph$guu
#' gLV <- ggm$graph$gLV
#' Ug <- ggm$graph$Ug; E(Ug)$weight <- 1
#' gplot(ig, main = "Input network")
#' plot(guu, main = "Covariance network")
#' plot(gLV, main = "Latent Variable network")
#' gplot(Ug, main = "Input + Covariance network")
#'
SEMbap<- function(graph, data, group=NULL, B=NULL, perm=0, type="dsep", alpha=0.05, method="BH", verbose=FALSE,...)
{
	# set SEM objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY<- data[,nodes] #colnames(dataY); head(dataY)
	ig<- simplify(induced_subgraph(graph, vids= which(V(graph)$name %in% nodes)))
	if( !is_dag(ig) ) cat("WARNING: input graph is not acyclic (DAG) !","\n\n")
	ug<- as.undirected(ig, mode="collapse") #V(ug)$name
	sem<- SEMfit(graph=graph, data=data, group=group, B=NULL, perm=1)

	# Covariance search :
	if(type == "dsep") {
	 SET<- SEMdsep(fit=sem$fit, psi=FALSE, alpha=alpha, method=method)
	}else if (type == "ggm") {
	 SET<- ggm.test(fit=sem, method="D-S_GL", alpha=alpha)
	}else if (type == "pc") {
	 SET<- pc.test(fit=sem, indepTest=pcalg::gaussCItest, alpha=alpha)
	}

	# Covariance and latent variables graphs (guu, gLV)
	guu<- graph_from_data_frame(SET[,1:2],directed=FALSE)
	V(guu)$name<- gsub("z", "", V(guu)$name)
	if( ecount(guu) == 0 ){
	 cat("NULL SEM fitting: No.covariances=0 !","\n\n")
	 return(list(fit=NULL, gest=NULL, model=NULL, graph=NULL, dataXY=NULL))
	}
	ftm<- as_edgelist(guu)
	ftmLV<- NULL
	for (i in 1:nrow(ftm)) ftmLV<- rbind(ftmLV, cbind(rep(paste0("L",i), 2),ftm[i,]))
	gLV<- graph_from_data_frame(ftmLV, directed=TRUE) #V(gLV)$name
	V(gLV)$color<- ifelse(substr(V(gLV)$name,1,1)=="L","yellow","lightblue")
	if(verbose) {
	 plot(guu); Sys.sleep(5)#readline(prompt="Press [enter] to continue")
	 plot(gLV); Sys.sleep(0)#readline(prompt="Press [enter] to continue")
	}
	graph<- list(ig=ig, guu=guu, gLV=gLV)

	# SEM fitting with interactome bow-free covariances:
	if( is.directed(ig) & !is.na(guu)[1] ) guu<- as.directed(guu, mode="mutual")
	Ug<- graph.union(g=list(ig, guu)[!is.na(list(ig, guu))])
	fit1<- SEMfit(graph=Ug, data=dataY, group=group, B=B, perm=perm)
	# summary(fit1$fit)
	graph<- list(ig=ig, guu=guu, gLV=gLV, Ug=Ug)

	return(list(fit=fit1$fit, gest=fit1$gest, model=fit1$model, graph=graph, dataXY=fit1$dataXY))
}

ggm.test<- function(fit, method="D-S_GL", alpha=0.05, ...)
{
	# Large-Scale Gaussian Graphical Model (D-S_GL) on residual data:
	A0<- as_adj(as.undirected(fit$graph), type="both", sparse=FALSE)
	group<- fit$dataXY[,1] #group
	B<- inspect(fit$fit, "est")$beta #colnames(B)
	p<- ncol(B)
	if(sum(group) == 0) {
	 Z<- fit$dataXY[,-1]
	 A0<- A0[colnames(Z), colnames(Z)]
	}else{
	 Z<-  cbind(fit$dataXY[,-1],group) #colnames(Z)
	 A0<- rbind(cbind(A0,rep(1,(p-1))),c(rep(1,(p-1)),0))
	 rownames(A0)[p]<-colnames(A0)[p]<- "group"
	 A0<- A0[colnames(Z), colnames(Z)] # colnames(A0)
	}
	E<- Z%*%(diag(p)- B) #colnames(E)

	wi<- SILGGM::SILGGM(E, method=method, global=TRUE, alpha=alpha)
	A1<- wi$global_decision[[1]] #sum(A1/2)
	rownames(A1)<- colnames(A1)<- colnames(E) #A1[1:10,1:10]

	#guu<- graph_from_adjacency_matrix(A1, mode="undirected");guu
	#as.undirected(fit$graph) %s% guu
	A1<- ifelse(A1 == A0, 0, A1) #sum(A1/2)
	del<- which(colSums(A1)==0)
	if( length(del) >  0 ) A1<- A1[-del,-del] #A1[1:10,1:10]
	SET<- as_data_frame(graph_from_adjacency_matrix(A1, mode="undirected"))
	return( SET )
}

pc.test<- function(fit, indepTest=pcalg::gaussCItest, alpha=0.05, ...)
{
	# PCalg skeleton search on residual correlations :
	A0<- as_adj(as.undirected(fit$graph), type="both", sparse=FALSE)
	group<- fit$dataXY[,1] #group
	B<- inspect(fit$fit, "est")$beta #colnames(B)
	p<- ncol(B)
	if(sum(group) == 0) {
	 Z<- fit$dataXY[,-1]
	 A0<- ifelse(A0 == 1, TRUE, FALSE)[colnames(Z), colnames(Z)]
	}else{
	 Z<-  cbind(fit$dataXY[,-1],group) #colnames(Z)
	 A0<- rbind(cbind(A0,rep(1,(p-1))),c(rep(1,(p-1)),0))
	 rownames(A0)[p]<-colnames(A0)[p]<- "group"
	 A0<- ifelse(A0 == 1, TRUE, FALSE)[colnames(Z), colnames(Z)]
	}
	E<- Z%*%(diag(p)- B) #colnames(E)

	ske<- pcalg::skeleton(suffStat=list(C=cor(E), n=nrow(E)),
		  indepTest=indepTest, alpha=alpha, labels=colnames(E),
		  method="stable.fast", fixedGaps=A0, fixedEdges=NULL, numCores=8)
	A1<- as(ske@graph, "matrix")

	del<- which(colSums(A1)==0)
	if( length(del) >  0 ) A1<- A1[-del,-del] #A1[1:10,1:10]
	SET<- as_data_frame(graph_from_adjacency_matrix(A1, mode="undirected"))
	return( SET )
}

SEMflip <- function(fit, data, group, perm, ...)
{
	# Set SEM objects:
	if (length(group) == 0 || perm == 1) {
		return(list(gest = NULL, pval = NA))
	}
	q <- ncol(data)
	p <- q - 1
	B <- inspect(fit, "est")$beta[-q, -q]
	#B[1:10, 1:10]
	psi <- inspect(fit, "est")$psi[-q, -q]
	#psi[1:10, 1:10]
	Y <- data[, colnames(B)]
	#Y[1:10, 1:10]

	if (sum(B) != 0) {
		Q <- t(diag(p) - B)%*%diag(1/sqrt(diag(psi)))
		D <- as.matrix(Y)%*%Q%*%rep(1/sqrt(p), p)
		#D[1:10]; diff
		A <- as.matrix(Y)%*%B%*%rep(1/sqrt(p), p)
		#A[1:10]; sink
		E <- as.matrix(Y)%*%t(B)%*%rep(1/sqrt(p), p)
		#E[1:10]; source
		df <- data.frame(D, A, E, group)
		#apply(df, 2, var)
		model <- paste0("D~group;A~group;E~group")
	} else if (sum(B) == 0) {
		E <- eigen(psi)    # Eigenvalues and eigenvectors of psi
		Q <- E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
		D <- as.matrix(Y)%*%Q%*%rep(1/sqrt(p), p)
		#D[1:10]; diff
		C <- as.matrix(Y)%*%(psi-diag(diag(psi)))%*%rep(1/sqrt(p), p)
		#A[1:10]; cov
		V <- as.matrix(Y)%*%diag(1/diag(psi))%*%rep(1/sqrt(p), p)
		#E[1:10]; var
		df <- data.frame(D, C, V, group)
		#apply(df, 2, var)
		model <- paste0("D~group;C~group;V~group")
	}

	if (perm > 1) {
		# permutation pvalues
		Z <- as.matrix(Y)%*%t((diag(p) - B))
		#head(Z)
		perm <- list(B = perm, seed = 123)
		gest <- flip::flip(cbind(df[, -4], Z), ~group,
		                   perms = perm,
		                   statTest = "coeff")
		# Chi-square approximate pvalues
		aveT2 <- apply(gest@permT[-1,]^2, 2, mean)
		varT2 <- apply(gest@permT[-1,]^2, 2, var)
		ni <- 2*(aveT2^2/varT2)
		T2 <- (ni*gest@permT[1,]^2)/aveT2
		gest@res$pchisq <- 1 - stats::pchisq(T2, ni)
		# Fisher's combined pvalue
		#pval <- flip::npc(gest[-c(1:3)], "fisher")@res$p
		pval <- gest@res$pchisq[1]
	} else if (perm == 0) {
		# Wald's pvalues
		res1 <- parameterEstimates(fit)[1:p, 1:7]
		fit2 <- sem(model, data = df, se = "standard", fixed.x = TRUE)
		res2 <- parameterEstimates(fit2)[1:3, 1:7]
		#summary(fit2)
		gest <- rbind(res2, res1)
		# Fisher's combined pvalue
		#pval <- 1 - pchisq(-2*sum(log(res1$pvalue)), 2*length(res1$pvalue))
		pval <- gest$pvalue[1]
	}

	return(list(gest = gest, pval = pval))
}

#' @title Shortest path and tree-based network reduction methods
#'
#' @description Uses different shortest path and tree-based strategies
#' for reducing a network to an optimal subset of nodes and edges.
#' @param graph An igraph object.
#' @param type Network reduction method. If subgraph = "mst", the input
#' network will be reduced to a MST, using the \code{\link[igraph]{mst}}
#' igraph implementation. This will produce a tree with the same node as
#' the input graph (N), and N-1 edges, such that the sum of edges is the
#' minimum among all the possible subgraphs. If subgraph = "core", an
#' induced subgraph of the input seed set will be generated (see
#' \code{\link[SEMgraph]{seedweight}} for seed list creation).
#' If subgraph = "kou", the fast Steiner tree algorithm from Kou et al. (1981)
#' will be applied starting from a set of seed nodes. If subgraph = "usp",
#' the resulting network will be the union of the shortest paths.
#' If subgraph = "isn" a subnetwork will be induced from the shortest
#' path nodes.
#' @param seed A vector containing seed node identifiers, either user-defined
#' or yielded by \code{\link[SEMgraph]{seedweight}}.
#' @param eweight Edge weight type derived from
#' \code{\link[SEMgraph]{edgeweight.cfa}},
#' \code{\link[SEMgraph]{edgeweight.r2z}},
#' \code{\link[SEMgraph]{edgeweight.sem}},
#' or from user-defined distances. This option determines the
#' weight-to-distance transform used to reduce the input graph. If set to
#' "none" (default), edge weights will be set to 1. If eweight = "kegg",
#' repressing interactions (-1) will be set to 1 (maximum distance),
#' neutral interactions (0) will be set to 0.5, and activating
#' interactions (+1) will be set to 0 (minimum distance).
#' If eweight = "zsign", all significant interactions will be set to 0
#' (minimum distance), while non-significant ones will be set to 1.
#' If eweight = "pvalue", weights (p-values) will be transformed to the
#' inverse of negative base-10 logarithm. If eweight = "custom", the
#' algorithm will use the distance measure specified by the user as "weight"
#' edge attribute.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A reduced graph of class \code{\link{igraph}}.
#'
#' @import igraph
#' @export
#'
#' @examples
#' # Data loading
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#'
#' # Graph weighting
#' graph1 <- edgeweight.sem(graph = graph, data = data, group = group)
#'
#' # Graph reduction
#' G1 <- reduceGraph(graph = graph1, type = "mst", seed = "none",
#'                   eweight = "pvalue")
#' gplot(G1)
#'
reduceGraph <- function(graph, type, seed = "none", eweight = "none", ...)
{
	# Define positive weights and the graph reduction method
	ig <- graph
	#V(ig)$name
	#E(ig)$weight
	if (eweight == "kegg") edgew <- (1 - E(ig)$weight)/2
	else if (eweight == "zsign") edgew <- 1 - abs(E(ig)$zsign)
	else if (eweight == "pvalue") edgew <- 1/(-log10(E(ig)$pv))
	else if (eweight == "none") edgew <- rep(1, ecount(ig))
	else if (eweight == "custom") edgew <- E(ig)$weight

	if (seed == "pvlm") seed <- V(ig)$name[V(ig)$spvlm == 1]
	else if (seed == "mst") seed <- V(ig)$name[V(ig)$smst == 1]
	else if (seed == "proto") seed <- V(ig)$name[V(ig)$sprot == 1]
	else if (seed == "none") seed <- vector()

	if (type == "mst" & length(seed) == 0) {
		return(minimum.spanning.tree(ig, weights = edgew, algorithm = "prim"))
	}
	else if (length(seed) != 0) {
		if (type == "core" & length(seed) != 0) {
			return(induced_subgraph(graph, seed))
		}
		else if (type == "kou" & length(seed) != 0) {
			return(SteinerTree(graph, seed, eweight = edgew))
		}
		else if (type == "isn" & length(seed) != 0) {
			return(USPT(graph, seed, eweight = edgew, test = "fisher",
			            alpha = 1)[[1]])
		}
		else if (type == "usp" & length(seed) != 0) {
			return(USPT(graph, seed, eweight = edgew, test = "fisher",
			            alpha = 1)[[2]])
		}
	} else {
		cat("WARNING: seed list or seed extraction method required with",
		    type, "method")
	}
}

#' @title Steiner tree via Kou's algorithm
#'
#' @description Find the Steiner tree connecting a set of seed nodes,
#' using the shortest path heuristic from Kou et al. (1981).
#' @param graph An igraph object.
#' @param seed A vector containing seed node identifiers, either user-defined
#' or yielded by \code{\link[SEMgraph]{seedweight}}.
#' @param eweight Edge weight type derived from
#' \code{\link[SEMgraph]{edgeweight.cfa}},
#' \code{\link[SEMgraph]{edgeweight.r2z}},
#' \code{\link[SEMgraph]{edgeweight.sem}},
#' or from user-defined distances. This option determines the
#' weight-to-distance transform used to reduce the input graph. If set
#' to "none" (default), edge weights will be set to 1.
#' If eweight = "kegg", repressing interactions (-1) will be set to 1
#' (maximum distance), neutral interactions (0) will be set to 0.5, and
#' activating interactions (+1) will be set to 0 (minimum distance).
#' If eweight = "zsign", all significant interactions will be set to 0
#' (minimum distance), while non-significant ones will be set to 1.
#' If eweight = "pvalue", weights (p-values) will be transformed to the
#' inverse of negative base-10 logarithm. If eweight = "custom", the
#' algorithm will use the distance measure specified by the user as "weight"
#' edge attribute.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return The Steiner tree as an object of class igraph.
#'
#' @import igraph
#' @export
#' @references
#' Kou L, Markowsky G, Berman L (1981). A fast algorithm for Steiner trees.
#' Acta Informatica, 15(2): 141-145. https://doi.org/10.1007/BF00288961
#'
#' Palluzzi F, Ferrari R, Graziano F, Novelli V, Rossi G, Galimberti D,
#' Rainero I, Benussi L, Nacmias B, Bruni AC, Cusi D, Salvi E, Borroni B,
#' Grassi M. (2017). A novel network analysis approach reveals DNA damage,
#' oxidative stress and calcium/cAMP homeostasis-associated biomarkers
#' in frontotemporal dementia. PLoS ONE 12(10): e0185797.
#' https://doi.org/10.1371/journal.pone.0185797
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' # Return graph properties, take the largest component, and convert
#' # grapNEL to igraph
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' # Transpose data matrix: 32 subjectx (rows) x 19726 genes (columns)
#' data <- t(FTLDu_GSE13162)
#'
#' # Generating edge weights for the KEGG interactome
#' graph1 <- edgeweight.r2z(graph = kegg, data = data, group = group)
#' graph1
#'
#' # Generating seed weights
#' graph2 <- seedweight(graph = graph1, data = data, group=group, alpha = 5E-05)
#' graph2
#' seed <- V(graph2)$name[V(graph2)$spvlm == 1]
#' length(seed)
#'
#' # Using -log10(P-values) as edge weights
#' eweight <- 1/(-log10(E(graph2)$pv))
#' graph3 <- SteinerTree(graph = graph2, seed = seed, eweight = eweight)
#' graph3
#' V(graph3)$color <- ifelse(V(graph3)$name %in% seed, "lightblue", "yellow")
#' V(graph3)$seed <- ifelse(V(graph3)$name %in% seed, 1, 0)
#'
#' # KEGG-derived Steiner tree
#' gplot(graph3, main = "Extracted KEGG Steiner tree")
#' iplot(graph3)    # interactive plot
#'
SteinerTree <- function(graph, seed, eweight)
{
	# Define graph, edge weights, and distance matrix
	ig <- graph
	E(ig)$weight <- eweight
	#E(ig)$zsign
	#E(ig)$pv
	D <- igraph::distances(ig, v = seed, to = seed,
	                       mode = "all",
	                       weights = eweight)
	#D[1:10,1:10]
	#dim(D)

	# Step 1: complete undirected distance graph Gd for terminal nodes
	Gd <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE)
	Gd <- Gd - igraph::edges(E(Gd)[which(E(Gd)$weight == Inf)])
	#plot(as_graphnel(Gd)); Gd
	#E(Gd)$weight

	# Step 2: MST T1 of the complete distance graph Gd
	T1 <- minimum.spanning.tree(Gd, weights = NULL, algorithm = "prim")
	#plot(as_graphnel(T1))

	# Step 3: for each edge in T1, replace it with the shortest path in ig
	edge_list <- as_edgelist(T1)
	N <- nrow(edge_list)
	subgraph <- vector()

	for (n in 1:N) {
		# print(paste("i = ", n, "/", N, " shortest path", sep = ""))
		i <- edge_list[n, 1]
		j <- edge_list[n, 2]

		# Extract from ig all nodes of the shortest paths between edges of T1
		path <- shortest_paths(ig, from = V(ig)[i], to = V(ig)[j],
		                       mode = "all",
		                       weights = eweight,
		                       output = "both")
		vpath <- V(ig)$name[path$vpath[[1]]]
		subgraph <- igraph::union(subgraph, vpath)
	}

	# Step 4: MST Ts of the extracted (induced) sub-graph Gs of ig
	Gs <- induced_subgraph(ig, unique(subgraph))
	#plot(as_graphnel(Gs)); Gs
	#E(Gs)$weight
	Ts <- minimum.spanning.tree(Gs, weights = NULL, algorithm = "prim")
	#plot(as_graphnel(Ts)); Ts

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
	#plot(as_graphnel(St)); St
	#E(St)$weight

	return(St)
}

#' @title Union of the Shortest Path Tree (USPT) algorithm
#'
#' @description Generate a subnetwork as the union of the most significant
#' shortest paths from an input network. Shortest paths are chosen on
#' the base of the perturbation (group difference) of their connections.
#' @param graph An igraph object.
#' @param seed A vector containing seed node identifiers, either user-defined
#' or yielded by \code{\link[SEMgraph]{seedweight}}.
#' @param eweight Edge weight type derived from
#' \code{\link[SEMgraph]{edgeweight.cfa}},
#' \code{\link[SEMgraph]{edgeweight.r2z}},
#' \code{\link[SEMgraph]{edgeweight.sem}},
#' or from user-defined distances. This option determines the
#' weight-to-distance transform used to reduce the input graph. If set
#' to "none" (default), edge weights will be set to 1. If eweight = "kegg",
#' repressing interactions (-1) will be set to 1 (maximum distance), neutral
#' interactions (0) will be set to 0.5, and activating interactions (+1)
#' will be set to 0 (minimum distance).
#' If eweight = "zsign", all significant interactions will be set to 0
#' (minimum distance), while non-significant ones will be set to 1.
#' If eweight = "pvalue", weights (p-values) will be transformed to the
#' inverse of negative base-10 logarithm. If eweight = "custom", the
#' algorithm will use the distance measure specified by the user as "weight"
#' edge attribute.
#' @param alpha Significance level to assess shortest paths significance.
#' By default, alpha is set to 1 (i.e., all the shortest paths will be
#' retained).
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of two igraph objects, namely "Gs" and "sG".
#' The former is the induced subnetwork of the shortest path nodes.
#' The latter is the union of the shortest paths (i.e., the USPT).
#'
#' @import igraph
#' @importFrom stats pchisq
#' @export
#'
#' @references
#' Meier J, Tewarie P, Van Mieghem P (2015). The Union of Shortest Path
#' Trees of Functional Brain Networks. Brain Connectivity 5(9); 575-581.
#' https://doi.org/10.1089/brain.2014.0330
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#'
#' # Generating edge weights
#' graph1 <- edgeweight.sem(graph = graph, data = data, group = group)
#'
#' # Generating seed weights
#' graph2 <- seedweight(graph = graph1, data = data, group = group)
#' seed <- V(graph2)$name[V(graph2)$spvlm == 1]  # select seed nodes
#'
#' # Using -log10(P-values) as edge weights
#' edgeweight <- 1/(-log10(E(graph2)$pv))
#' U <- USPT(graph = graph2, seed = seed, eweight = edgeweight)
#' graph3 <- U[[1]]    # induced subnetwork of the shortest path nodes
#' graph4 <- U[[2]]    # union of the shortest paths tree (USPT)
#' # Plotting output graph
#' V(graph3)$color <- ifelse(V(graph3)$name %in% seed, "yellow", "lightblue")
#' gplot(graph3)
#' V(graph4)$color <- ifelse(V(graph4)$name %in% seed, "yellow", "lightblue")
#' gplot(graph4)
#'
USPT <- function(graph, seed, eweight, alpha = 1, ...)
{
	# Define graph, edge weights, and distance matrix:
	ig <- graph #E(ig)$weight; E(ig)$zsign; E(ig)$pv
	E(ig)$weight <- eweight
	D <- igraph::distances(ig, v = seed, to = seed,
	                       mode = "all",
	                       weights = eweight)
	#D[1:10,1:10]; dim(D)

	# Complete directed distance graph for terminal nodes:
	Gd <- graph_from_adjacency_matrix(D, mode = "undirected",
	                                  weighted = TRUE)
	Gd <- Gd - igraph::edges(E(Gd)[which(E(Gd)$weight == Inf)])

	# For each edge in Gd, replace it with the shortest path:
	edge_list <- as_edgelist(Gd)
	N <- nrow(edge_list)
	if (alpha == 1) alpha <- N
	ftm <- NULL
	subgraph <- vector()

	for (n in 1:N) { #n=1
		#print(paste("i= ",n,"/",N, " shortest path", sep=""))
		i <- edge_list[n, 1]
		j <- edge_list[n, 2]

		# extract from ig all nodes of the shortest paths between edges of Gd
		path <- shortest_paths(ig, from = V(ig)[i], to = V(ig)[j],
		                       mode = "all",
		                       weights = eweight,
		                       output = "both")
		vpath <- V(ig)$name[path$vpath[[1]]]
		pvalue <- E(ig)$pv[path$epath[[1]]]

		# significance test of the shortest path
		ppath <- 1 - pchisq(-2*sum(log(pvalue)), df = 2*length(pvalue))
		#if(test == "fisher") ppath <- 1-pchisq(-2*sum(log(pvalue)), df = 2*length(pvalue))
		#if(test == "brown") ppath <- Brown.test(y = data[, vpath], p = pvalue)
		#if(test == "wald") ppath <- sp.test(vpath, data, group)[[3]]
		ppath[is.na(ppath)] <- 0.5

		if(ppath < alpha/N) {
			for (k in 1:(length(vpath)-1)) {
				ftm <- rbind(ftm, c(vpath[k], vpath[k+1]))
			}
			subgraph <- igraph::union(subgraph, vpath)
		} else {
			ftm <- ftm
			subgraph <- subgraph
		}
	}

	# Merging the shortest paths
	if(length(unique(subgraph)) > 0) {
		sG <- induced_subgraph(ig, unique(subgraph)) #sG
		del <- which(duplicated(ftm) == TRUE)
		if(length(del) > 0) ftm <- ftm[-del,]
		Gs <- simplify(graph_from_edgelist(ftm, directed = FALSE)) #Gs

		if(is.directed(ig)) Gs <- arrowDirection(ug = Gs, dg = ig) #Gs
		# plot(as_graphnel(sG))
		# plot(as_graphnel(Gs))
	} else {
		sG <- Gs <- make_empty_graph(0)
	}

	return(list(sG, Gs))
}

#' @title Compute Average Causal Effect (ACE) for a given DAG
#'
#' @description Compute total effects as ACEs of source variables X on
#' target variables in a DAG (i.e., when paths connecting X to Y are
#' directed and acyclic). The ACE will be estimated as the path coefficient
#' of X (i.e., theta) in the linear equation Y ~ X + pa(X), where pa(X)
#' is the set of parent variables of X (Pearl, 1998).
#' @param model A directed acyclic graph as an igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "theta", data.frame containing total effect estimates;
#' \item "gD", graph representation of significant total effects.
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom dagitty adjustmentSets
#' @export
#'
#' @references
#'
#' Pearl J (1998). Graphs, Causality, and Structural Equation Models.
#' Sociological Methods & Research. https://doi.org/10.1177/0049124198027002004
#'
#' @seealso \code{\link[SEMgraph]{SEMfit}}
#'
#' @examples
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#' group <- c(rep(0, 17), rep(1, 15))
#' ace <- SEMace(graph, data, group)
#' ace$theta
#'
SEMace <- function(graph, data, group = NULL, ...)
{
	# Set igraph and dagitty graph objects
	nodes <- colnames(data)
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	dag <- graph2dagitty(ig)
	#cat(dag)

	# Set distance matrix and distance graph
	D <- igraph::distances(ig, mode = "out", weights = NA)
	D <- ifelse(D == Inf, 0, D)
	#sum(D>0)
	gD <- simplify(graph_from_adjacency_matrix(D, mode = "directed",
	               weighted = TRUE))
	# gD; table(E(gD)$weight

	# Compute total effect (theta = ACE) from DAG
	theta <- NULL
	ftm <- as_edgelist(gD)

	for (i in 1:nrow(ftm)) { #i=482
		cat("\rACE =", i, "of", nrow(ftm), " ")
		flush.console()
		x <- ftm[i, 1]
		y <- ftm[i, 2]
		Z <- dagitty::adjustmentSets(dag, x, y,
		                             type = "minimal",
		                             effect = "total")
		#str(Z)
		#Z <- gRbase::parents(x, as_graphnel(ig))
		#Z <- igraph::neighbors(ig, x, mode = "in")

		ftz <- NULL
		if (length(as.character(unlist(Z))) > 0) {
			for(k in 1:length(Z[[1]])) {
				ftz <- rbind(ftz, c(Z[[1]][k], y))
			}
			df <- data[, c(x, y, Z[[1]])]
			#head(df)
		} else {
			df <- data[, c(x, y)]
		}
		ig1 <- graph_from_data_frame(rbind(c(x, y), ftz))

		# SEM fitting
		if (length(group) == 0) {
			try(quiet(sem <- SEMfit(graph = ig1, data = df, group = NULL)$fit))
			try(t <- parameterEstimates(sem)[1,])
		} else {
			try(quiet(sem <- SEMfit2(graph = ig1, data = df, group = group)))
			try(t <- sem$dest[1,])
		}
		if (!is.null(t)) {
			t$op <- "<-"
			colnames(t)[1:3] <- c("sink", "op", "source")
			theta <- rbind(theta, t)
		}
	}

	return(list(theta = theta, gD = gD))
}

#' @title Test (shortest) paths between two nodes of a network
#'
#' @description Find and fit all the (shortest) paths between two nodes
#' of a network.
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param from Starting node name.
#' @param to Ending node name.
#' @param path If path = "all", all the paths between the two nodes will
#' be included in the fitted model. If path = "short", only shortest paths
#' will be considered.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A list of two objects: a fitted model object of class
#' \code{\link[lavaan]{lavaan}} ("fit"), and the extracted subnetwork as
#' an igraph object ("graph").
#'
#' @import igraph
#' @import lavaan
#' @importFrom dagitty paths
#' @export
#'
#' @examples
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#' group <- c(rep(0, 17), rep(1, 15))
#' spt <- SEMpath(graph, data, group,
#'                from = "5155",
#'                to = "7082",
#'                path = "short")
#'
SEMpath <- function(graph, data, group, from, to, path, ...)
{
	# Set igraph and dagitty graph objects
	nodes <- colnames(data)
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	dag <- graph2dagitty(ig)
	if (distances(ig, from, to, mode = "out", weights = NA) == Inf ) {
		cat("Inf distance from =", from, "to =", to, "! \n\n")
		return(NULL)
	}

	if( path == "short" ) {
		# Set shortest path graph
		paths <- all_shortest_paths(ig, from, to, mode = "out", weights = NA)
		#path <- shortest_paths(ig, from, to,
		#                       mode = "out",
		#                       weights = NULL,
		#                       output = "both")
		nodes <- unique(names(unlist(paths$res)))
	} else if (path == "all") {
		# Set directed path graph
		paths <- dagitty::paths(dag, from, to, directed = TRUE)$paths
		nodes <- unique(unlist(strsplit(gsub("->", "", paths), "  ")))
	}

	# Plot selected paths
	ig1 <- induced_subgraph(ig, vids = which(V(ig)$name %in% nodes))
	gplot(x = ig1, main = "")
	Sys.sleep(5)
	V(ig)$color <- "white"
	V(ig)$color[V(ig)$name %in% V(ig1)$name] <- "aquamarine2"
	E(ig)$color <- "black"
	E(ig)$color[attr(E(ig), "vnames") %in% attr(E(ig1), "vnames")] <- "aquamarine2"
	E(ig)$width <- 1
	E(ig)$width[attr(E(ig), "vnames") %in% attr(E(ig1), "vnames")] <- 2.5
	gplot(x = ig, main = "")

	# SEM fitting
	if (length(group) == 0) {
		sem <- SEMfit(graph = ig1, data = data, group = NULL)
		print(parameterEstimates(sem$fit)[1:length(sem$model),])
	} else {
		sem <- SEMfit2(graph = ig1, data = data, group = group)
		print(sem$dest)
	}

	return(list(fit = sem$fit, graph = ig1))
}

#' @title Search and test for causal relationships in a graph model
#'
#' @description Test the presence of causal relationships in a graphical
#' model. The set of conditionally-independent variables are computed
#' according to either Shipley's (Shipley, ...) or Pearl's (Pearl, ...)
#' d-separation test, for directed acyclic graphs (DAGs) and non-DAGs,
#' respectively. SEM modification indices (MIs) can be used to assign
#' directions (i.e., causality) to the significantly interacting variables
#' found in the first step.
#' @param fit A SEM fitted model object of class \code{\link[lavaan]{lavaan}}.
#' @param psi Logical value. If FALSE (default), covariances will be excluded
#' from the output graph. If TRUE, covariances will be converted into
#' bidirected graph edges.
#' @param k Numeric. Maximum number of conditioning variables in the
#' dependence set. By default, k is set to the degrees of freedom of the
#' graphical model. By default, k is set to 3. Increasing k may significantly
#' increase computational time.
#' @param MI Logical value. If TRUE, new causal relationships are filtered
#' by significant modification indices. By default, MI = FALSE. Enabling MIs
#' will sensibly increase computational time.
#' @param alpha Numeric. Significance level for edge inclusion.
#' @param method P-value correction method. It may take the p.adjust.methods
#' from \code{\link[stats]{p.adjust}} function. By default, it is set to
#' "BH" ("Benjamini-Hochberg").
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A data.frame of causal interactions between d-separated variables.
#'
#' @import igraph
#' @import lavaan
#' @importFrom dagitty lavaanToGraph impliedConditionalIndependencies localTests
#' @importFrom stats qchisq pchisq p.adjust
#' @export
#'
#' @references
#'
#' REF1
#'
#' REF2
#'
#' @examples
#' ...
#'
SEMdsep <- function(fit, psi = FALSE, k = 3, MI = FALSE, alpha = 0.05, method = "BH", ...)
{
	# Set parameters:
	S<- lavInspect(fit, "sampstat")$cov
	n<- lavInspect(fit, "nobs")
	bap <- dagitty::lavaanToGraph(x=parTable(fit), fixed.x=FALSE) #cat(bap)
	#plot( dagitty::graphLayout(bap) )
	dag <- lavaan2graph(toString(bap, "lavaan"), psi = psi) #dag

	if( !is_dag(dag) ){
	 imp<- Filter(function(x) length(x$Z)<=k, dagitty::impliedConditionalIndependencies(bap))
	 res<- dagitty::localTests(bap, type="cis", tests=imp, sample.cov=S, sample.nobs=n)
	 dsep<- NULL
	 for(j in 1:nrow(res)) { #j=1
	  s<- strsplit(rownames(res)[j],"\\_") #s
	  X<- gsub(" ","",s[[1]][1])
	  Y<- gsub(" ","", strsplit( s[[1]][3],"\\|" )[[1]][1] )
	  SET<- gsub(" ","", strsplit( s[[1]][3],"\\|" )[[1]][2] )
	  dsep<- rbind(dsep, data.frame(X, Y, SET, pvalue=res$p.value[j]))
	 }
	}else{ dsep<- dsep.test(dag=dag, S=S, n=n) }

	d_sep<- subset(dsep, p.adjust(dsep$pvalue, method=method) < alpha)
	if( nrow(d_sep) == 0 ) {
	 cat("ALL adjusted pvalues >", alpha, "!", "\n\n")
	 return( dsep=NULL )
	}

	cat("n. significant local test =", nrow(d_sep), "of", nrow(dsep), "\n\n")
	if( MI==TRUE ) return( dsep=MI.test(fit=fit, dsep=d_sep, alpha=alpha, method=method) )

	return( dsep=d_sep )
}

dsep.test<- function(dag, S, n, ...)
{
 	# Basis set of a DAG
	idx <- as.numeric(topo_sort(dag, mode="out"))
	#idx <- unlist(dagitty::topologicalOrdering(dag))
	A <- as_adj(dag, sparse=FALSE)[idx,idx]
	nodes <- colnames(A)

	B <- NULL
	for (r in 1:length(nodes)) {
	 for (s in r:length(nodes)) {
	  if ((A[r,s] != 0) | (s == r)) next
	   else{
		ed <- nodes[c(r,s)]
		pa.r <- nodes[A[,r] == 1]
		pa.s <- nodes[A[,s] == 1]
		dsep <- union(pa.r, pa.s)
		dsep <- setdiff(dsep, ed)
		b <- list(c(ed, dsep))
		B <- c(B, b)
	  }
	 }
	}

	SET <- NULL
	for (i in 1:length(B)) { #i=1
	 cat("\r","BasiSet=", i, "of", length(B))
     flush.console()
	 if (length(B[[i]]) > (n-3)) next
	 X <- B[[i]][1]
	 Y <- B[[i]][2]
	 Set <- paste(B[[i]][-c(1:2)], collapse=",")
	 try(pvalue<- pcor.test(S, B[[i]], n, H0=0.05))
	 SET <- rbind(SET, data.frame(X, Y, SET=Set, pvalue))
	}

	cat("\n done! \n")
	return( SET )
}

MI.test<- function(fit, dsep, alpha=0.05, method="BH", ...)
{
	MI<- modindices(fit)[,1:4]
	pvalue<- sapply(MI$mi, function(x) pchisq(x, df=1, lower.tail=FALSE))
	MI<- subset(MI, p.adjust(pvalue, method=method) < alpha )
	mi<- paste0(MI$lhs,MI$op,MI$rhs) #mi

	dxy<- paste0(dsep$X,"~",dsep$Y) #dxy
	dyx<- paste0(dsep$Y,"~",dsep$X) #dyx
	arrowL<- rep("<->", nrow(dsep))
	arrowL[dxy %in% mi]<- " <- "
	arrowR<- rep("<->", nrow(dsep))
	arrowR[dyx %in% mi]<- " -> "

	return( dsep=cbind(d_sep, arrowL, arrowR) )
}

pcor.test<- function(S, B, n, H0=0, ...)
{
	k <- solve(S[B,B])
    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
	q <- length(B)-2
	if( H0 == 0 ) {
	 df <- n - 2 - q
	 tval <- r * sqrt(df)/sqrt(1 - r * r)
	 pval <- 2 * pt(-abs(tvalue), df)
	}else{
	 z <- atanh(r)
	 se <-  1/sqrt(n - 3 - q)
	 pval <- pchisq((z/se)^2, df=1, ncp=(atanh(.05)/se)^2, lower.tail=FALSE)
	}
	return( pval )
}

#' @title SEM-based gene set analysis
#'
#' @description Test group effect on biological pathways as SEM, evaluating
#' overall pathway perturbation, perturbation emission from source nodes,
#' and perturbation accumulation on target nodes.
#' @param g A list of pathways to be tested.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes.
#' @param group A binary vector. This vector must be as long as the number
#' of subjects. Each vector element must be 1 for cases and 0 for control
#' subjects.
#' @param alpha Gene set test significance level. Alpha is set to 0.05
#' by default.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return A data.frame reporting group effect estimates over nodes
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats pchisq p.adjust na.omit
#' @export
#' @references
#'
#' Grassi M, Pepe D, Palluzzi F (in preparation). Structural Equation
#' Models (SEM)-based assessment of signalling pathway topology and
#' perturbation. Journal ---- (xxxx)
#'
#' @seealso \code{\link[SEMgraph]{SEMfit}}
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' graph <- properties(kegg.pathways$hsa04540_Gap_junction)[[1]]
#' data <- t(FTLDu_GSE13162)
#'
#' ids <- which(names(kegg.pathways) %in% c("hsa04010_MAPK_signaling_pathway",
#'                                          "hsa04310_Wnt_signaling_pathway",
#'                                          "hsa04330_Notch_signaling_pathway",
#'                                          "hsa04722_Neurotrophin_signaling_pathway",
#'                                          "hsa04141_Protein_processing_in_endoplasmic_reticulum",
#'                                          "hsa04144_Endocytosis"))
#'
#' gsa <- SEMgsa(kegg.pathways[ids], data, group, alpha = 0.05)
#' gsa
#'
SEMgsa <- function(g = list(), data, group, alpha, ...)
{
	# Set SEM objects:
	pX2 <- function(x) 1 - stats::pchisq(-2*sum(log(x)), 2*length(x))
	gs <- names(g)
	K <- length(g)
	res.tbl <- NULL
	#DE <- NULL

	for (k in 1:K) {
		cat( "k =", k, gs[k], "\n" )
		if(!is_igraph(g[[k]])) {
			ig <- graph_from_graphnel(g[[k]])
		} else {
			ig <- g[[k]]
		}
		#E(ig)$weight; plot(as_graphnel(ig))

		# SEM fitting
		fit <- SEMfit(graph = ig, data, group, B = 0.1, perm = 500)
		if(length(fit[[1]]) == 0) {
			res.tbl <- rbind(res.tbl, rep(NA, 7))
			next
		}
		p <- ncol(fit$dataXY)
		B <- inspect(fit[[1]], "est")$beta[-p, -p]
		if(sum(B) == 0) {
			res.tbl <- rbind(res.tbl, rep(NA, 7))
			next
		}
		pval<- fit$gest@res$pchisq[1:(p-1)]
		#DE <- c(DE, list())

		# data.frame of combined SEM results :
		df <- data.frame(No.nodes = vcount(ig),
		                 #No.edges= ecount(ig),
		                 No.DEGs = sum(p.adjust(pval, method = "BH") < alpha),
		                 pD = round(fit$gest@res[p, 5], 6),
		                 pA = round(fit$gest@res[p + 1, 5], 6),
		                 pE = round(fit$gest@res[p + 2, 5], 6))
		P.VALUE <- round(pX2(x = df[1, 3:5]), 6)
		res.tbl <- rbind(res.tbl, cbind(df, P.VALUE))
	}

	#colnames(res.tbl) <- c(...)
	rownames(res.tbl) <- gs
	return(na.omit(res.tbl))
}

quiet <- function(x) {
	sink(tempfile())
	on.exit(sink())
	invisible(force(x))
}

