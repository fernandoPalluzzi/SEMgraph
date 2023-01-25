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

#' @title Fit a graph as a Structural Equation Model (SEM)
#'
#' @description \code{SEMrun()} converts a (directed, undirected, or mixed)
#' graph to a SEM and fits it. If a binary group variable (i.e., case/control)
#' is present, node-level or edge-level perturbation is evaluated.
#' This function can handle loop-containing models, although multiple
#' links between the same two nodes (including self-loops and mutual
#' interactions) and bows (i.e., a directed and a bidirected link between
#' two nodes) are not allowed.
#'
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If \code{NULL} (default), group influence will
#' not be considered.
#' @param fit A numeric value indicating the SEM fitting mode.
#' If \code{fit = 0} (default), no group effect is considered.
#' If \code{fit = 1}, a "common" model is used to evaluate group effects
#' on graph nodes.
#' If \code{fit = 2}, a two-group model is used to evaluate group effects
#' on graph edges.
#' @param algo MLE method used for SEM fitting. If \code{algo = "lavaan"}
#' (default), the SEM will be fitted using the NLMINB solver from
#' \code{lavaan} R package, with standard errors derived from the expected
#' Fisher information matrix. If \code{algo = "ricf"}, the model is fitted
#' via residual iterative conditional fitting (RICF; Drton et al. 2009).
#' If \code{algo = "cggm"}, model fitting is based on constrained Gaussian
#' Graphical Modeling (CGGM; Hastie et al. 2009, p. 446).
#' @param start Starting value of SEM parameters for \code{algo = "lavaan"}.
#' If start is \code{NULL} (default), the algorithm will determine the
#' starting values. If start is a numeric value, it will be used as a
#' scaling factor for the edge weights in the graph object (graph attribute
#' \code{E(graph)$weight}).
#' For instance, a scaling factor is useful when weights have fixed values
#' (e.g., 1 for activated, -1 for repressed, and 0 for unchanged interaction).
#' Fixed values may compromise model fitting, and scaling them is a safe
#' option to avoid this problem. As a rule of thumb, to our experience,
#' \code{start = 0.1} generally performs well with {-1, 0, 1} weights.
#' @param SE If "standard" (default), with \code{algo = "lavaan"},
#' conventional standard errors are computed based on inverting the observed
#' information matrix. If "none", no standard errors are computed.
#' @param n_rep Number of randomization replicates (default = 1000),
#' for permutation flip or boostrap samples, if \code{algo = "ricf"}.
#' @param limit An integer value corresponding to the network size
#' (i.e., number of nodes). Beyond this limit, the execution under
#' \code{algo = "lavaan"} will run with \code{SE = "none"}, if 
#' \code{fit = 0}, or will be ridirected to \code{algo = "ricf"}, if
#' \code{fit = 1}, or to \code{algo = "cggm"}, if \code{fit = 2}.
#' This redirection is necessary to reduce the computational demand of
#' standard error estimation by lavaan. Increasing this number will
#' enforce lavaan execution when \code{algo = "lavaan"}.
#' @param ... Currently ignored.
#'
#' @details SEMrun maps data onto the input graph and converts it into a
#' SEM. Directed connections (X -> Y) are interpreted as direct causal
#' effects, while undirected, mutual, and bidirected connections are
#' converted into model covariances. SEMrun output contains different sets
#' of parameter estimates. Beta coefficients (i.e., direct effects) are
#' estimated from directed interactions and residual covariances (psi
#' coefficients) from bidirected, undirected, or mutual interactions.
#' If a group variable is given, exogenous group effects on nodes (gamma
#' coefficients) or edges (delta coefficients) will be estimated.
#' By default, maximum likelihood parameter estimates and P-values for
#' parameter sets are computed by conventional z-test (= estimate/SE),
#' and fits it through the \code{\link[lavaan]{lavaan}} function, via
#' Maximum Likelihood Estimation (estimator = "ML", default estimator in
#' \code{\link[lavaan]{lavOptions}}).
#' In case of high dimensionality (n.variables >> n.subjects), the covariance
#' matrix could not be semi-definite positive and thus parameter estimates
#' could not be done. If this happens, covariance matrix regularization
#' is enabled using the James-Stein-type shrinkage estimator implemented
#' in the function \code{\link[corpcor]{pcor.shrink}} of corpcor R package.
#' Argument \code{fit} determines how group influence is evaluated in the
#' model, as absent (\code{fit = 0}), node perturbation (\code{fit = 1}),
#' or edge perturbation (\code{fit = 2}). When \code{fit = 1}, the group
#' is modeled as an exogenous variable, influencing all the other graph
#' nodes. When \code{fit = 2}, SEMrun estimates the differences
#' of the beta and/or psi coefficients (network edges) between groups.
#' This is equivalent to fit a separate model for cases and controls,
#' as opposed to one common model perturbed by the exogenous group effect.
#' Once fitted, the two models are then compared to assess significant
#' edge (i.e., direct effect) differences (d = beta1 - beta0).
#' P-values for parameter sets are computed by z-test (= d/SE), through
#' \code{\link[lavaan]{lavaan}}. As an alternative to standard P-value
#' calculation, SEMrun may use either RICF (randomization P-values) or
#' GGM (de-sparsified P-values) methods. These algorithms are much faster
#' than \code{\link[lavaan]{lavaan}} in case of large input graphs.
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "fit", SEM fitted lavaan, ricf, or ggmncv object,
#' depending on the MLE method specified by the \code{algo} argument;
#' \item "gest" or "dest", a data.frame of node-specific
#' ("gest") or edge-specific ("dest") group effect estimates and P-values;
#' \item "model", SEM model as a string if \code{algo = "lavaan"},
#' and \code{NULL} otherwise;
#' \item "graph", the induced subgraph of the input network mapped
#' on data variables. Graph edges (i.e., direct effects) with P-value < 0.05
#' will be highlighted in red (beta > 0) or blue (beta < 0). If a group
#' vector is given, nodes with significant group effect (P-value < 0.05)
#' will be red-shaded (beta > 0) or lightblue-shaded (beta < 0);
#' \item "dataXY", input data subset mapping graph nodes, plus
#' group at the first column (if no group is specified, this column will
#' take NA values).
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom boot boot
#' @importFrom graphics abline curve hist legend par
#' @importFrom mgcv gam
#' @importFrom stats approx as.dist coefficients cor cov
#'             cutree dchisq dnorm formula hclust lm lm.fit
#' 			   na.omit p.adjust pchisq pnorm pt qchisq
#'             qnorm quantile rnorm runif sd var
#' @importFrom utils flush.console
#' @export
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Pearl J (1998). Graphs, Causality, and Structural Equation Models.
#' Sociological Methods & Research., 27(2):226-284.
#' <https://doi.org/10.1177/0049124198027002004>
#'
#' Yves Rosseel (2012). lavaan: An R Package for Structural Equation
#' Modeling. Journal of Statistical Software, 48(2): 1-36.
#' <https://www.jstatsoft.org/v48/i02/>
#'
#' Pepe D, Grassi M (2014). Investigating perturbed pathway modules
#' from gene expression data via Structural Equation Models. BMC
#' Bioinformatics, 15: 132.
#' <https://doi.org/10.1186/1471-2105-15-132>
#'
#' Drton M, Eichler M, Richardson TS (2009). Computing Maximum Likelihood
#' Estimated in Recursive Linear Models with Correlated Errors.
#' Journal of Machine Learning Research, 10(Oct): 2329-2348.
#' <https://www.jmlr.org/papers/volume10/drton09a/drton09a.pdf>
#'
#' Jankova J, van de Geer S (2015). Confidence intervals for high-dimensional
#' inverse covariance estimation. Electronic Journal of Statistics,
#' 9(1): 1205-1229. <https://doi.org/10.1214/15-EJS1031>
#'
#' Hastie T, Tibshirani R, Friedman J. (2009). The Elements of Statistical
#' Learning (2nd ed.). Springer Verlag: New York. ISBN: 978-0-387-84858-7
#'
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
#'
#' @seealso See \code{\link[ggm]{fitAncestralGraph}} and \code{\link[ggm]{fitConGraph}}
#' for RICF algorithm and constrained GGM algorithm details, respectively.
#'
#' @examples
#'
#' #### Model fitting (no group effect)
#'
#' sem0 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc))
#' summary(sem0$fit)
#' head(parameterEstimates(sem0$fit))
#'
#' # Graphs
#' gplot(sem0$graph, main = "significant edge weights")
#' plot(sem0$graph, layout = layout.circle, main = "significant edge weights")
#'
#'
#' #### Model fitting (common model, group effect on nodes)
#'
#' sem1 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc),
#'                group = sachs$group)
#'
#' # Fitting summaries
#' summary(sem1$fit)
#' print(sem1$gest)
#' head(parameterEstimates(sem1$fit))
#'
#' # Graphs
#' gplot(sem1$graph, main = "Between group node differences")
#' plot(sem1$graph, layout = layout.circle, main = "Between group node differences")
#'
#'
#' #### Two-group model fitting (group effect on edges)
#'
#' sem2 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc),
#'                group = sachs$group,
#'                fit = 2)
#'
#' # Summaries
#' summary(sem2$fit)
#' print(sem2$dest)
#' head(parameterEstimates(sem2$fit))
#'
#' # Graphs
#' gplot(sem2$graph, main = "Between group edge differences")
#' plot(sem2$graph, layout = layout.circle, main = "Between group edge differences")
#'
#' \donttest{
#'
#' # Fitting and visualization of a large pathway:
#'
#' g <- kegg.pathways[["MAPK signaling pathway"]]
#' G <- properties(g)[[1]]; summary(G)
#' 
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' g1 <- SEMrun(G, als.npn, alsData$group, algo = "cggm")$graph
#' g2 <- SEMrun(g1, als.npn, alsData$group, fit = 2, algo = "cggm")$graph
#'
#' # extract the subgraph with between group node and edge differences
#' g2 <- g2 - E(g2)[-which(E(g2)$color != "gray50")]
#' g <- properties(g2)[[1]]
#'
#' # plot graph
#' E(g)$color<- E(g2)$color[E(g2) %in% E(g)]
#' gplot(g, l = "fdp", main="node and edge group differences")
#'
#' }
#'
SEMrun <- function(graph, data, group = NULL, fit = 0, algo = "lavaan",
                   start = NULL, SE = "standard", n_rep = 1000,
				   limit = 100, ...)
{
	if (is.null(group) & fit != 0) fit <- 0
	if (!is.null(group) & fit == 0) fit <- 1

	if (fit == 0) {
		if (algo == "lavaan") {
			return(fit = SEMfit(graph = graph, data = data, group = NULL,
			                    start = start, SE = SE, limit = limit))
		} else if (algo == "cggm") {
			return(fit = SEMggm(graph = graph, data = data, group = NULL))
		} else if (algo == "ricf") {
			return(fit = SEMricf(graph = graph, data = data, group = NULL,
								 n_rep = n_rep))
		}
	}

	if (fit == 1) {
		if (algo == "lavaan") {
			return(fit = SEMfit(graph = graph, data = data, group = group,
			                    start = start, SE = SE, limit = limit))
		} else if (algo == "cggm") {
			return(fit = SEMggm(graph = graph, data = data, group = group))
		} else if( algo == "ricf" ) {
			return(fit = SEMricf(graph = graph, data = data, group = group,
			                     n_rep = n_rep))
		}
	}

	if (fit == 2) {
		if (algo == "lavaan") {
			return(fit = SEMfit2(graph = graph, data = data, group = group,
			                     start = start, SE = SE, limit = limit))
		} else if (algo == "cggm") {
			return(fit = SEMggm2(graph = graph, data = data, group = group))
		} else if (algo == "ricf") {
			return(fit = SEMricf2(graph = graph, data = data, group = group,
			                      n_rep = n_rep))
		}
	}
}

SEMmodel <- function(ig, nodes, group, ...)
{
	# Set from-to-matrix representation of gene-gene links
	ftm <- as_data_frame(ig)
	if (is.directed(ig) & sum(which_mutual(ig)) > 0) {
		dg <- ig - E(ig)[which_mutual(ig)]
		ug <- as.undirected(ig-E(ig)[!which_mutual(ig)])
		ftm <- as_data_frame(dg)
		ftb <- as_data_frame(ug)
	} else {
		ftb <- data.frame(NULL)
	}

	modelY <- modelV <- vector()
	if (is.directed(ig)) {
		if (nrow(ftm) > 0) {
			for(j in 1:nrow(ftm)) {
				modelY[j] <- paste0("z", ftm[j, 2], "~z", ftm[j, 1])
			}
		}
		if (nrow(ftb) > 0) {
			for(k in 1:nrow(ftb)) {
				modelV[k] <- paste0("z", ftb[k, 2], "~~z", ftb[k, 1])
			}
		}

	} else {
		for(j in 1:nrow(ftm)) {
			modelY[j] <- paste0("z", ftm[j, 2], "~~z", ftm[j, 1])
		}
	}

	# Group mean differences effect
	modelC <- sort(paste0("z", nodes, "~", "group"))

	# Equal residual variance
	#modelV <- sort(paste0("z", nodes, "~~v*z", nodes))

	# Unequal residual variance
	#modelV <- sort(paste0("z", nodes, "~~z", nodes))

	if (is.null(group)) {
		model <- paste0(c(sort(modelY), modelV))
	} else {
		model <- paste0(c(modelC, sort(modelY), modelV))
	}

	return(model)
}

SEMstart <- function(ig, data, group, a, ...)
{
	if (is.numeric(a)) {

		if (is.null(E(ig)$weight)) {
			ig <- weightGraph(ig, data, group = NULL, method = "r2z",
		                      seed = "none")
		    E(ig)$weight <- E(ig)$zsign
		}

		dg <- ig - E(ig)[which_mutual(ig)]
		ug <- as.undirected(ig - E(ig)[!which_mutual(ig)])
		ftm <- as_data_frame(dg)
		ftb <- as_data_frame(ug)

		if(nrow(ftm) > 0) {
			Reg <- paste0("z", ftm[, 2], "~start(", a*ftm[, 3], ")*z",
			              ftm[, 1])
		} else {
			Reg <- NULL
		}

		if (nrow(ftb) > 0) {
			Cov <- paste0("z", ftb[, 2], "~~start(", a*ftb$weight/2,
			              ")*z", ftb[,1])
		} else {
			Cov <- NULL
		}

		if (!is.null(group)) {
			ReG <- paste0("z", V(ig)$name, "~", "group")
		} else {
			ReG <- NULL
		}
		return(model = paste0(c(sort(ReG), sort(Reg), sort(Cov))))
	}

	if (!is.numeric(a)) {
		est <- quiet(SEMggm(ig, data, group)$fit$parameterEstimates)
		est <- rbind(est$Reg[, 1:4], est$Cov)

		if (is.null(group)) {
			Reg <- paste0("z", est[, 1], est[, 2], "start(", est[, 4],
			              ")*z", est[, 3])
		return(model = sort(Reg))

		} else {
			G <- which(est$rhs == "group")
			Reg <- paste0("z", est[-G, 1], est[-G, 2], "start(", est[-G, 4],
			              ")*z", est[-G, 3])
			ReG <- paste0("z", est[G, 1], est[G, 2], "start(", est[G, 4],
			              ")*", est[G, 3])
			return(model=paste0(c(sort(ReG), sort(Reg))))
		}
	}
}

SEMfit <- function(graph, data, group = NULL, start = NULL, fit = 0,
					SE = "standard", limit = 100,...)
{
	# Change SEM fitting if n.nodes > limit
	if (vcount(graph) > limit & is.null(group)) {
		message("WARNING: very large input graph (>", limit, " nodes) !
		 SEs are not computed...\n")
		SE <- "none"
	}
	if (vcount(graph) > limit & !is.null(group)) {
		message("WARNING: very large input graph (>", limit, " nodes) !
		 RICF solver activated...\n")
		return(fit = SEMricf(graph = graph, data = data, group = group))
	}

	# Set data and model objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	colnames(dataY) <- paste0("z", nodes)
	if (is.null(group)) {
		dataXY <- dataY
	} else {
		dataXY <- cbind(group, dataY)
	}
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	if (corpcor::is.positive.definite(cor(dataXY)[1:p, 1:p])) {
	 covXY <- cor(dataXY)[1:p, 1:p]
	} else {
	 covXY<- corpcor::cor.shrink(dataXY,verbose=TRUE)
	 if (attributes(covXY)$lambda > 0.5) message(
	  "WARNING: lambda > 0.5, correlation matrix is 'shrinked' near the identity matrix !\n")
	 covXY<- covXY[1:p,1:p]
	}

	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	if (is.null(start)) {
		model <- SEMmodel(ig, nodes, group)
	} else if (!is.numeric(start)) {
		model <- SEMstart(ig, data, group, a = start)
	} else if (is.numeric(start)) {
		model <- SEMstart(ig, data, group, a = start)
	}

	# SEM fitting based on lavaan
	suppressWarnings(
	fit <- lavaan(model, sample.cov = covXY, sample.nobs = n, se = SE,
	              fixed.x = TRUE, int.ov.free = TRUE, auto.var = TRUE,
	              information = "observed", observed.information = "hessian",
	              auto.cov.y = FALSE, control = list(abs.tol = 1e-20,
	              rel.tol = 1e-10))
    )
	if (fit@Fit@converged == TRUE) {
		srmr <- fitMeasures(fit, "srmr")
		dev <- fitMeasures(fit, "chisq")
		df <- fitMeasures(fit, "df")
		cat(paste0("NLMINB solver ended normally after ", fit@Fit@iterations,
		           " iterations"), "\n\n")
		cat("deviance/df:", dev/df, " srmr:", srmr, "\n\n")

	} else {
		cat("Model converged:", fit@Fit@converged, "\n\n")
		return(fit = NULL)
	}

	est <- parameterEstimates(fit)
	if (!is.null(group) & SE != "none") {
		gest <- est[1:(p - 1),]
		gest$lhs <- sub("z", "", gest$lhs)
		pval1 <- Brown.test(p = gest$pvalue, x = dataY, theta = gest$est,
		                    tail = "positive")
		pval2 <- Brown.test(p = gest$pvalue, x = dataY, theta = gest$est,
		                    tail = "negative")
		cat("Brown's combined P-value of node activation:", pval1, "\n\n")
		cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	} else {
		gest <- NULL
	}

	# Output objects
	if (SE != "none") ig <- colorGraph(est = est, graph = ig, group = group, alpha = 0.05)
	if (is.null(group)) dataXY <- cbind(group = rep(NA, n), dataXY)
	colnames(dataXY) <- gsub("z", "", colnames(dataXY))

	return(list(fit = fit, gest = gest, model = model, graph = ig,
	            dataXY = dataXY))
}

SEMfit2 <- function(graph, data, group, start = NULL, SE = "standard",
                    limit = 100, ...)
{
	# Model fitting with GGM algo if n.nodes > limit
	if (vcount(graph) > limit) {
		message("WARNING: input graph is very large ( >", limit, " nodes ) !
		 GGM (constrained) solver activated...\n")
		return(fit = SEMggm2(graph = graph, data = data, group = group))
	}

	# Set data and model objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	colnames(dataY) <- paste0("z", nodes)
	p <- ncol(dataY)

	# covariances for cases (GROUP 1)
	data1<- dataY[group==1,]
	n1<- nrow(data1)
	if( corpcor::is.positive.definite(cor(data1)[1:p,1:p]) ){
	 cov1<- cor(data1)[1:p,1:p]
	}else{
	 cov1<- corpcor::cor.shrink(data1,verbose=TRUE)
	 if (attributes(cov1)$lambda > 0.5) message(
	  "WARNING: lambda > 0.5, correlation matrix is 'shrinked' near the identity matrix !\n")
	 cov1<- cov1[1:p,1:p]
	}
	
	# covariances for controls (GROUP 0)
	data0<- dataY[group==0,]
	n0<- nrow(data0)
	if( corpcor::is.positive.definite(cor(data0)[1:p,1:p]) ){
	 cov0<- cor(data0)[1:p,1:p]
	}else{
	 cov0<- corpcor::cor.shrink(data0,verbose=TRUE)
	 if (attributes(cov0)$lambda > 0.5) message(
	  "WARNING: lambda > 0.5, correlation matrix is 'shrinked' near the identity matrix !\n")
	 cov0<- cov0[1:p,1:p]
	}

	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))

	if (is.null(start)) {
		model <- SEMmodel(ig, nodes, group = NULL)
	} else if (!is.numeric(start)) {
		model <- SEMstart(ig, data, group = NULL, a = start)
	} else if (is.numeric(start)) {
		model <- SEMstart(ig, data, group = NULL, a = start)
	}

	# SEM fitting with lavaan
	suppressWarnings(
	fit <- lavaan(model, sample.cov = list(cov0,cov1),
	              sample.nobs = list(n0, n1), se = SE, fixed.x = TRUE,
	              int.ov.free = TRUE, auto.var = TRUE, auto.cov.y = FALSE,
	              information = "observed", observed.information = "hessian",
	              control = list(abs.tol = 1e-20, rel.tol = 1e-10))
	)
	if (fit@Fit@converged == TRUE) {
		srmr <- fitMeasures(fit, "srmr")
		dev <- fitMeasures(fit, "chisq")
		df <- fitMeasures(fit, "df")
		cat(paste0("NLMINB solver ended normally after ", fit@Fit@iterations,
		           " iterations"),"\n\n")
		cat("deviance/df:", dev/df, " srmr:", srmr, "\n\n")
	} else {
		cat("Model converged:", fit@Fit@converged, "\n\n")
		return(fit = NULL)
	}

	est <- parameterEstimates(fit)
	est0 <- est[est$group == 1 & est$op == "~",]
	est1 <- est[est$group == 2 & est$op == "~",]
	d_est <- est1$est - est0$est

	if (sum(d_est != 0)) {
		d_se <- sqrt(est1$se^2 + est0$se^2)
		pvalue <- 2*(1-pnorm(abs(d_est/d_se)))
		d_lower <- d_est - 1.96*d_se
		d_upper <- d_est + 1.96*d_se
		pval1 <- Brown.test(p = pvalue, x = NULL, theta = d_est, tail = "positive")
		pval2 <- Brown.test(p = pvalue, x = NULL, theta = d_est, tail = "negative")
		cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
		cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
		# Output objects
		dest <- cbind(est0[, 1:3], d_est, d_se, d_z = d_est/d_se, pvalue, d_lower, d_upper)
		dest$lhs <- sub("z", "", dest$lhs)
		dest$rhs <- sub("z", "", dest$rhs)
		class(dest)<- c("lavaan.data.frame" ,"data.frame")
		if (SE != "none") ig <- colorGraph(est = dest, graph = ig, group = NULL, alpha = 0.05)
	} else {
		dest <- NULL
	}
	dataXY <- cbind(c(rep(1, n1), rep(0, n0)), rbind(data1, data0))
	colnames(dataXY) <- gsub("z", "", colnames(dataXY))

	return(list(fit = fit, dest = dest, model = model, graph = ig,
	            dataXY = dataXY))
}

SEMricf<- function (graph, data, group = NULL, random.x = FALSE, n_rep = 1000, ...) 
{
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	if (is.null(group)) {
	 dataXY <- dataY
	}else{
	 dataXY <- cbind(group, dataY)
	}
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	if (corpcor::is.positive.definite(cor(dataXY)[1:p, 1:p])) {
	 covXY <- cor(dataXY)[1:p, 1:p]
	}else{
	 covXY<- corpcor::cor.shrink(dataXY,verbose=TRUE)
	 if (attributes(covXY)$lambda > 0.5) message(
	  "WARNING: lambda > 0.5, correlation matrix is 'shrinked' near the identity matrix !\n")
	 covXY<- covXY[1:p,1:p]
	}
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	E(ig)$weight <- ifelse(which_mutual(ig), 100, 1)
	A <- as_adj(ig, attr = "weight", sparse = FALSE)[nodes,nodes]
	if (is.null(group) & random.x == TRUE) {
		Vx <- which(colSums(A) == 0)
		A[Vx, Vx] <- 100
		diag(A) <- 0
	}
	if (!is.null(group)) {
		A <- cbind(rep(0, p), rbind(rep(1, p - 1), A))
		colnames(A)[1] <- rownames(A)[1] <- "group"
	}
	fit <- ggm::fitAncestralGraph(amat = A, S = covXY, n = n, tol = 1e-06)
	cat(paste0("RICF solver ended normally after ", fit$it, 
		" iterations"), "\n\n")
	idx <- fitIndices(n, fit$df, covXY, fit$Shat)
	cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3], 7), "\n\n")
	
	est <- parameterEstimates.RICF(fit)
	gest <- list(NULL, NULL)
	pval <- NULL
	
	if (!is.null(group) & n_rep != 0) {
	 gest <- flip.RICF(fit = fit, data = dataY, group = group, n_rep = n_rep)
	 pval1 <- Brown.test(p = gest[[1]][,4], x = dataY,
		theta = gest[[1]][,2], tail = "positive")
	 pval2 <- Brown.test(p = gest[[1]][,4], x = dataY, 
		theta = gest[[1]][,2], tail = "negative")
	 cat("Brown's combined P-value of node activation:", pval1, "\n\n")
	 cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	 ig <- colorGraph(est = gest[[1]], ig, group, alpha = 0.05)
	 pval <- c(pval1, pval2)
	}
	if (is.null(group) & n_rep != 0) {
	 est <- boot.RICF(A, dataXY, group=NULL, est, n_rep)
	 ig <- colorGraph(est, ig, group=NULL, alpha=0.05)
	}
	if (is.null(group)) dataXY<- cbind(group=rep(NA, n), dataXY)
	
	fit <- list(ricf = fit, fitIdx = idx, parameterEstimates = est)
	class(fit) <- "RICF"
	
	return(list(fit = fit, gest = gest[[1]], model = NULL, graph = ig, 
		dataXY = dataXY, r_gest = gest[[2]], pval = pval))
}

SEMricf2 <- function(graph, data, group, random.x = FALSE, n_rep = 1000, ...)
{
	# Set graph and data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	dataY <- data[, nodes]
	data1 <- dataY[group == 1,]
	data0 <- dataY[group == 0,]
	n1 <- nrow(data1)
	n0 <- nrow(data0)
	n <- n1 + n0

	# Fitting RICF for group = 1
	fit1 <- quiet(SEMricf(graph, data1, group = NULL, random.x = FALSE, n_rep = 0))
	est1 <- fit1$fit$parameterEstimates

	# Fitting RICF for group = 0
	fit0 <- quiet(SEMricf(graph, data0, group = NULL, random.x = FALSE, n_rep = 0))
	est0 <- fit0$fit$parameterEstimates

	# Two-group fit indices
	it <- fit1$fit$ricf$it + fit0$fit$ricf$it
	cat(paste0("RICF solver ended normally after ", it, " iterations"), "\n\n")
	srmr <- (n1/n)*fit1$fit$fitIdx[3] + (n0/n)*fit0$fit$fitIdx[3]
	dev <- fit1$fit$fitIdx[1] + fit0$fit$fitIdx[1]
	df <- fit1$fit$fitIdx[2] + fit0$fit$fitIdx[2]
	cat("deviance/df:", dev/df , " srmr:", round(srmr, 7), "\n\n")

	# edge differences based on group randomization with n_rep=1000
	d_est<- est1$est - est0$est
	dest<- cbind(est1[,c(1:3)], d_est)[est1$op == "~",]
	
	if (n_rep != 0) {
	 A<- ifelse(abs(fit1$fit$ricf$Bhat) > 0, 1, 0)
	 diag(A)<- 0
	 dest<- boot.RICF(t(A), dataY, group, dest, R=n_rep)
	 pval1<- Brown.test(p=dest$pvalue, x=NULL, theta=dest$d_est, tail="positive")
	 pval2<- Brown.test(p=dest$pvalue, x=NULL, theta=dest$d_est, tail="negative")
	 cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
     cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
	 ig<- colorGraph(est=dest, graph=ig, group=NULL, alpha=0.05)
	}

	# output objects:
	fit<- list(Group_1 = fit1[[1]], Group_0 = fit0[[1]])	
	dataXY<- cbind(c(rep(1,n1),rep(0,n0)),rbind(data1,data0))
		
	return(list(fit = fit, dest = dest, model = NULL, graph = ig, dataXY = dataXY))
}

flip.RICF<- function(fit, data, group, n_rep, ...)
{
	# Permutation pvalues of group -> nodes:
	p <- ncol(data)
	B <- (diag(p+1)-fit$Bhat)[-1,-1] #B[1:10,1:10]
	Y <- data[,colnames(B)] # Y[1:10,1:10]
	Z <- as.matrix(Y)%*%t(diag(p)-B)*sd(group) # head(Z)
	perm <- list(B = n_rep + 1, seed = 123)
	gest <- flip::flip(Z, ~group, perms = perm, statTest = "t")
	
	# N(0,1) approximate pvalues:
	aveT <- apply(gest@permT[-1,], 2, mean)
	sdT <- apply(gest@permT[-1,], 2, sd)
	z <- abs(gest@permT[1,]- aveT)/sdT
	gest@res[,4]<- 2*(1-pnorm(z))
	colnames(gest@res)[4] <- "pvalue"
	rownames(gest@res)<- sub("X", "", rownames(gest@res))
	
	return( list(gest=gest@res, r_gest=gest@permT) )
}

boot.RICF <- function(A, Z, group, L, R, ...)
{
	est<- function(Z, i){
	 Z <- Z[i,]
	 if (corpcor::is.positive.definite(cor(Z))) {
	  covXY <- cor(Z)
	 }else{
	  covXY<- corpcor::cor.shrink(Z, verbose = FALSE)
	 }
	 fit <- ggm::fitAncestralGraph(A, covXY, n=nrow(Z), tol = 1e-06)
	 B <- parameterEstimates.RICF(fit)$est
	}
	dest<- function(Z, i){
	 Zi<- Z[i,]
	 D <- est(Zi[group == 1,]) - est(Zi[group == 0,])
	}

	message(paste0("Model randomization with B = ", R, " bootstrap samples ...\n\n"))
	ncpus<- parallel::detectCores(logical = FALSE)
	if (is.null(group)){
	 xboot<- boot::boot(Z, est, R, parallel="snow", ncpus=ncpus)
	 #xboot<- boot::boot(Z, FUN1, R=R)
	}else{
	 xboot<- boot::boot(Z, dest, strata=group, R=R, parallel="snow", ncpus=ncpus)
	 #xboot<- boot::boot(Z, est, strata=group, R=R)
	}
	t0<- xboot$t0[1:nrow(L)]
	se<- apply(xboot$t[,1:nrow(L)], MARGIN=2, sd)
	z <- t0/se
	
	est<- data.frame(L, se = se, z=z,
	 pvalue = 2*(1-pnorm(abs(z))),
	 ci.lower = (t0-1.96*se),
	 ci.upper = (t0+1.96*se))
	class(est)<- c("lavaan.data.frame" ,"data.frame")

	return(est)
}

parameterEstimates.RICF <- function(object, ...)
{
	# Convert into a vector the output of parameter estimates
	p <- nrow(object$Bhat)
	B <- gdata::unmatrix(diag(p) - object$Bhat, byrow = FALSE)
	B <- B[which(B != 0)]
	O <- ifelse(lower.tri(object$Ohat, diag = TRUE), object$Ohat, 0)
	diag(O) <- ifelse(diag(O) == 0, 1, diag(O))
	rownames(O) <- colnames(O) <- rownames(object$Ohat)
	O <- gdata::unmatrix(O, byrow = FALSE)
	O <- O[which(O != 0)]
	P <- c(B, O)

	est <- NULL
	for(j in 1:length(P)) {
		s <- strsplit(names(P)[j], ":")
		lhs <- s[[1]][1]
		rhs <- s[[1]][2]
		if (j <= length(B)) {
			est <- rbind(est, data.frame(lhs, op = "~", rhs, est = P[j],
			             stringsAsFactors = FALSE))
		} else {
			est <- rbind(est, data.frame(lhs, op = "~~", rhs, est = P[j],
			             stringsAsFactors = FALSE))
		}
	}
	rownames(est)<- NULL
	
	return(est)
}

SEMggm<- function(graph, data, group = NULL, method = "none", alpha = 0.05, ...) 
{
	# Set data objects :
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY<- data[,nodes]
	if (is.null(group)){
     dataXY <- dataY
    }else{ dataXY <- cbind(group, dataY)}
	n<- nrow(dataXY)
	p<- ncol(dataXY)
	if( corpcor::is.positive.definite(cor(dataXY)[1:p,1:p]) ){
	 covXY<- cor(dataXY)[1:p,1:p]
	}else{ 
	 covXY<- corpcor::cor.shrink(dataXY,verbose=TRUE)
	 if (attributes(covXY)$lambda > 0.5) message(
	  "WARNING: lambda > 0.5, correlation matrix is 'shrinked' near the identity matrix !\n")
	 covXY<- covXY[1:p,1:p]
	}

	# Set graph objects :
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	E(ig)$weight<- ifelse(which_mutual(ig), 100, 1)
	dadj<- as_adj(ig, attr="weight", sparse=FALSE)[nodes,nodes]
	diag(dadj) <- 100
	adj<- as_adj(as.undirected(ig), sparse=FALSE)[nodes,nodes]
	if( !is.null(group) ){
	 adj<- cbind(c(0,rep(1, p-1)), rbind(rep(1, p-1),adj))
	 colnames(adj)[1]<- rownames(adj)[1]<- "group"
	 dadj<- cbind(rep(0, p), rbind(rep(1, p-1),dadj))
	 colnames(dadj)[1]<- rownames(dadj)[1]<- "group"
	}
	
	# constrained GGM 
	cggm<- ggm::fitConGraph(amat=adj, S=covXY, n=n)
	Sigma<- cggm$Shat
	Theta<- solve(cggm$Shat)
	#cggm<- GGMncv::constrained(covXY, adj)
	#Sigma<- cggm$Sigma
	#Theta<- cggm$Theta
	rownames(Sigma)<- colnames(Sigma) <- colnames(dataXY)
	rownames(Theta)<- colnames(Theta) <- colnames(dataXY)	
	
	# Beta, Psi & Sigma matrices:
	betas<- function(Theta){
	 sapply(1:p, function(x) Theta[x,]/Theta[x,x]) * -1
	}
	Beta<- ifelse(dadj == 1, betas(Theta), 0)
	Psi<- ifelse(dadj == 100, Sigma, 0)
	diag(Psi)<- apply(scale(dataXY) %*% (diag(p)- Beta), 2, var)
	df<- p*(p+1)/2 - (sum(Beta != 0) + (sum(Psi != 0)-p)/2 + p)
	cat(paste0("GGM (constrained) solver ended normally after ", 0, " iterations"),"\n\n")
	idx<- fitIndices(n, df, covXY, Sigma, Theta)
	cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3],7), "\n\n")

	est <- parameterEstimates.GGM(dadj, Beta, Psi, Theta, R=covXY, n=n)
	if (!is.null(group)) {
	 gest<- est[est$op == "~" & est$rhs == "group",]
	 pval1<- Brown.test(p=gest$pvalue, x=dataY, theta=gest$est, tail="positive")
	 pval2<- Brown.test(p=gest$pvalue, x=dataY, theta=gest$est, tail="negative")
	 cat("Brown's combined P-value of node activation:", pval1, "\n\n")
     cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	}else{gest<- NULL}

	# output objects:
	Reg <- est[which(est$op == "~"),]
	ig<- colorGraph(est=Reg, graph=ig, group=group, alpha=0.05)
	if (is.null(group)) dataXY<- cbind(group=rep(NA, n), dataXY)
	fit<- list(Theta=Theta, Beta=Beta, Psi=Psi, fitIdx=idx, parameterEstimates=est)
	class(fit)<- "GGM"
	
	return( list(fit=fit, gest=gest, model=NULL, graph=ig, dataXY=dataXY) )	
}

SEMggm2<- function(graph, data, group, method = "none", alpha = 0.05, ...) 
{
	# Set graph and data objects:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	ig<- induced_subgraph(graph, vids= which(V(graph)$name %in% nodes))
	E(ig)$weight<- ifelse(which_mutual(ig), 100, 1)
	dadj<- as_adj(ig, attr="weight", sparse=FALSE)[nodes,nodes]
	dataY<- data[,nodes]
	data1<- dataY[group == 1,]
	data0<- dataY[group == 0,]
	n1<- nrow(data1)
	n0<- nrow(data0)
	n<- n1 + n0
	p<- ncol(dataY)
		
	# constrained GGM for group = 1
	cggm1<- quiet(SEMggm(graph, data1, group=NULL, method = method, alpha = alpha))
	est1<- parameterEstimates(cggm1$fit)
	# constrained GGM for group = 0
	cggm0<- quiet(SEMggm(graph, data0, group=NULL, method = method, alpha = alpha))
	est0<- parameterEstimates(cggm0$fit)
	# two-group fit indices
	cat(paste0("GGM (constrained) solver ended normally after ", 0, " iterations"),"\n\n")
	srmr<- (n1/n)*cggm1$fit$fitIdx[3] + (n0/n)*cggm0$fit$fitIdx[3]
	dev<- cggm1$fit$fitIdx[1] + cggm0$fit$fitIdx[1]
	df<- cggm1$fit$fitIdx[2] + cggm0$fit$fitIdx[2]
	cat("deviance/df:", dev/df , " srmr:", round(srmr,7), "\n\n")
		
	# Edge differences based on the de-sparsified precision matrix	
	est1<- est1[est1$op == "~",] #case
	est0<- est0[est0$op == "~",] #control
	d_omega<- est1$omega - est0$omega #case - control
	d_est <- est1$est - est0$est
	if (sum(d_est) != 0) {
	 d_se<- sqrt(est1$se^2 + est0$se^2)
	 d_z<- d_omega/d_se
	 pvalue<- 2*(1-pnorm(abs(d_z)))
	 d_lower<- d_omega - 1.96*d_se
	 d_upper<- d_omega + 1.96*d_se
	 pval1<- Brown.test(p=pvalue, x=NULL, theta=d_est, tail="positive")
	 pval2<- Brown.test(p=pvalue, x=NULL, theta=d_est, tail="negative")
	 cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
     cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
	 dest<- cbind(est0[,1:3], d_est, d_omega, d_se, d_z, pvalue, d_lower, d_upper)
	 class(dest)<- c("lavaan.data.frame" ,"data.frame")
	 ig<- colorGraph(est=dest, graph=ig, group=NULL, alpha=0.05)
	}else{ dest<- NULL }

	# output objects :
	fit<- list(Group_1 = cggm1$fit, Group_0 = cggm0$fit)
	dataXY<- cbind(c(rep(1,n1),rep(0,n0)),rbind(data1,data0))
		
	return( list(fit=fit, dest=dest, model=NULL, graph=ig, dataXY=dataXY) )	
}

parameterEstimates.GGM <- function(dadj, Beta, Psi, Theta, R, n, ...)
{
	# Convert into a vector the output of parameter estimates	
	A <- gdata::unmatrix(dadj, byrow=TRUE)
	B <- gdata::unmatrix(Beta, byrow=TRUE)
	B <- na.omit(ifelse(A == 1, B, NA))
	O <- ifelse(lower.tri(Psi, diag=TRUE), Psi, 0)
	rownames(O)<-colnames(O)<- rownames(Psi)
    O <- gdata::unmatrix(O, byrow=TRUE)
	O <- O[-which(O == 0)]
	
	L <- NULL
	P <- c(B,O)
	for(j in 1:length(P)) {
	 s<- strsplit(names(P)[j],":")
	 lhs<- s[[1]][2]
	 rhs<- s[[1]][1]
	 if(j <= length(B)){ 
	  L<- rbind(L, data.frame(lhs, op= "~", rhs))
	 }else{
	  L<- rbind(L, data.frame(lhs, op= "~~", rhs))}
	}
	
	# Equation 7 in
	# Jankova, J., & Van De Geer, S. (2015). Confidence intervals for high-dimensional
	# inverse covariance estimation. Electronic Journal of Statistics, 9(1), 1205-1229.
	Tstar <- 2 * Theta - Theta %*% R %*% Theta
	Tstar <- gdata::unmatrix(Tstar, byrow=TRUE)
	Bstar <- Tstar[names(Tstar) %in% names(B)]
	Ostar <- Tstar[names(Tstar) %in% names(O)]
	# standard error
	sds <- sqrt((tcrossprod(diag(Theta)) + Theta^2))
	sds <- gdata::unmatrix(sds, byrow=TRUE)
	Bse <- sds[names(sds) %in% names(B)]/sqrt(n)
	Ose <- sds[names(sds) %in% names(O)]/sqrt(n)
	# z test
	z <- c(Bstar,Ostar)/c(Bse,Ose)

	est <- data.frame(L, est = c(B,O),
	 omega = c(Bstar,Ostar), se = c(Bse,Ose), z = z,
	 pvalue = 2*(1-pnorm(abs(z))),
	 ci.lower = (c(Bstar,Ostar) - 1.96*c(Bse,Ose)),
	 ci.upper = (c(Bstar,Ostar) + 1.96*c(Bse,Ose)))
	rownames(est) <- NULL
	#class(est) <- c("lavaan.data.frame" ,"data.frame")

	return(est)
}

fitIndices <- function(n, df, S, Sigma, Theta = NULL, ...)
{
	p <- nrow(S)
	t <- p*(p + 1)/2 - df
	if (corpcor::is.positive.definite(Sigma) == FALSE){
	 message(" WARNING: Sigma is not a definite positive matrix...\n")
	}
	E <- S - Sigma

	# Deviance and df for model 1 (fitted model)
	if (is.null(Theta)) {
		ST <- S %*% solve(Sigma)
	} else {
		ST <- S %*% Theta
	}
	dev <- n*(sum(diag(ST)) - log(det(ST)) - p)

	# Deviance and df for model 0 (null model)
	dev0 <- n*(sum(diag(S)) - log(det(S)) - p) # n*(-log(det(R)))
	df0 <- p*(p + 1)/2 - p

	# Standardized Root Mean Square Residual (SRMR)
	SRMR <- sqrt(mean(E[lower.tri(E, diag = TRUE)]^2))

	# Root Mean Square Error of Approximation (RMSEA)
	RMSEA <- sqrt(max((dev - df), 0)/(df*(n - 1)))

	# Comparative Fit Index (CFI)
	CFI <- 1 - (dev - df)/(dev0 - df0)

	# Tucker-Lewis Index (TLI)
	TLI <- (dev0/df0 - dev/df)/(dev0/df0 - 1)

	# ULS Goodness of Fit Index (GFI)
	ULS <- 1 - sum(diag(t(E)%*%E))/sum(diag(t(S)%*%S))

	return(c(dev = dev, df = df, srmr = SRMR, rmsea = RMSEA, n = n, t = t))
}

Brown.test <- function (p, x = NULL, theta = NULL, tail = "both", ...)
{
	# from two-sided to one-sided (positive or negative) tests
	del <- which(is.na(p) | p <= 0 | p >= 1)
    if (length(del) > 0) {
     p <- p[-del]
     theta <- theta[-del]
    }
	if (tail == "positive") p <- ifelse(theta > 0, p/2, 1 - p/2)
	if (tail == "negative") p <- ifelse(theta > 0, 1 - p/2, p/2)
	
	#Fisher's (1932, 4th ed.) combined X2 test
	if (is.null(x)) {
	 pX2<- 1 - pchisq(q = -2 * sum(log(p)), df = 2 * length(p))
	 if (pX2 == 0) pX2 <- max(0, min(length(p)*min(p),1))
	 return(pX2)
	}
	
	#Brown's (1975) combined X2 test
	tmp <- c(-2.59, -2.382, -2.17, -1.946, -1.709, -1.458, -1.194,
			-0.916, -0.625, -0.32, 0, 0.334, 0.681, 1.044, 1.421,
			1.812, 2.219, 2.641, 3.079, 3.531, 4)
	s2X2 <- 4 * ncol(x) + 2 * sum(approx(seq(-1, 1, 0.1), tmp,
				xout = cor(x)[which(as.vector(lower.tri(cor(x))))])$x)
	EX2 <- 2 * ncol(x)
	
	# X2-correction, "c" = s2X2/(2 * 2 * k)
	# df-correction, "f" = 2 * (2 * k)^2 /s2X2
	fX2 <- -2 * sum(log(p))
	c <- s2X2/(2 * EX2)
	f <- 2 * EX2^2/s2X2
	if (f > 2 * length(p)) {
        f <- 2 * length(p)
        c <- 1.0
    }
	pX2 <- 1 - pchisq(q = fX2/c, df = f)
	if (pX2 == 0) pX2 <- max(0, min(length(p)*min(p),1))

	return(pX2)
}

#' @title Conditional Independence (CI) local tests of an acyclic graph
#'
#' @description P-values of one minimal testable implication (with the
#' smallest possible conditioning set) is returned per missing edge
#' given an acyclic graph (DAG or BAP) using the function
#' \code{\link[dagitty]{impliedConditionalIndependencies}} plus the
#' function \code{\link[dagitty]{localTests}} from package \code{dagitty}.
#' Without assuming any particular dependence structure, the p-values of
#' every CI test, in a DAG (BAP), is then combined using the Bonferroni’s
#' statistic in an overall test of the fitted model, B = K*min(p1,...,pK),
#' as reviewed in Vovk & Wang (2020).
#'
#' @param graph A directed graph as an igraph object.
#' @param data A data matrix with subjects as rows and variables as
#' columns.
#' @param bap If TRUE, the input graph is trasformend in a BAP, if FALSE
#' (defult) the input graph is reduced in a DAG.
#' @param verbose If TRUE, LocalCI results will be showed to
#' screen (default = TRUE).
#' @param limit An integer value corresponding to the number of missing
#' edges of the extracted acyclic graph. Beyond this limit, switch to
#' Shipley's C-test (Shipley 2000) is enabled to reduce the computational
#' burden.
#' By default, \code{limit = 10000}.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of three objects: (i) "dag": the DAG used to perform the localCI
#' test (ii) "msep": the list of all m-separation tests over missing edges in the
#' input graph and (iii) "mtest":the overall Bonferroni's P-value.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Vovk V, Wang R (2020). Combining p-values via averaging. Biometrika
#' 107(4): 791-808. <https://doi.org/10.1093/biomet/asaa027>
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs.
#' Structural Equation Modeling, 7(2): 206-218.
#' <https://doi.org/10.1207/S15328007SEM0702_4>
#'
#' @examples
#'
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' sem <- SEMrun(alsData$graph, als.npn)
#' B_test <- localCI.test(sem$graph, als.npn, verbose = TRUE)
#'
localCI.test<- function(graph, data, bap=FALSE, limit=10000, verbose=TRUE, ...)
{
	# Set graph and covariance matrix:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	dataY<- as.matrix(data[, nodes])
	df1<- vcount(graph)*(vcount(graph)-1)/2-ecount(as.undirected(graph))
	if (df1 > limit){
	 message("too many degree of fredoom, switch to Shipley.test...\n")
	 CI<- Shipley.test(graph, dataY)
	 return( list(bap=CI$dag, msep=CI$dsep, mtest=CI$ctest) )
	}

	# graph to DAG (or BAP) conversion :
	if (!is_dag(graph)){
	 cat("WARNING: input graph is not acyclic !\n")
	 if (bap){
	  cat(" Applying graph -> BAP conversion...\n")
	  bap<- graph2dag(graph, dataY, bap=bap) #del all <->
	 }else{
	  cat(" Applying graph -> DAG conversion...\n")
	  bap<- graph2dag(graph, dataY, bap=bap) #del cycles & all <->
	 }
	 df2<- vcount(bap)*(vcount(bap)-1)/2-ecount(as.undirected(bap))
	  cat(" \nDegrees of freedom:\n Input graph  =", 
             df1, "\n Output graph =", df2, "\n\n")
	}else{
	 bap <- graph
	 df2 <- df1
	}

	# d-separation local tests (minimal set) & Bonferroni's overall pvalue
	cat("d-separation test (minimal set) of", df2, "edges...\n")
	msep<- msep.test(bap=bap, S=cor(dataY), n=nrow(dataY), verbose=FALSE)
	#Bonferroni's multiple tests procedure: k*min(p1,...,pk) < alpha:
	mtest<- max(atanh(msep$estimate)^2*(nrow(dataY)-msep[,3]-3))[1]
	df<- 1
	pv<- min(p.adjust(msep$p.value, method="bonferroni"))[1]
	if (verbose) print(data.frame(B_test=mtest, df=df, B_pvalue=round(pv,6)))
	
	return( list(bap=bap, msep=msep, mtest=c(mtest, df, pv)) )
}

msep.test<- function(bap, S, n, verbose, ...)
{
	# conversion from igraph to dagitty
	if (is.dag(bap)) {
	 dagi<- graph2dagitty(bap, verbose=FALSE)
	}else{
	 dagi<- graph2dagitty(bap, verbose=FALSE)
	 dagi<- dagitty::canonicalize(dagi)$g
	}
	if (verbose) plot(dagitty::graphLayout(dagi))
	
	imp0<- dagitty::impliedConditionalIndependencies(dagi,
	 type="missing.edge", max.results=Inf)
	imp<- Filter(function(x) length(x$Z) <= 10, imp0)
	XY<- t(sapply(1:length(imp), function(x) imp[[x]][1:2]))
	K<- sapply(1:length(imp), function(x) length(imp[[x]][3]$Z))
	del<- which(duplicated(XY[,1:2]) == TRUE)
	if (length(del)>0) imp<- imp[-del]
	res<- dagitty::localTests(dagi, type = "cis", tests = imp,
	 sample.cov=S, sample.nobs=n, max.conditioning.variables=NULL, tol=0.05)
	if (length(del)>0) {
	 SET<- cbind(XY[-del,], "|Z|"=K[-del], res)
	}else{
	 SET<- cbind(XY, "|Z|"=K, res)
	}
	rownames(SET)<- NULL
	
	return( SET=na.omit(SET) )
}

#' @title Parameter Estimates of a fitted SEM
#'
#' @description Wrapper of the lavaan parameterEstimates() function
#' for RICF and CGGM algorithms
#'
#' @param fit A RICF or constrained GGM fitted model object.
#' @param ... Currently ignored.
#'
#' @return A data.frame containing the estimated parameters
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @examples
#' ricf1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "ricf")
#' parameterEstimates(ricf1$fit)
#'
#' cggm1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "cggm")
#' parameterEstimates(cggm1$fit)
#'
#' @export
parameterEstimates<- function(fit, ...)
{
	if (inherits(fit, "lavaan")){
	 est <- lavaan::parameterEstimates(fit)
	 est$lhs <- gsub("z", "", est$lhs)
	 est$rhs <- gsub("z", "", est$rhs)
	 if ("group" %in% colnames(est)){
	  est <- cbind(est[,-c(4,5)], group=est[,5])
	 }
	}
	else if (!inherits(fit, "lavaan")){
	 if (length(fit) != 2){
	  est0<- fit$parameterEstimates
	  sel<- which(est0$lhs == "group")
	  if(length(sel) != 1){
	   est <- est0
	  }else{
	   est <- est0[-sel,]
	   est <- rbind(est, est0[sel,])
	  }
	 }else{
	  est0 <- fit$Group_0$parameterEstimates
	  est1 <- fit$Group_1$parameterEstimates
	  group <- c(rep(1, nrow(est0)), rep(2,nrow(est1)))
	  est <- cbind(rbind(est0,est1),group)
	 } 
	}
	class(est)<- c("lavaan.data.frame" ,"data.frame")
	return(est)
}

#' @title RICF model summary
#'
#' @description Generate a summary for a RICF solver similar to
#' lavaan-formatted summary
#'
#' @param object A RICF fitted model object.
#' @param ... Currently ignored.
#'
#' @method summary RICF
#' @return Shown the lavaan-formatted summary to console
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @seealso \code{\link[SEMgraph]{SEMrun}}.
#'
#' @examples
#' sem1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "ricf")
#' summary(sem1$fit)
#'
#' @export
summary.RICF <- function(object, ...)
{
	.local <- function(object) {
		it <- object$ricf$it
		t <- object$fitIdx[6]
		n <- object$fitIdx[5]
		dev <- round(object$fitIdx[1], 3)
		df <- object$fitIdx[2]
		srmr <- round(object$fitIdx[3], 3)

		cat(paste0("RICF solver ended normally after ", it, " iterations"),
		           "\n\n")
		cat(paste0("  Estimator                                       ML"),
		           "\n")
		cat(paste0("  Optimization method                             RICF"),
		           "\n\n")
		cat(paste0("  Number of free parameters                       ",
		           t), "\n\n")
		cat(paste0("  Number of observations                          ",
		           n), "\n")
		cat("\nModel Test User Model\n\n")
		cat(paste0("  Test statistic (Deviance)                       ",
		           dev), "\n")
		cat(paste0("  Degrees of freedom (df)                         ",
		           df), "\n")
		cat(paste0("  Deviance/df                                     ",
		           round(dev/df, 3)), "\n")
		cat(paste0("  Standardized Root Mean Square Residual (srmr)   ",
		           srmr), "\n")
		cat("\nParameter Estimates:\n\n")

		#print(fit$parameterEstimates)
		K <- c("Regressions:", "Covariances:", "Variances:")
		L <- object$parameterEstimates
		reg <- L[which(L$op == "~"),]
		cov <- L[which(L$op == "~~"),]
		sel <- which(cov$lhs == cov$rhs)
		L <- list(reg, cov[-sel,], cov[sel,])
		for (l in 1:3) {
			cat(K[l], "\n\n")
			print(data.frame(lapply(L[[l]], function(y) {
					if(is.numeric(y)) round(y, 3) else y
				}
			)))
			cat("\n")
		}
	}
	.local(object)
}

#' @title GGM model summary
#'
#' @description Generate a summary for a constrained Gaussian Graphical
#' Model (GGM) similar to lavaan-formated summary
#'
#' @param object A constrained GGM fitted model object.
#' @param ... Currently ignored.
#'
#' @method summary GGM
#' @return Shown the lavaan-formatted summary to console
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @seealso \code{\link[SEMgraph]{SEMrun}}.
#'
#' @examples
#' sem0 <- SEMrun(sachs$graph, log(sachs$pkc), algo = "cggm")
#' summary(sem0$fit)
#'
#' @export
summary.GGM <- function(object, ...)
{
	.local <- function(object) {
		it <- 0
		t <- object$fitIdx[6]
		n <- object$fitIdx[5]
		dev <- round(object$fitIdx[1], 3)
		df <- object$fitIdx[2]
		srmr <- round(object$fitIdx[3], 3)

		cat(paste0("GGM (constrained) solver ended normally after ", it,
		           " iterations"), "\n\n")
		cat(paste0("  Estimator                                       ML"),
		           "\n")
		cat(paste0("  Optimization method                             CGGM"),
		           "\n\n")
		cat(paste0("  Number of free parameters                       ",
		           t), "\n\n")
		cat(paste0("  Number of observations                          ",
		           n), "\n")
		cat("\nModel Test User Model\n\n")
		cat(paste0("  Test statistic (Deviance)                       ",
		           dev), "\n")
		cat(paste0("  Degrees of freedom (df)                         ",
		           df), "\n")
		cat(paste0("  Deviance/df                                     ",
		           round(dev/df, 3)), "\n")
		cat(paste0("  Standardized Root Mean Square Residual (srmr)   ",
		           srmr), "\n")
		cat("\nParameter Estimates:\n\n")

		#print(object$parameterEstimates)
		K <- c("Regressions:", "Covariances:", "Variances:")
		L <- object$parameterEstimates
		reg <- L[which(L$op == "~"),]
		cov <- L[which(L$op == "~~"),]
		sel <- which(cov$lhs == cov$rhs)
		L <- list(reg, cov[-sel,], cov[sel,])
		for (l in 1:3) {
			cat(K[l], "\n\n")
			print(data.frame(lapply(L[[l]], function(y) {
					if(is.numeric(y)) round(y, 3) else y
				}
			)))
			cat("\n")
		}
	}
	.local(object)
}
