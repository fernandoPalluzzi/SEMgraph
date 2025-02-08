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
#' via residual iterative conditional fitting (RICF; Drton et al. 2009),
#' with standard error derived from randomization or bootstrap procedures.
#' If \code{algo = "cggm"}, model fitting is based on constrained Gaussian
#' Graphical Modeling (CGGM), with DAG nodewise Lasso procedure and
#' de-biasing asymptotic inference (Jankova & Van De Geer, 2019).
#' @param start Starting value of SEM parameters for \code{algo = "lavaan"}.
#' If start is \code{NULL} (default), the algorithm will determine the
#' starting values. If start is a numeric value, it will be used as a
#' scaling factor for the edge weights in the graph object (graph attribute
#' \code{E(graph)$weight}).
#' For instance, a scaling factor is useful when weights have fixed values
#' (e.g., 1 for activated, -1 for repressed, and 0 for unchanged interaction).
#' Fixed values may compromise model fitting, and scaling them is a safe
#' option to avoid this problem. As a rule of thumb, to our experience,
#' \code{start = 0.1} generally performs well with (-1, 0, 1) weights.
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
#' calculation, SEMrun may use either RICF (randomization or bootstrap
#' P-values) or GGM (de-biased asymptotically normal P-values) methods.
#' These algorithms are much faster than \code{\link[lavaan]{lavaan}}
#' in case of large input graphs.
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "fit", SEM fitted lavaan, ricf, or cggm object,
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
#' \item "data", input data subset mapping graph nodes, plus
#' group at the first column (if no group is specified, this column will
#' take NA values).
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom boot boot
#' @importFrom graphics abline curve hist legend par polygon
#' @importFrom mgcv gam
#' @importFrom stats approx as.dist coefficients cor cov cutree 
#'             dchisq dnorm density formula hclust lm lm.fit
#' 			   median na.omit p.adjust pchisq pnorm pt qchisq
#'             qnorm quantile rnorm runif sd var
#' @importFrom utils flush.console tail
#'
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
#' Jankova, J., & Van De Geer, S (2019). Inference in high-dimensional
#' graphical models. In Handbook of Graphical Models (2019).
#' Chapter 14 (sec. 14.2): 325-349. Chapman & Hall/CRC. ISBN: 9780429463976
#'
#' Hastie T, Tibshirani R, Friedman J. (2009). The Elements of Statistical
#' Learning (2nd ed.). Springer Verlag. ISBN: 978-0-387-84858-7
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
#' # Fitting and visualization of a large pathway:
#'
#' g <- kegg.pathways[["Neurotrophin signaling pathway"]]
#' G <- properties(g)[[1]]
#' summary(G)
#' 
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
#'
#' g1 <- SEMrun(G, als.npn, alsData$group, algo = "cggm")$graph
#' g2 <- SEMrun(g1, als.npn, alsData$group, fit = 2, algo = "cggm")$graph
#'
#' # extract the subgraph with node and edge differences
#' g2 <- g2 - E(g2)[-which(E(g2)$color != "gray50")]
#' g <- properties(g2)[[1]]
#'
#' # plot graph
#' E(g)$color<- E(g2)$color[E(g2) %in% E(g)]
#' gplot(g, l="fdp", psize=40, main="node and edge group differences")
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
			fit <- SEMfit(graph = graph, data = data, group = NULL,
							start = start, SE = SE, limit = limit)
		} else if (algo == "cggm") {
			fit <- SEMggm(graph = graph, data = data, group = NULL)
		} else if (algo == "ricf") {
			fit <- SEMricf(graph = graph, data = data, group = NULL,
							n_rep = n_rep)
		}
		class(fit) <- "SEM"
		return(fit)
	}

	if (fit == 1) {
		if (algo == "lavaan") {
			fit <- SEMfit(graph = graph, data = data, group = group,
							start = start, SE = SE, limit = limit)
		} else if (algo == "cggm") {
			fit <- SEMggm(graph = graph, data = data, group = group)
		} else if( algo == "ricf" ) {
			fit <- SEMricf(graph = graph, data = data, group = group,
							n_rep = n_rep)
		}
		class(fit) <- "SEM"
		return(fit)
	}

	if (fit == 2) {
		if (algo == "lavaan") {
			fit <- SEMfit2(graph = graph, data = data, group = group,
							start = start, SE = SE, limit = limit)
		} else if (algo == "cggm") {
			fit <- SEMggm2(graph = graph, data = data, group = group)
		} else if (algo == "ricf") {
			fit <- SEMricf2(graph = graph, data = data, group = group,
							n_rep = n_rep)
		}
		class(fit) <- "SEM"
		return(fit)
	}
}

SEMmodel <- function(ig, nodes, group, ...)
{
	# Set from-to-matrix representation of gene-gene links
	ftm <- igraph::as_data_frame(ig)
	if (is.directed(ig) & sum(which_mutual(ig)) > 0) {
		dg <- ig - E(ig)[which_mutual(ig)]
		ug <- as.undirected(ig-E(ig)[!which_mutual(ig)])
		ftm <- igraph::as_data_frame(dg)
		ftb <- igraph::as_data_frame(ug)
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
		ftm <- igraph::as_data_frame(dg)
		ftb <- igraph::as_data_frame(ug)

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
		cat(paste0("NLMINB solver ended normally after ", fit@Fit@iterations,
		           " iterations"), "\n\n")
		#srmr <- fitMeasures(fit, "srmr")
		#dev <- fitMeasures(fit, "chisq")
		#df <- fitMeasures(fit, "df")
		#cat("deviance/df:", dev/df, " srmr:", srmr, "\n\n")
		df <- fitMeasures(fit, "df")
		npar <- fitMeasures(fit, "npar")
		Shat <- fitted(fit)$cov
		Sobs <- cor(dataXY)[colnames(Shat), colnames(Shat)]
		idx <- fitIndices(n, df, npar, Sobs, Shat)
		#if (idx[7] > 0) {
		# message(paste0(" WARNING: deviance is estimated by removing ", idx[7], " singular values < 1e-12...\n"))
		#}
		cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3], 7), "\n\n")

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
	colnames(dataXY) <- sub(".", "", colnames(dataXY))

	return(list(fit = fit, gest = gest, model = model, graph = ig,
	            data = dataXY))
}

SEMfit2 <- function(graph, data, group, start = NULL, SE = "standard",
                    limit = 100, ...)
{
	# Model fitting with GGM algo if n.nodes > limit
	if (vcount(graph) > limit) {
		message("WARNING: input graph is very large ( >", limit, " nodes ) !
		 GGM (de-biased nodewise L1) solver activated...\n")
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
		cat(paste0("NLMINB solver ended normally after ", fit@Fit@iterations,
		           " iterations"),"\n\n")
		#srmr <- fitMeasures(fit, "srmr")
		#dev <- fitMeasures(fit, "chisq")
		#df <- fitMeasures(fit, "df")
		#cat("deviance/df:", dev/df, " srmr:", srmr, "\n\n")
		df <- fitMeasures(fit, "df")/2
		npar <-  fitMeasures(fit, "npar")/2
		Shat0 <- fitted(fit)$"Group 1"$cov
		Sobs0 <- cor(data0)[colnames(Shat0), colnames(Shat0)]
		Idx0 <- fitIndices(n0, df, npar, Sobs0, Shat0)
		Shat1 <- fitted(fit)$"Group 2"$cov
		Sobs1 <- cor(data1)[colnames(Shat1), colnames(Shat1)]
		Idx1 <- fitIndices(n1, df, npar, Sobs1, Shat1)
		srmr <- (n0/(n0+n1))*Idx0[3] + (n1/(n0+n1))*Idx1[3]
		dev <- Idx0[1] + Idx1[1]
		r <- Idx0[7] + Idx1[7]
		#if (r > 0) {
		# message(paste0(" WARNING: deviance is estimated by removing ", r, " singular values < 1e-12...\n"))
		#}
		cat("deviance/df:", dev/(2*df), " srmr:", round(srmr, 7), "\n\n")	

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
	colnames(dataXY) <- c("group", sub(".", "", colnames(dataXY)[-1]))

	return(list(fit = fit, dest = dest, model = model, graph = ig,
	            data = dataXY))
}

SEMricf<- function (graph, data, group = NULL, n_rep = 1000, ...) 
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
	A <- as_adjacency_matrix(ig, attr = "weight", sparse = FALSE)[nodes,nodes]
	if (!is.null(group)) {
		A <- cbind(rep(0, p), rbind(rep(1, p - 1), A))
		colnames(A)[1] <- rownames(A)[1] <- "group"
	}
	Vx<- which(colSums(A) == 0)
	Vy <- which(colSums(A) != 0)
	px <- length(Vx)
	py <- length(Vy)
	
	fit <- ggm::fitAncestralGraph(amat = A, S = covXY, n = n, tol = 1e-06)
	cat(paste0("RICF solver ended normally after ", fit$it, 
		" iterations"), "\n\n")
	B <- diag(p)- t(fit$Bhat[c(Vx,Vy), c(Vx,Vy)])
	O <- fit$Ohat[c(Vx,Vy), c(Vx,Vy)]
	if (is.null(group) & px > 1) {
	 O[1:px, 1:px]<- cor(dataXY[, Vx])
	}else{ O[1,1] <- 1 }
	npar<- sum(A == 1) + sum(A == 100)/2 + py
	df <- p*(p+1)/2 - npar - (px*(px+1)/2)
	Sobs <- cor(dataXY[, c(Vx,Vy)])
	#Shat <- solve((diag(p)-B) %*% solve(O) %*% t(diag(p)-B))
	Shat <- t(solve(diag(p)-B)) %*% O %*% solve(diag(p)-B)
	#Sobs <- covXY[c(Vx,Vy), c(Vx,Vy)]
	#Shat <- fit$Shat[c(Vx,Vy), c(Vx,Vy)]
	idx <- fitIndices(n, df, npar, Sobs, Shat)
	#if (idx[7] > 0) {
	# cat(paste0(" WARNING: deviance is estimated by removing ", idx[7], " singular values < 1e-12...\n\n"))
	#}
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

	fit <- list(Sigma=Shat, Beta=B, Psi=O, fitIdx=idx, parameterEstimates=est,
				it = fit$it, r_gest = gest[[2]], pval = pval)
	class(fit) <- "RICF"

	return(list(fit = fit, gest = gest[[1]], model = NULL, graph = ig, data = dataXY))
}

SEMricf2 <- function(graph, data, group, n_rep = 1000, ...)
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
	fit1 <- quiet(SEMricf(graph, data1, group = NULL, n_rep = 0))
	est1 <- fit1$fit$parameterEstimates
	# Fitting RICF for group = 0
	fit0 <- quiet(SEMricf(graph, data0, group = NULL, n_rep = 0))
	est0 <- fit0$fit$parameterEstimates

	# Two-group fit indices
	it <- fit1$fit$it + fit0$fit$it
	cat(paste0("RICF solver ended normally after ", it, " iterations"), "\n\n")
	srmr <- (n1/n)*fit1$fit$fitIdx[3] + (n0/n)*fit0$fit$fitIdx[3]
	dev <- fit1$fit$fitIdx[1] + fit0$fit$fitIdx[1]
	df <- fit1$fit$fitIdx[2] + fit0$fit$fitIdx[2]
	r <- fit1$fit$fitIdx[7] + fit0$fit$fitIdx[7]
	#if (r > 0) {
	# cat(paste0(" WARNING: deviance is estimated by removing ", r, " singular values < 1e-12...\n\n"))
	#}
	cat("deviance/df:", dev/df , " srmr:", round(srmr, 7), "\n\n")

	# edge differences based on group randomization with n_rep=1000
	d_est<- est1$est - est0$est
	dest<- cbind(est1[,c(1:3)], d_est)[est1$op == "~",]
	
	if (n_rep != 0) {
	 A<- ifelse(abs(fit1$fit$Beta) > 0, 1, 0)
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
	colnames(dataXY)[1]<- "group"
		
	return(list(fit = fit, dest = dest, model = NULL, graph = ig, data = dataXY))
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

	message(paste0("Model randomization with B = ", R, " bootstrap samples ...\n"))
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

parameterEstimates.RICF <- function(fit, ...)
{
	# Convert into a vector the output of parameter estimates
	p <- nrow(fit$Bhat)
	B <- gdata::unmatrix(diag(p) - fit$Bhat, byrow = FALSE)
	B <- B[which(B != 0)]
	O <- ifelse(lower.tri(fit$Ohat, diag = TRUE), fit$Ohat, 0)
	diag(O) <- ifelse(diag(O) == 0, 1, diag(O))
	rownames(O) <- colnames(O) <- rownames(fit$Ohat)
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
	class(est) <- c("lavaan.data.frame","data.frame")
	
	return(est)
}

SEMggm <- function(graph, data, group = NULL, ...) 
{
	# Set data and graph objects :
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	dag<- graph2dag(graph, data, bap=FALSE) #del cycles & all <->
	dataY<- data[, V(dag)$name]
	if (is.null(group)){
     dataXY<- dataY
    }else{
	 dataXY<- cbind(group, dataY)
	}
	n<- nrow(dataXY)
	p<- ncol(dataXY)

	din<- igraph::degree(dag, mode= "in")
	Vx<- V(dag)$name[din == 0]
	Vy<- V(dag)$name[din != 0]
	px<- length(Vx)
	py<- length(Vy)
	dadj<- as_adjacency_matrix(dag, sparse=FALSE)
		
	if( !is.null(group) ){
	 dadj<- cbind(rep(0, p), rbind(rep(1, p-1),dadj))
	 colnames(dadj)[1]<- rownames(dadj)[1]<- "group"
	 Vx<- "group"
	 Vy<- V(dag)$name
	 px<- 1
	 py<- p - 1
	}

	# de-biased (de-sparsified) DAG
	Ay<- as.matrix(dadj[,Vy])
	colnames(Ay)<- Vy
	res<- parameterEstimates.DAG(Ay, Z=scale(dataXY))
	est<- res$beta
	sigma<- res$sigma

	# Beta, Psi & Sigma matrices :
	gB<- graph_from_data_frame(est[,c(3,1,4)])
	V<- V(dag)$name[V(dag)$name %in% V(gB)$name == FALSE]
	gB<- gB + vertices(V)
	B<- as_adjacency_matrix(gB, attr="est", sparse=FALSE)
	B<- B[c(Vx,Vy), c(Vx,Vy)]
	O<- diag(sigma[c(Vx,Vy)])
	rownames(O)<- colnames(O)<- names(sigma)
	if(is.null(group) & px >1) O[1:px, 1:px]<- cor(dataXY[, Vx])
	Sobs <- cor(dataXY[, c(Vx,Vy)])
	npar<- sum(dadj != 0) + py
	df <- p*(p+1)/2 - npar - (px*(px+1)/2)
	Shat <- solve((diag(p)-B) %*% solve(O) %*% t(diag(p)-B))
	#Shat <- t(solve(diag(p)-B))%*%O%*%solve(diag(p)-B)
	
	cat(paste0("GGM (de-biased nodewise L1) solver ended normally after ", py, " iterations"),"\n\n")
	idx<- fitIndices(n, df, npar, Sobs, Shat)
	#if (idx[7] > 0) {
	# cat(paste0(" WARNING: deviance is estimated by removing ", idx[7], " singular values < 1e-12...\n\n"))
	#}
	cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3],7), "\n\n")

	if (!is.null(group)) {
	 Reg<- est[which(est$op == "~"),]
	 sel<- which(Reg$rhs == "group")
	 Reg<- rbind(Reg[sel,], Reg[-sel,])
	 gest<- Reg[sel,]
	 pval1<- Brown.test(x=dataY, p=gest$pvalue, theta=gest$est, tail="positive")
	 pval2<- Brown.test(x=dataY, p=gest$pvalue, theta=gest$est, tail="negative")
	 cat("Brown's combined P-value of node activation:", pval1, "\n\n")
     cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	}else{gest<- NULL}

	# output objects
	Reg<- est[which(est$op == "~"),]
	fit<- list(Sigma=Shat, Beta=B, Psi=O, fitIdx=idx, parameterEstimates=Reg)
	dag<- colorGraph(est=Reg, graph=dag, group=group, alpha=0.05)
	if (is.null(group)) dataXY<- cbind(group=rep(NA, n), dataXY)
	class(fit)<- "GGM"

	return( list(fit=fit, gest=gest, model=NULL, graph=dag, data=dataXY) )	
}

SEMggm2 <- function(graph, data, group, ...) 
{
	# Set data and graph objects :
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	dag<- graph2dag(graph, data, bap=FALSE) #del cycles & all <->
	dataY<- data[, V(dag)$name]
	data1<- dataY[group == 1,]
	data0<- dataY[group == 0,]
	n1<- nrow(data1)
	n0<- nrow(data0)
	n<- n1 + n0
	p<- ncol(dataY)

	# de-biased (de-sparsified) DAG for group = 1
	dag1<- quiet(SEMggm(dag, data1, group=NULL))
	est1<- dag1$fit$parameterEstimates
	# de-biased (de-sparsified) DAG for group = 0
	dag0<- quiet(SEMggm(dag, data0, group=NULL))
	est0<- dag0$fit$parameterEstimates

	# two-group fit indices
	cat(paste0("GGM (de-biased nodewise L1) solver ended normally after ", 0, " iterations"),"\n\n")
	srmr<- (n1/n)*dag1$fit$fitIdx[3] + (n0/n)*dag0$fit$fitIdx[3]
	dev<- dag1$fit$fitIdx[1] + dag0$fit$fitIdx[1]
	df<- dag1$fit$fitIdx[2] + dag0$fit$fitIdx[2]
	r<- dag1$fit$fitIdx[7] + dag0$fit$fitIdx[7]
	#if (r > 0) {
	# cat(paste0(" WARNING: deviance is estimated by removing ", r, " singular values < 1e-12...\n\n"))
	#}
	cat("deviance/df:", dev/df , " srmr:", round(srmr,7), "\n\n")

	# edge differences based on the de-sparsified Beta matrix
	d_est<- est1$est - est0$est
	if (sum(d_est) != 0) {
	 d_se<- sqrt(est1$se^2 + est0$se^2)
	 pvalue<- 2*(1-pnorm(abs(d_est/d_se)))
	 d_lower<- d_est - 1.96*d_se
	 d_upper<- d_est + 1.96*d_se
	 pval1<- Brown.test(p=pvalue, x=NULL, theta=d_est, tail="positive")
	 pval2<- Brown.test(p=pvalue, x=NULL, theta=d_est, tail="negative")
	 cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
     cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
	 dest<- cbind(est0[,1:3], d_est, d_se, d_z=d_est/d_se, pvalue, d_lower, d_upper)
	 class(dest)<- c("lavaan.data.frame" ,"data.frame")
	 dag<- colorGraph(est=dest, graph=dag, group=NULL, alpha=0.05) #gplot(dag)
	}else{ dest<- NULL }

	# output objects :
	fit<- list(Group_1 = dag1$fit, Group_0 = dag0$fit)
	dataXY<- cbind(c(rep(1,n1),rep(0,n0)),rbind(data1,data0)) 
	colnames(dataXY)[1]<- "group"

	return( list(fit=fit, dest=dest, model=NULL, graph=dag, data=dataXY) )	
}

parameterEstimates.DAG <- function(Ay, Z, ...)
{
	# Asymptotic inference in sparse DAGs with the procedure of:
	# Jankova, J., & Van De Geer, S. Inference in high-dimensional
	# graphical models. In Handbook of Graphical Models (2019).
	# Chapter 14: 325-349. Chapman & Hall/CRC. ISBN 9780429463976

	est <- NULL
	n <- nrow(Z)
	sigma <- rep(1, nrow(Ay)-ncol(Ay))
	
	for (j in 1:ncol(Ay)) { #j=1
	 yy <- colnames(Ay)[j]
	 xx <- rownames(Ay)[Ay[,j] == 1]
	 y <- as.matrix(Z[,yy])
	 x <- as.matrix(Z[,xx])
	 if (ncol(y) == 1) colnames(y) <- yy
	 if (ncol(x) == 1) colnames(x) <- xx
	 if (ncol(x) > 1) {
	    p <- ncol(x) + 1
		l <- ifelse(p > n, sqrt(log(p)/n), 1e-3)
		fit <- glmnet::glmnet(x=x, y=y, family="gaussian", lambda=l)
		b <- coef(fit, s = NULL)[-1,]
		wi <- glasso::glasso(s=cov(x), rho=l)$wi/n
		rownames(wi) <- colnames(wi) <- colnames(x)
	 }else{
		fit <- lm(y ~ x)
		b <- coef(fit)[2]
		names(b) <- colnames(x)
		wi <- 1/(n*var(x))
	 }
	 e <- y - x%*%b
	 B <- b + wi%*%t(x)%*%e/n
	 colnames(B) <- NULL
	 SE <- sd(e)*sqrt(diag(wi))
	 estj <- data.frame(
		lhs = rep(colnames(y), length(B)),
		op = "~",
		rhs = names(b),
		est = B,
		se = SE,
		z_test = B/SE,
		pvalue = 2*(1-pnorm(abs(B/SE))),
		ci.lower = B - 1.96*SE,
		ci.uppper = B + 1.96*SE
	  )
	 est <- rbind(est, estj)
	 sigma <- c(sigma, var(e))
	}
	rownames(est) <- NULL
	names(sigma) <- rownames(Ay)
	class(est) <- c("lavaan.data.frame","data.frame")

	return(list(beta = est, sigma = sigma))
}

fitIndices <- function(n, df, npar, S, Sigma, ...)
{
	p <- ncol(S)
	df <- as.numeric(df)
	t_free <- npar
	t_fixed <- p*(p+1)/2 - npar - df

	# remove latent variables (LV)
	colnames(S) <- sub(".", "", colnames(S))
	lv <- which(substr(colnames(S),1,2) == "LV")
	q <- length(lv)
	#q <- 0
	if (q > 0) {
	 S <- S[-lv,-lv]
	 Sigma <- Sigma[-lv,-lv]
	 p <- ncol(S)
	}

	if (corpcor::is.positive.definite(Sigma) == FALSE){
	 message(" WARNING: Sigma is not a definite positive matrix...\n")
	}
	# Deviance and df for model 1 (fitted model)
	ST <- S %*% solve(Sigma)
	#dev <- n * (sum(diag(ST)) - log(det(ST)) - p) # ggm:::likGau
	d <- svd(S, nu=0, nv=0)$d
	r <- length(which(d < 1e-18))
	d <- d[d > 1e-18]
	dev <- n * (log(det(Sigma)) + sum(diag(ST)) - sum(log(d)) - p)
    devN <- (100/n)*dev
	# Deviance and df for model 0 (null model)
	dev0 <- n * (-sum(log(d))) # n*(-log(det(R)))
	df0 <- p*(p + 1)/2 - p

	# Standardized Root Mean Square Residual (SRMR)
	E <- S - Sigma
	SRMR <- sqrt(mean(E[lower.tri(E, diag = TRUE)]^2))
	# Root Mean Square Error of Approximation (RMSEA)
	RMSEA <- sqrt(max((dev - df), 0)/(df*(n - 1)))
	ncp <- df*n*RMSEA^2
	pRMSEA <- pchisq(dev, df, ncp, lower.tail = FALSE)
	# Comparative Fit Index (CFI)
	CFI <- 1 - (dev - df)/(dev0 - df0)
	# Tucker-Lewis Index (TLI)
	TLI <- (dev0/df0 - dev/df)/(dev0/df0 - 1)
	# ULS Goodness of Fit Index (GFI)
	ULS <- 1 - sum(diag(t(E)%*%E))/sum(diag(t(S)%*%S))

	return(c(dev=dev, df=df, srmr=SRMR, rmsea=RMSEA, n=n, npar=npar, r=r))
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

#' @title Missing edge testing implied by a DAG with Shipley's basis-set
#'
#' @description Compute all the P-values of the d-separation tests
#' implied by the missing edges of a given acyclic graph (DAG).
#' The conditioning set Z is represented, in a DAG, by the union of the
#' parent sets of X and Y (Shipley, 2000). 
#' The results of every test, in a DAG, is then combined using the
#' Fisher’s statistic in an overall test of the fitted model
#' C = -2*sum(log(P-value(k))), where C is distributed as a chi-squared
#' variate with df = 2k, as suggested by Shipley (2000).
#'
#' @param graph A directed graph as an igraph object.
#' @param data A data matrix with subjects as rows and variables as
#' columns.
#' @param MCX2 If TRUE, a Monte Carlo P-value of the combined C test is
#' enabled using the R code of Shipley extracted from 
#' <https://github.com/BillShipley/CauseAndCorrelation>.
#' @param cmax Maximum number of parents set, C. This parameter can be
#' used to perform only those tests where the number of conditioning
#' variables does not exceed the given value. High-dimensional conditional
#' independence tests can be very unreliable. By default, cmax = Inf.
#' @param limit An integer value corresponding to the graph size (vcount)
#' tolerance. Beyond this limit, multicore computation is enabled to
#' reduce the computational burden. By default, \code{limit = 100}.
#' @param verbose If TRUE, Shipley's test results will be showed to
#' screen (default = TRUE).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of three objects: (i) "dag":  the DAG used to perform the Shipley
#' test (ii) "dsep": the data.frame of all d-separation tests over missing edges in
#' the DAG and (iii) "ctest": the overall Shipley's' P-value.
#'
#' @author Mario Grassi \email{mario.grassi@unipv.it}
#'
#' @references
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs.
#' Structural Equation Modeling, 7(2): 206-218.
#' <https://doi.org/10.1207/S15328007SEM0702_4>
#'
#' @examples
#'
#' #\donttest{
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
#' 
#' sem <- SEMrun(alsData$graph, als.npn)
#' C_test <- Shipley.test(sem$graph, als.npn, MCX2 = FALSE)
#' #MC_test <- Shipley.test(sem$graph, als.npn, MCX2 = TRUE)
#' #}
#'
Shipley.test<- function(graph, data, MCX2=FALSE, cmax=Inf, limit=100, verbose=TRUE,...)
{
	# graph to DAG conversion :
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	df1<- vcount(graph)*(vcount(graph)-1)/2-ecount(as.undirected(graph))
	dataY<- as.matrix(data[, nodes])
	S<- cov(dataY)
	n<- nrow(dataY)
	
	if (!is_dag(graph)){
	 cat("WARNING: input graph is not acyclic !\n")
	 cat(" Applying graph -> DAG conversion...\n")
	 dag<- graph2dag(graph, dataY, bap=FALSE) #del cycles & all <->
	 df2<- vcount(dag)*(vcount(dag)-1)/2-ecount(as.undirected(dag))
	 cat(" \nDegrees of freedom:\n Input graph  =", 
            df1, "\n Output graph =", df2, "\n\n")
	}else{
	 dag <- graph
	 df2 <- df1
	}

	# d-separation local tests (B_U) & Shipley's overall pvalue
	cat("d-separation test (basis set) of", df2, "edges...\n")
	dsep<- dsep.test(dag, S, n, cmax=cmax, limit=limit)
	colnames(dsep)[c(3,6,7)]<- c("|Z|", "2.5%", "97.5%")
	#Combining p-values with Fisher's procedure:
	X2<- -2 * sum(log(dsep$p.value + 1E-16))
	df<- 2 * nrow(dsep)
	if (MCX2) {
	 pv<- MCX2(model.df=df, n.obs=n, model.chi.square=X2)[[1]]
	}else{
	 pv<- 1 - pchisq(q = X2, df = df)
	}
	if (verbose) print(data.frame(C_test=X2, df=df, pvalue=round(pv,6)))
			
	return( list(dag=dag, dsep=dsep, ctest=c(X2, df, pv)) )
}

dsep.test<- function(dag, S, n, cmax, limit, ...)
{
	# d-sep (basis set) testing of a DAG
	A <- ifelse(as_adjacency_matrix(as.undirected(dag), sparse=FALSE) == 1, 0, 1)
	ug <- graph_from_adjacency_matrix(A, mode="undirected", diag=FALSE)
	M <- attr(E(ug), "vnames")
	
	local<- function(x) {
	 s <- strsplit(x,"\\|")
	 ed <- c(s[[1]][1], s[[1]][2])
	 pa.r <- V(dag)$name[SEMgraph::parents(dag, ed[1])]
	 pa.s <- V(dag)$name[SEMgraph::parents(dag, ed[2])]
	 dsep <- union(pa.r, pa.s)
     dsep <- setdiff(dsep, ed)
	 B <- c(ed, dsep)
	 pcor <- pcor.test(S, B, n, cmax, H0=0.05)
	 #set <- paste(B[-c(1:2)], collapse=",")
	 return(data.frame(X=B[1], Y=B[2], pcor))
	}

	#message("d-separation test (basis set) of ", length(M), " edges...")
	op<- pbapply::pboptions(type = "timer", style = 2)
	#df<- vcount(dag)*(vcount(dag)-1)/2 - ecount(dag)
	if (vcount(dag) > limit){
	 n_cores <- parallel::detectCores(logical = FALSE)
	 cl <- parallel::makeCluster(n_cores)
	 #cl <- parallel::makeCluster(3)
	 parallel::clusterExport(cl,
	  c("local", "dag", "S", "n"),  envir = environment())
	 SET<- pbapply::pblapply(M, local, cl=cl)
	 parallel::stopCluster(cl)
	}else{
	 SET<- pbapply::pblapply(M, local, cl=NULL)
	}
	SET<- do.call(rbind, lapply(SET, as.data.frame))

	return( SET=na.omit(SET) )
}

pcor.test<- function(S, B, n, cmax, H0, ...)
{
	#Set objects
	q <- length(B)-2
	if (cmax == Inf) cmax <- n - 3
	if (q > cmax) {
	 pcor<- cbind(K=NA, estimate=NA, p.value=NA, lower=NA, upper=NA)
	 return(pcor)
	}
	k <- tryCatch(solve(S[B,B]), error = function(err) NA)
	if (!is.matrix(k)) {
	 pcor <- cbind(K=NA, estimate=NA, p.value=NA, lower=NA, upper=NA)
	 return(pcor)
	}
    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
	if( H0 == 0 ) {
	#Test null H0: r=abs(r(X,Y|Z))=0
	 df <- n - 2 - q
	 tval <- r * sqrt(df)/sqrt(1 - r * r)
	 pval <- 2 * pt(-abs(tval), df)
	}else{
	#Test of not-close fit, H0: r=abs(r(X,Y|Z)) vs. r<.05
	 z <- atanh(r)
	 se <-  1/sqrt(n - 3 - q)
	 pval <- pchisq((z/se)^2, df=1, ncp=(atanh(H0)/se)^2, lower.tail=FALSE)
	}
	lower<- (exp(2*(z - 1.96*se))-1)/(exp(2*(z - 1.96*se))+1)
	upper<- (exp(2*(z + 1.96*se))-1)/(exp(2*(z + 1.96*se))+1)
	pcor<- cbind(K=q, estimate=r, p.value=pval, lower, upper)
	return( pcor )
}

MCX2 <- function (model.df, n.obs, model.chi.square, n.sim = 10000, ...)
{
	#Monte Carlo Chi-square simulator (Author: Bill Shipley) from:
	#devtools::install_github("BillShipley/CauseAndCorrelation")
	# All rights reserved.  See the file COPYING for license terms.
	x <- (-1 + sqrt(1 + 8 * model.df))/2
	if ((x - as.integer(x)) == 0)
	v <- x
	if ((x - as.integer(x)) > 0 & (x - as.integer(x)) < 1) 
	v <- as.integer(x) + 1
	if ((x - as.integer(x)) > 1)return(message("ERROR: check model df !"))
	c.value <- v * (v + 1)/2 - model.df
	MCX2 <- rep(NA, n.sim)
	for (i in 1:n.sim) {
	 dat <- matrix(rnorm(n.obs * v), ncol = v)
	 obs.VCV <- var(dat)
	 model.VCV <- diag(v)
	 diag(model.VCV)[1:c.value] <- diag(obs.VCV)[1:c.value]
	 MCX2[i] <- (n.obs - 1) * (log(det(model.VCV)) + sum(diag(obs.VCV) * 
							(1/diag(model.VCV))) - log(det(obs.VCV)) - v)
	}
	MCprob <- sum(MCX2 >= model.chi.square)/n.sim
	x <- seq(0, max(MCX2))
	theoretical.prob <- dchisq(x, model.df)
	MLprob<- pchisq(model.chi.square, model.df, lower.tail=FALSE)
	
	return(list(MCprob = MCprob, MLprob = MLprob))
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
#' @param limit An integer value corresponding to the size of the
#' extracted acyclic graph. Beyond this limit, switch to Shipley's
#' C-test (Shipley 2000) is enabled to reduce the computational burden.
#' By default, \code{limit = 100}.
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
#' # Nonparanormal(npn) transformation
#' als.npn <- transformData(alsData$exprs)$data
#' 
#' sem <- SEMrun(alsData$graph, als.npn)
#' B_test <- localCI.test(sem$graph, als.npn, verbose = TRUE)
#'
localCI.test<- function(graph, data, bap=FALSE, limit=100, verbose=TRUE, ...)
{
	# Set graph and covariance matrix:
	nodes<- colnames(data)[colnames(data) %in% V(graph)$name]
	graph<- induced_subgraph(graph, vids=which(V(graph)$name %in% nodes))
	dataY<- as.matrix(data[, nodes])
	df1<- vcount(graph)*(vcount(graph)-1)/2-ecount(as.undirected(graph))
	if (vcount(graph) > limit){
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
	imp<- Filter(function(x) length(x$Z) <= 3, imp0)
	XY<- t(sapply(1:length(imp), function(x) imp[[x]][1:2]))
	K<- sapply(1:length(imp), function(x) length(imp[[x]][3]$Z))
	del<- which(duplicated(XY[,1:2]) == TRUE)
	if (length(del)>0) imp<- imp[-del]
	res<- dagitty::localTests(dagi, type = "cis", tests = imp,
	 sample.cov=S, sample.nobs=n, max.conditioning.variables=3, tol=0.05)
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
#'
parameterEstimates<- function(fit, ...)
{
	if (inherits(fit, "lavaan")){
	 est <- lavaan::parameterEstimates(fit)
	 est$lhs <- sub("z", "", est$lhs)
	 est$rhs <- sub("z", "", est$rhs)
	 if ("group" %in% colnames(est)){
	  est <- cbind(est[,-c(4,5)], group=est[,5])
	 }
	}
	else if (!inherits(fit, "lavaan")){
	 if (length(fit) != 2){
	  est0<- fit$parameterEstimates
	  sel<- which(est0$rhs == "group")
	  if(length(sel) == 0){
	   est <- est0
	  }else{
	   est <- est0[sel,]
	   est <- rbind(est, est0[-sel,])
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
#'
summary.RICF <- function(object, ...)
{
	.local <- function(object) {
		it <- object$it
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
			print(L[[l]])
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
#'
summary.GGM <- function(object, ...)
{
	.local <- function(object) {
		it <- 0
		t <- object$fitIdx[6]
		n <- object$fitIdx[5]
		dev <- round(object$fitIdx[1], 3)
		df <- object$fitIdx[2]
		srmr <- round(object$fitIdx[3], 3)

		cat(paste0("GGM (de-biased nodewise L1) solver ended normally after ", it,
		           " iterations"), "\n\n")
		cat(paste0("  Estimator                                       ML"),
		           "\n")
		cat(paste0("  Optimization method                             GGM"),
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

		B<- object$parameterEstimates
		cat("Regressions:\n\n")
		print(B)
		V<- diag(object$Psi)
		names(V)<- rownames(object$Psi)
		cat("\nVariances:\n\n")
		print(round(V,3))
	}
	.local(object)
}
