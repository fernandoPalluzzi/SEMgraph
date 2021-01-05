#  SEMgraph library
#  Copyright (C) 2019-2021 Fernando Palluzzi; Mario Grassi
#  e-mail: <fernando.palluzzi@gmail.com>
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
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #


SEMmodel<- function(ig, nodes, group, ...)
{
	# Set from-to-matrix representation of gene-gene links
	ftm <- as_data_frame(ig)
	if (is.directed(ig) & sum(which_mutual(ig)) > 0) {
		dg <- ig - E(ig)[which_mutual(ig)]
		ug <- as.undirected(ig-E(ig)[!which_mutual(ig)])
		xx <- V(dg)$name[igraph::degree(dg, mode="in") == 0]
		ug <- ug - V(ug)[V(ug)$name %in% xx]
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

#' @title Fit a graph as a Structural Equation Model (SEM)
#'
#' @description SEMrun converts a (directed, undirected, or mixed) graph 
#' to a SEM and fits it. If a binary group variable (i.e., case/control) 
#' is present, node-level or edge-level perturbation is evaluated. SEMrun 
#' can handle loop-containing models, although multiple links between 
#' the same two nodes (including self-loops and mutual interactions) and 
#' bows (i.e., a directed and a bidirected link between two nodes) are 
#' not allowed.
#'
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param fit A numeric value indicating the SEM fitting mode. 
#' If fit = 0 (default), no group effect is considered. If fit = 1, a 
#' "common" model is used to evaluate group effects on graph nodes. 
#' If fit = 2, a two-group model is used to evaluate group effects on 
#' graph edges.
#' @param algo MLE method used for SEM fitting. If algo = "lavaan" (default), 
#' the SEM will be fitted using the NLMINB solver from \code{lavaan} R 
#' package, with standard errors derived from the expected Fisher information 
#' matrix. If algo = "ricf", the model is fitted via residual iterative 
#' conditional fitting (RICF; Drton et al. 2009). 
#' If algo = "ggm", model fitting is based on constrained Gaussian Graphical 
#' Modeling (GGM) and de-sparsified glasso estimator (Williams, 2020).
#' @param start Starting value of SEM parameters for \code{algo} = "lavaan". 
#' If start is NULL (default), the algorithm will determine the starting 
#' values. If start is a numeric value, it will be used as a scaling factor 
#' for the edge weights in the graph object (graph attribute E(graph)$weight). 
#' For instance, a scaling factor is useful when weights have fixed values 
#' (e.g., 1 for activated, -1 for repressed, and 0 for unchanged interaction). 
#' Fixed values may compromise model fitting, and scaling them is a safe 
#' option to avoid this problem. As a rule of thumb, to our experience, 
#' start = 0.1 generally performs well with {-1, 0, 1} weights.
#' @param limit An integer value corresponding to the network size
#' (i.e., number of nodes). Beyond this limit, the execution under 
#' \code{algo} = "lavaan" will be ridirected to \code{algo} = "ricf", if 
#' fit is either 0 or 1, or to algo = "ggm", if fit = 2. This redirection 
#' is necessary to reduce the computational demand of standard error 
#' estimation by lavaan. Increasing this number will enforce lavaan 
#' execution when \code{algo} = "lavaan".
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
#' coefficients) will be estimated. This will also lead to the estimation 
#' of a set of aggregated group effects, if \code{algo} = "ricf" (see 
#' \code{\link[SEMgraph]{SEMgsa}}).
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
#' model, as absent (\code{fit} = 0), node perturbation (\code{fit} = 1), 
#' or edge perturbation (\code{fit} = 2). When \code{fit} = 1, the group 
#' is modeled as an exogenous variable, influencing all the other graph 
#' nodes. When \code{fit} = 2, SEMrun estimates the differences 
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
#' \item "fit", SEM fitted lavaan, ricf, or ggmncv object, depending 
#' on the MLE method specified by the \code{algo} argument;
#' \item "gest" or "dest", a data.frame of node-specific ("gest") 
#' or edge-specific ("dest") group effect estimates and P-values;
#' \item "model", SEM model as a string if \code{algo} is "lavaan", and 
#' NULL otherwise;
#' \item "graph", the induced subgraph of the input network mapped on
#' data variables. Graph edges (i.e., direct effects) with P-value < 0.05 
#' will be highlighted in red (beta > 0) or blue (beta < 0). If a group 
#' vector is given, nodes with significant group effect (P-value < 0.05) 
#' will be red-shaded (beta > 0) or lightblue-shaded (beta < 0);
#' \item "dataXY", input data subset mapping graph nodes, plus group at
#' the first column (if no group is specified, this column will take NA
#' values).
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom stats cor
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom GGMncv constrained inference ggm_compare
#' @importFrom ggm fitAncestralGraph
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom gdata unmatrix
#' @export
#'
#' @references
#'
#' Pearl J (1998). Graphs, Causality, and Structural Equation Models.
#' Sociological Methods & Research., 27(2):226-284. 
#' https://doi.org/10.1177/0049124198027002004
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
#' Drton M, Eichler M, Richardson TS (2009). Computing Maximum Likelihood 
#' Estimated in Recursive Linear Models with Correlated Errors. 
#' Journal of Machine Learning Research, 10(Oct): 2329-2348. 
#' http://www.jmlr.org/papers/volume10/drton09a/drton09a.pdf
#'
#' Larson JL and Owen AB (2015). Moment based gene set tests. BMC
#' Bioinformatics, 16: 132. https://doi.org/10.1186/s12859-015-0571-7
#'
#' Grassi M, Palluzzi F (2021). SEMgraph: An R Package for Causal Network 
#' Analysis of High-Throughput Data with Structural Equation Models. 
#' xxxxx x(x): xxxxx. https://doi.org/xxxxx
#' 
#' Williams D (2020). GGMncv: Gaussian Graphical Models with Non-Convex 
#' Penalties. R package version 1.1.0.
#' https://CRAN.R-project.org/package=GGMncv
#'
#' @seealso See \code{\link[ggm]{fitAncestralGraph}} for RICF algorithm 
#' details, \code{\link[flip]{flip}} for randomization P-values, and 
#' \code{\link[GGMncv]{constrained}} for constrained GGM, and 
#' \code{\link[GGMncv]{inference}} for de-sparsified P-values.
#'
#' @examples
#' 
#' #### Model fitting (no group effect)
#' 
#' sem <- SEMrun(graph = sachs$graph, data = log(sachs$pkc))
#' 
#' # Fitting summaries
#' summary(sem$fit)
#' print(sem$gest)
#' head(parameterEstimates(sem$fit))
#' 
#' # Graphs
#' gplot(sem$graph, main="node differences")
#' plot(sem$graph, layout=layout.circle, main="node differences")
#' 
#' 
#' #### Model fitting (common model, group effect on nodes)
#' 
#' sem <- SEMrun(graph = sachs$graph, data = log(sachs$pkc), group = sachs$group)
#' 
#' # Fitting summaries
#' summary(sem$fit)
#' print(sem$gest)
#' head(parameterEstimates(sem$fit))
#' 
#' # Graphs
#' gplot(sem$graph, main="node differences")
#' plot(sem$graph, layout=layout.circle, main="node differences")
#'
#' 
#' #### Two-group model fitting (group effect on edges)
#' 
#' sem2 <- SEMrun(graph = sachs$graph, data = log(sachs$pkc),
#'                group = sachs$group,
#'                fit = 2)
#' 
#' # Summaries
#' summary(sem2$fit) # class lavaan
#' print(sem2$dest)
#' head(parameterEstimates(sem2$fit))
#' 
#' # Graphs
#' gplot(sem2$graph, main = "Edge differences")
#' plot(sem2$graph, layout = layout.circle, main = "Edge differences")
#' 
#' ## NOT run ## {
#' 
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#'
#' # Multicore computation
#' 
#' j <- which(names(kegg.pathways) == "Wnt signaling pathway")
#' 
#' sem2 <- SEMrun(graph = kegg.pathways[j], data = als.npn,
#'                group = alsData$group,
#'                fit = 2,
#'                algo = "ricf")
#' 
#' head(parameterEstimates(sem2$fit))
#' head(sem2$dest)
#' # gplot(sem2$graph, main = "Edge differences")
#' # plot(sem2$graph, layout = layout.circle, main = "Edge differences")
#' 
#' # Edges null randomization distribution
#' r_dest <- t(na.omit(sem2$r_dest[,-c(1:3)]))
#' beta <- scale(r_dest)
#' colnames(beta) <- paste0(sem2$r_dest$lhs, " <- ", sem2$r_dest$rhs)
#' pairwiseMatrix(beta)
#' 
#' ## }
#'
SEMrun <- function(graph, data, group = NULL, fit = 0, algo = "lavaan",
                   start = NULL, limit = 100, ...)
{
	if (is.null(group) & fit != 0) fit <- 0
	if (!is.null(group) & fit == 0) fit <- 1
	
	if (fit == 0) {
		if (algo == "lavaan") {
			return(fit = SEMfit(graph = graph, data = data, group = NULL,
			                    start = start,
			                    limit = limit))
		} else if (algo == "ggm") {
			return(fit = SEMggm(graph = graph, data = data, group = NULL))
		} else if (algo == "ricf") {
			return(fit = SEMricf(graph = graph, data = data, group = NULL))
		}
	}
	
	if (fit == 1) {
		if (algo == "lavaan") {
			return(fit = SEMfit(graph = graph, data = data, group = group,
			                    start = start,
			                    limit = limit))
		} else if (algo == "ggm") {
			return(fit = SEMggm(graph = graph, data = data, group = group))
		} else if( algo == "ricf" ) {
			return(fit = SEMricf(graph = graph, data = data, group = group,
			                     n_rep = 5000))
		}
	}
	
	if (fit == 2) {
		if (algo == "lavaan") {
			return(fit = SEMfit2(graph = graph, data = data, group = group,
			                     start = start,
			                     limit = limit))
		} else if (algo == "ggm") {
			return(fit = SEMggm2(graph = graph, data = data, group = group))
		} else if (algo == "ricf") {
			return(fit = SEMricf2(graph = graph, data = data, group = group,
			                      n_rep = 1000))
		}
	}
}

SEMfit <- function(graph, data, group = NULL, start = NULL, limit = 100, ...)
{
	# Model fitting with GGM algo if n.nodes > limit
	if (vcount(graph) > limit) {
		cat("WARNING: very large input graph (>", limit, "nodes) !\n")
		cat(" RICF solver activated...\n\n")
		return(fit = SEMricf(graph = graph, data = data, group = group))
	}
	
	# Set data and model objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- as.matrix(data[, nodes])
	colnames(dataY) <- paste0("z", nodes)
	dataXY <- cbind(group, dataY)
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	if (corpcor::is.positive.definite(cor(dataXY)[1:p, 1:p])) {
		covXY <- cor(dataXY)[1:p, 1:p]
	} else {
		covXY <- corpcor::cor.shrink(dataXY, verbose = TRUE)[1:p, 1:p]
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
	suppressWarnings(fit <- lavaan(model, sample.cov = covXY,
	                               sample.nobs = n, se = "standard",
	                               fixed.x = TRUE, int.ov.free = TRUE,
	                               auto.var = TRUE, auto.cov.y = FALSE,
	                               control = list(abs.tol = 1e-20,
	                                              rel.tol = 1e-10)))
	
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
	if (!is.null(group)) {
		gest <- est[1:(p - 1),]
		pval1 <- Brown.test(x = dataY, p = gest$pvalue, theta = gest$est,
		                    tail = "positive")
		pval2 <- Brown.test(x = dataY, p = gest$pvalue, theta = gest$est,
		                    tail = "negative")
		cat("Brown's combined P-value of node activation:", pval1, "\n\n")
		cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	} else {
		gest <- NULL
	}
	
	# Output objects
	ig <- colorGraph(est = est, graph = ig, group = group, alpha = 0.05)
	if (is.null(group)) dataXY <- cbind(group = rep(NA, n), dataXY)
	colnames(dataXY) <- gsub("z", "", colnames(dataXY))
	
	return(list(fit = fit, gest = gest, model = model, graph = ig,
	            dataXY = dataXY))
}

SEMfit2 <- function(graph, data, group, start = NULL, limit = 100, ...)
{
	# Model fitting with GGM algo if n.nodes > limit
	if (vcount(graph) > limit) {
		cat("WARNING: input graph is very large ( >", limit, "nodes ) !\n")
		cat(" GGM (constrained) solver activated...\n\n")
		return(fit = SEMggm2(graph = graph, data = data, group = group))
	}
	
	# Set data and model objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- as.matrix(data[, nodes])
	colnames(dataY) <- paste0("z", nodes)
	p <- ncol(dataY)
	
	# Covariances for cases (GROUP 1)
	data1 <- dataY[group == 1,]
	n1 <- nrow(data1)
	if (corpcor::is.positive.definite(cor(data1)[1:p,1:p])) {
		cov1 <- cor(data1)[1:p, 1:p]
	} else {
		cov1 <- corpcor::cor.shrink(data1, verbose = TRUE)[1:p, 1:p]
	}
	
	# Covariances for controls (GROUP 0)
	data0 <- dataY[group == 0,]
	n0 <- nrow(data0)
	if (corpcor::is.positive.definite(cor(data0)[1:p, 1:p])) {
		cov0 <- cor(data0)[1:p, 1:p]
	} else {
		cov0 <- corpcor::cor.shrink(data0, verbose = TRUE)[1:p, 1:p]
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
	suppressWarnings(fit <- lavaan(model, sample.cov = list(cov0,cov1),
	                               sample.nob = list(n0, n1),
	                               se = "standard", fixed.x = TRUE,
	                               int.ov.free = TRUE, auto.var = TRUE,
	                               auto.cov.y = FALSE,
	                               control = list(abs.tol = 1e-20,
	                                              rel.tol = 1e-10)))
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
	d_se <- sqrt(est1$se^2 + est0$se^2)
	pvalue <- 2*(1-pnorm(abs(d_est/d_se)))
	d_lower <- d_est - 1.96*d_se
	d_upper <- d_est + 1.96*d_se
	pval1 <- Brown.test(x = NULL, p = pvalue, theta = d_est, tail = "positive")
	pval2 <- Brown.test(x = NULL, p = pvalue, theta = d_est, tail = "negative")
	cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
	cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
	
	# Output objects
	dest <- cbind(est0[, 1:3], d_est, d_se, d_z = d_est/d_se, pvalue,
	              d_lower, d_upper)
	dest <- data.frame(lapply(dest,
	                   function(y) if(is.numeric(y)) round(y, 3) else y))
	ig <- colorGraph(est = dest, graph = ig, group = NULL, alpha = 0.05)
	dataXY <- cbind(c(rep(1, n1), rep(0, n0)), rbind(data1, data0))
	colnames(dataXY) <- gsub("z", "", colnames(dataXY))
	
	return(list(fit = fit, dest = dest, model = model, graph = ig,
	            dataXY = dataXY))
}

SEMricf <- function(graph, data, group = NULL, random.x = FALSE, n_rep = 5000, ...)
{
	# Set data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	dataXY <- cbind(group,dataY)
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	if (corpcor::is.positive.definite(cor(dataXY)[1:p, 1:p])) {
		covXY <- cor(dataXY)[1:p, 1:p]
	} else {
		covXY <- corpcor::cor.shrink(dataXY, verbose = TRUE)[1:p, 1:p]
	}
	
	# Set graph objects
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	E(ig)$weight <- ifelse(which_mutual(ig), 100, 1)
	A <- igraph::as_adj(ig, attr = "weight", sparse = FALSE)[nodes, nodes]
	if (is.null(group) & random.x == TRUE) {
		Vx <- which(colSums(A) == 0)
		A[Vx, Vx] <- 100
		diag(A) <- 0
	}
	if (!is.null(group)) {
		A <- cbind(rep(0, p), rbind(rep(1, p - 1), A))
		colnames(A)[1] <- rownames(A)[1] <- "group"
	}
	
	# Fit model with RICF algorithm
	fit <- ggm::fitAncestralGraph(amat = A, S = covXY, n = n, tol = 1e-6)
	cat(paste0("RICF solver ended normally after ", fit$it, " iterations"),
	           "\n\n")
	idx <- fitIndices(ig, n, covXY, fit$Shat)
	cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3], 7), "\n\n")
	
	est <- parameterEstimates.RICF(fit)
	if (!is.null(group)) {
		gest <- gest.RICF(fit = fit, data = dataY, group = group,
		                  n_rep = n_rep)
		pval1 <- Brown.test(x = dataY, p = gest[[1]][-c(1:3), 4],
							theta = gest[[1]][-c(1:3), 2],
							tail = "positive")
		pval2 <- Brown.test(x = dataY, p = gest[[1]][-c(1:3), 4],
							theta = gest[[1]][-c(1:3), 2],
							tail = "negative")
		cat("Brown's combined P-value of node activation:", pval1, "\n\n")
		cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
		ig <- colorGraph(est = gest[[1]], ig, group, alpha = 0.05)
	} else {
		gest <- list(NULL, NULL)
	}
	
	# Output objects
	fit <- list(ricf = fit, fitIdx = idx, parameterEstimates = est)
	if (is.null(group)) {
		dataXY <- cbind(group = rep(NA, n), dataXY)
	}
	
	return(list(fit = fit, gest = gest[[1]], model = NULL, graph = ig,
	            dataXY = dataXY, r_gest = gest[[2]]))
}

SEMricf2 <- function(graph, data, group, random.x = FALSE, n_rep = 1000,
                     n_cores = NULL, ...)
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
	fit1 <- quiet(SEMricf(graph, data1, group = NULL, random.x = FALSE, n_rep))
	est1 <- fit1$fit$parameterEstimates$Reg
	
	# Fitting RICF for group = 0
	fit0 <- quiet(SEMricf(graph, data0, group = NULL, random.x = FALSE, n_rep))
	est0 <- fit0$fit$parameterEstimates$Reg
	
	# Two-group fit indices
	it <- fit1$fit$ricf$it + fit0$fit$ricf$it
	cat(paste0("RICF solver ended normally after ", it, " iterations"), "\n\n")
	srmr <- (n1/n)*fit1$fit$fitIdx[3] + (n0/n)*fit0$fit$fitIdx[3]
	dev <- (n1/n)*fit1$fit$fitIdx[1] + (n0/n)*fit0$fit$fitIdx[1]
	df <- fit1$fit$fitIdx[2]
	cat("deviance/df:", dev/df , " srmr:", round(srmr, 7), "\n\n")
	
	# Edge differences based on group randomization with n_rep flips
	est <- list(Group_1 = est1, Group_0 = est0)
	dest <- dest.RICF(ig, dataY, group, est, n_rep, n_cores)
	pval1 <- Brown.test(x = NULL, p = dest[[1]][, 5],
	                    theta = dest[[1]][, 4],
	                    tail = "positive")
	pval2 <- Brown.test(x = NULL, p = dest[[1]][, 5],
	                    theta = dest[[1]][, 4],
	                    tail = "negative")
	cat("Brown's combined P-value of node activation:", pval1, "\n\n")
	cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	
	# Output objects
	fit <- list(fit1$fit[[1]], fit0[[1]], parameterEstimates = est)
	ig <- colorGraph(est = dest[[1]], graph = ig, group = NULL, alpha = 0.05)
	dataXY <- cbind(c(rep(1, n1), rep(0, n0)), rbind(data1, data0))
	
	return(list(fit = fit, dest = dest[[1]], model = NULL, graph = ig,
	            dataXY = dataXY, r_dest = dest[[2]]))
}

gest.RICF <- function(fit, data, group, n_rep, ...)
{
	# Permutation pvalues of group effects on graph nodes
	p <- ncol(data)
	B <- (diag(p + 1) - fit$Bhat)[-1, -1]
	Psi <- fit$Ohat[-1, -1]
	Y <- data[, colnames(B)]
	EE <- eigen(Psi)
	R <- EE$vectors%*%diag(1/sqrt(EE$values))%*%t(EE$vectors)
	D <- as.matrix(Y)%*%R%*%rep(1/sqrt(p), p)
	#D <- as.matrix(Y)%*%t(diag(p) - B)%*%rep(1/sqrt(p), p)
	A <- as.matrix(Y)%*%B%*%rep(1/sqrt(p), p)
	E <- as.matrix(Y)%*%t(B)%*%rep(1/sqrt(p), p)
	Z <- as.matrix(Y)%*%t(diag(p) - B)*sd(group)
	perm <- list(B = n_rep + 1, seed = 123)
	gest <- flip::flip(cbind(D, A, E, Z), ~group, perms = perm,
	                   statTest = "t")
	#pval <- flip::npc(gest[-c(1:3)], "fisher")@res$p
	
	# N(0,1) approximate pvalues
	aveT <- apply(gest@permT[-1,], 2, mean)
	sdT <- apply(gest@permT[-1,], 2, sd)
	z <- abs(gest@permT[1,]- aveT)/sdT
	gest@res[, 4] <- 2*(1 - pnorm(z))
	colnames(gest@res)[4] <- "pvalue"
	rownames(gest@res)[1:3] <- c("D", "A", "E")
	
	return(list(gest = gest@res, r_gest = gest@permT))
}

dest.RICF <- function (graph, data, group, est, n_rep, n_cores, ...)
{
	# Randomized pvalues of the beta (psi) differences between groups
	r_dest <- modelRandomization.RICF(graph, data, group, n_rep, n_cores)
	r_names <- paste0(est[[1]]$lhs, ":", est[[1]]$rhs)
	r_dest <- r_dest[match(r_names, rownames(r_dest)),]
	d_est <- est[[1]]$est - est[[2]]$est
	pv <- vector()
	pvalue <- vector()
	for (j in 1:nrow(r_dest)) {
		R_dest <- r_dest[j,]
		O_dest <- d_est[j]
		pv1 <- (sum(R_dest >= O_dest) + 1)/(n_rep + 1)
		pv2 <- (sum(R_dest <= O_dest) + 1)/(n_rep + 1)
		pv[j] <- 2*min(pv1, pv2)
		
		# N(0,1) approximate pvalues
		aveT <- mean(as.numeric(R_dest))
		sdT <- sd(as.numeric(R_dest))
		z <- abs(O_dest - aveT)/sdT
		pvalue[j] <- 2*(1 - pnorm(z))
	}
	dest <- cbind(est[[1]][, c(1:3)], d_est, pvalue)
	r_dest <- cbind(est[[1]][, c(1:3)], r_dest)
	rownames(r_dest) <- NULL
	
	return(list(dest = dest, r_dest = r_dest))
}

modelRandomization.RICF <- function(graph, data, group, n_rep, n_cores, ...)
{
	# Set data and progress bar
	require(foreach)
	if (is.null(n_cores)) n_cores <- parallel::detectCores()
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	n <- nrow(dataY)
	p <- ncol(dataY)
	opts <- list(progress = pb(n_rep))
	
	# Model fitting
	message("Model Randomization ...")
	set.seed(1234)
	idx <- lapply(1:n_rep, function(x) sample(1:n, size = n, replace = FALSE))
	cl <- parallel::makeCluster(n_cores)
	doSNOW::registerDoSNOW(cl)
	#doParallel::registerDoParallel(cl)
	r <- foreach::foreach(i = 1:n_rep, .combine = cbind,
	                      .options.snow = opts) %dopar% {
	     dataY <- dataY[idx[[i]],]
	     data0 <- dataY[group == 0,]
	     data1 <- dataY[group == 1,]
	     fit0 <- quiet(try(SEMricf(graph, data = data0, random.x = FALSE,
	                       n_rep = 0)$fit))
	     fit1 <- quiet(try(SEMricf(graph, data = data1, random.x = FALSE,
	                       n_rep = 0)$fit))
	     B0 <- gdata::unmatrix(diag(p) - fit0$ricf$Bhat, byrow = FALSE)
	     B1 <- gdata::unmatrix(diag(p) - fit1$ricf$Bhat, byrow = FALSE)
	     return(B1[B1 != 0] - B0[B0 != 0])
	}
	parallel::stopCluster(cl)
	cat("\n")
	return(r)
}

parameterEstimates.RICF <- function(fit, ...)
{
	# Output  parameter estimates
	S <- fit$S
	H <- 1*lower.tri(S, diag = TRUE)
	B <- gdata::unmatrix(diag(nrow(S)) - fit$Bhat, byrow = FALSE)
	B <- as.matrix(B[which(B != 0)])
	O <- gdata::unmatrix(fit$Ohat*H, byrow = FALSE)
	O <- as.matrix(O[which(O != 0)])
	P <- round(rbind(B, O), 3)
	
	est <- NULL
	for(j in 1:length(P)) {
		s <- strsplit(rownames(P)[j], ":")
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
	
	reg <- est[which(est$op == "~"),]
	cov <- est[which(est$op == "~~"),]
	sel <- which(cov$lhs == cov$rhs)
	var <- cov[sel,]
	if (nrow(cov[-sel,]) != 0) cov <- cov[-sel,] else cov <- NULL
	reg <- reg[order(reg$lhs),]
	if(!is.null(cov)) cov <- cov[order(cov$lhs),]
	var <- var[order(var$lhs),]
	
	return(list(Reg = reg, Cov = cov, Var = var))
}

SEMggm <- function(graph, data, group = NULL, method = "none", alpha = 0.05,
                   ...)
{
	# Set data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[, nodes]
	dataXY <- cbind(group, dataY)
	n <- nrow(dataXY)
	p <- ncol(dataXY)
	
	if (corpcor::is.positive.definite(cor(dataXY)[1:p, 1:p])) {
		covXY <- cor(dataXY)[1:p, 1:p]
	} else {
		covXY <- corpcor::cor.shrink(dataXY, verbose = TRUE)[1:p, 1:p]
	}
	
	# Set graph objects
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	adj <- as_adj(as.undirected(ig), sparse = FALSE)[nodes, nodes]
	E(ig)$weight <- ifelse(which_mutual(ig), 100, 1)
	dadj <- as_adj(ig, attr = "weight", sparse = FALSE)[nodes, nodes]
	
	if (!is.null(group)) {
		adj <- cbind(c(0, rep(1, p - 1)), rbind(rep(1, p - 1), adj))
		colnames(adj)[1] <- rownames(adj)[1] <- "group"
		dadj <- cbind(rep(0, p), rbind(rep(1, p - 1), dadj))
		colnames(dadj)[1] <- rownames(dadj)[1] <- "group"
	}
	
	# Constrained GGM
	cggm <- GGMncv::constrained(covXY, adj)
	cat(paste0("GGM (constrained) solver ended normally after ", 0,
	           " iterations"), "\n\n")
	idx <- fitIndices(ig, n, covXY, cggm$Sigma)
	cat("deviance/df:", idx[1]/idx[2], " srmr:", round(idx[3], 7), "\n\n")
	Theta <- cggm$Theta
	rownames(Theta) <- colnames(Theta) <- colnames(dataXY)
	fit <- list(Theta = Theta, n = n, R = covXY)
	class(fit) <- c("ggmncv", "default")
	
	# Beta & Psi matrices
	betas <- function(Theta) {
		-1*sapply(1:p, function(x) Theta[x,]/Theta[x, x])
	}
	B <- ifelse(dadj == 1, betas(Theta), 0)
	colnames(B) <- colnames(dataXY)
	Dt <- diag(1/sqrt(diag(Theta)))
	O <- ifelse(dadj == 100, Dt%*%(-1*Theta)%*%Dt, 0)
	#O <- ifelse(dadj == 100, cggm$Sigma, 0)
	diag(O) <- apply(scale(dataXY) %*% (diag(p) - B), 2, var)
	colnames(O) <- colnames(dataXY)
	
	# Edge pvalues based on the de-sparsified precision matrix
	dggm <- suppressWarnings(GGMncv::inference(fit, method = method,
	                         alpha = alpha))
	pvB <- ifelse(dadj == 1, dggm$corrected, 0)
	rownames(pvB) <- colnames(pvB) <- colnames(dataXY)
	pvO <- ifelse(dadj == 100, dggm$corrected, 0)
	rownames(pvO) <- colnames(pvO) <- colnames(dataXY)
	
	est <- parameterEstimates.GGM(fit = list(B = B, O = O, pvB = pvB,
	                              pvO = pvO))
	if (!is.null(group)) {
		gest <- est$Reg[est$Reg$rhs == "group",]
		pval1 <- Brown.test(x = dataY, p = gest$pvalue, theta = gest$est,
		                    tail = "positive")
		pval2 <- Brown.test(x = dataY, p = gest$pvalue, theta = gest$est,
		                    tail = "negative")
		cat("Brown's combined P-value of node activation:", pval1, "\n\n")
		cat("Brown's combined P-value of node inhibition:", pval2, "\n\n")
	} else {
		gest <- NULL
	}
	
	# Output objects
	fit <- list(cggm = fit, Beta = B, Psi = O, fitIdx = idx,
	            parameterEstimates = est)
	ig <- colorGraph(est = est$Reg, graph = ig, group = group, alpha = 0.05)
	if (is.null(group)) dataXY <- cbind(group = rep(NA, n), dataXY)
	
	return(list(fit = fit, gest = gest, model = NULL, graph = ig,
	            dataXY = dataXY))
}

SEMggm2 <- function(graph, data, group, method = "none", alpha = 0.05, ...)
{
	# Set graph and data objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	E(ig)$weight <- ifelse(which_mutual(ig), 100, 1)
	dadj <- as_adj(ig, attr = "weight", sparse = FALSE)[nodes, nodes]
	dataY <- data[, nodes]
	data1 <- dataY[group == 1,]
	data0 <- dataY[group == 0,]
	n1 <- nrow(data1)
	n0 <- nrow(data0)
	n <- n1 + n0
	p <- ncol(dataY)
	
	# Constrained GGM for group = 1
	cggm1 <- quiet(SEMggm(graph, data1, group = NULL, method = method,
	                      alpha = alpha))
	fit1 <- cggm1$fit$cggm
	
	# constrained GGM for group = 0
	cggm0 <- quiet(SEMggm(graph, data0, group = NULL, method = method,
	                      alpha = alpha))
	fit0 <- cggm0$fit$cggm
	
	# Two-group fit indices
	cat(paste0("GGM (constrained) solver ended normally after ", 0,
	           " iterations"), "\n\n")
	srmr <- (n1/n)*cggm1$fit$fitIdx[3] + (n0/n)*cggm0$fit$fitIdx[3]
	dev <- (n1/n)*cggm1$fit$fitIdx[1] + (n0/n)*cggm0$fit$fitIdx[1]
	df <- cggm1$fit$fitIdx[2]
	cat("deviance/df:", dev/df , " srmr:", round(srmr, 7), "\n\n")
	
	# Edge differences based on the de-sparsified precision matrix
	dggms <- suppressWarnings(GGMncv::ggm_compare(fit1, fit0,
	                                              method = method,
	                                              alpha = alpha))
	d_est <- cggm1$fit$Beta - cggm0$fit$Beta
	d_pv <- ifelse(dadj == 1, dggms$corrected, 0)
	rownames(d_pv) <- colnames(d_pv) <- colnames(dataY)
	
	dest <- parameterEstimates.GGM(fit = list(B = d_est, O = diag(p),
	                               pvB = d_pv, pvO = diag(p)))[[1]]
	pval1 <- Brown.test(x = NULL, p = dest$pvalue, theta = dest$d_est,
	                    tail = "positive")
	pval2 <- Brown.test(x = NULL, p = dest$pvalue, theta = dest$d_est,
	                    tail = "negative")
	cat("Brown's combined P-value of edge activation:", pval1, "\n\n")
    cat("Brown's combined P-value of edge inhibition:", pval2, "\n\n")
	
	# Output objects
	est <- list(Group_1 = cggm1$fit[[5]], Group_0 = cggm0$fit[[5]])
	fit <- list(fit1 = cggm1[-5], fit0 = cggm0[-5], parameterEstimates = est)
	ig <- colorGraph(est = dest, graph = ig, group = NULL, alpha = 0.05)
	dataXY <- cbind(c(rep(1, n1), rep(0, n0)), rbind(data1, data0))
	
	return(list(fit = fit, dest = dest, model = NULL, graph = ig,
	            dataXY = dataXY))
}

parameterEstimates.GGM <- function(fit, ...)
{
	# Output parameter estimates
	B <- gdata::unmatrix(fit$B, byrow = TRUE)
	B <- as.matrix(B[which(B != 0)])
	var <- diag(fit$O)
	O <- ifelse(lower.tri(fit$O), fit$O, 0)
	rownames(O) <- colnames(O) <- rownames(fit$O)
	O <- gdata::unmatrix(O, byrow=TRUE)
	O <- as.matrix(O[which(O != 0)])
	P <- round(rbind(B, O), 3)
	PvB <- gdata::unmatrix(fit$pvB, byrow = TRUE)
	PvB <- as.matrix(PvB[which(PvB != 0)])
	PvO <- gdata::unmatrix(fit$pvO, byrow = TRUE)
	PvO <- as.matrix(PvO[which(PvO != 0)])
	Pv <- rbind(PvB, PvO)
	
	est <- NULL
	for(j in 1:length(P)) {
		s <- strsplit(rownames(P)[j], ":")
		lhs <- s[[1]][2]
		rhs <- s[[1]][1]
		if (j <= length(B)) {
			est <- rbind(est, data.frame(lhs, op = "~", rhs, est = P[j],
			             pvalue = Pv[j],
			             stringsAsFactors = FALSE))
		} else {
			est <- rbind(est, data.frame(lhs, op = "~~", rhs, est = P[j],
			             pvalue = Pv[j],
			             stringsAsFactors = FALSE))
		}
	}
	
	reg <- est[which(est$op == "~"),]
	reg <- reg[order(reg$lhs),]
	if (sum(var) == length(var)) {
		colnames(reg)[4] <- "d_est"
		cov <- NULL
		var <- NULL
	} else {
		cov <- est[which(est$op == "~~"),]
		if (nrow(cov) == 0) cov <- NULL
		var <- data.frame(lhs = names(var), op = "~~", rhs = names(var),
		                  est = round(var, 3))
		var <- var[order(var$lhs),]
		rownames(var) <- NULL
	}
	
	return(list(Reg = reg, Cov = cov, Var = var))
}

fitIndices <- function(graph, n, S, Sigma, ...)
{
	p <- nrow(S)
	q <- ecount(graph) - sum(which_mutual(graph))
	E <- S - Sigma
	ST <- S %*% solve(Sigma)
	
	# Deviance and df for model 1 (fitted model)
	dev <- n*(sum(diag(ST)) - log(det(ST)) - p)
	df <- p*(p + 1)/2 - (p + q)
	
	# Deviance and df for model 0 (null model)
	dev0 <- n*(sum(diag(S)) - log(det(S)) - p)
	df0 <- p*(p + 1)/2 - p
	
	# Standardized Root Mean Square Residual (SRMR)
	SRMR <- sqrt(mean(E[lower.tri(E, diag = TRUE)]^2))
	
	# Root Mean Square Error of Approximation(RMSEA)
	RMSEA <- sqrt(max((dev - df),0)/(df*(n - 1)))
	
	# Comparative Fit Index (CFI)
	CFI <- 1 - (dev - df)/(dev0 - df0)
	
	# Tucker-Lewis Index (TLI)
	TLI <- (dev0/df0 - dev/df)/(dev0/df0 - 1)
	
	#ULS Goodness of Fit Index (GFI) 
	#ULS <- 1 - sum(diag(t(E)%*%E))/sum(diag(t(S)%*%S))
	
	return(c(dev = dev, df = df, srmr = SRMR, rmsea = RMSEA, n = n, t = p + q))
}

#' @title Bow-free covariance search
#'
#' @description Search for new bow-free interactions. SEMbap does an 
#' exhaustive missing interaction search, using Shipley's d-separation test 
#' (Shipley, 2000) to add a bidirected edge (i.e., a bow-free covariance) 
#' between pairs of disconnected nodes (X, Y) if there is a significant 
#' association between them, given a set Z of conditioning nodes. 
#' Here, the "bow-free" condition states that there cannot be both 
#' directed and bidirected connections between the same pair of nodes 
#' (Nowzohour, 2017).
#'
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes (variables).
#' @param gnet Reference "global" network as an igraph object.
#' @param d An integer value indicating the maximum length of indirect
#' interactions between pairs of nodes. If d = 1, direct interactions 
#' between nodes will be searched in the reference interactome (if given). 
#' If d > 1, indirect interactions of length d or shorter (i.e., with at 
#' most d - 1 connectors) between bow-free nodes will be searched. 
#' Setting d = 0, is equivalent to gnet = NULL.
#' @param alpha Significance level for Bonferroni multiple testing correction 
#' used for the independence tests (by default, alpha = 0.05), as suggested 
#' by Shipley (2000).
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen, saving execution time.
#' @param ... Currently ignored.
#'
#' @details Shipley's d-separation test (Shipley, 2000) evaluates if two 
#' variables (X, Y) in a DAG are conditionally independent for a given 
#' conditioning set Z, assuming that the underlying statistical model is 
#' recursive (i.e., acyclic) and can be represented as a graph in which 
#' nodes (i.e., variables) are connected through directed edges (i.e., 
#' causal effects), with no self-loops. The conditioning set Z is 
#' represented by the union of the parent sets of X and Y (i.e., the 
#' "basis set"). SEMbap algorithm makes an exhaustive search of all possible 
#' missing edges, corresponding to the number of degrees of freedom (df) 
#' of the DAG. A new bow-free covariance is added if there is a significant 
#' X-Y association, after Bonferroni multiple testing correction. In addition, 
#' if a reference network is given, new interactions will be retained 
#' only if present in the reference. 
#' SEMbap function may work on any kind of directed network, by converting 
#' mixed graphs into directed acyclic graphs. 
#' See \code{\link[SEMgraph]{graph2dag}} for conversion method details.
#'
#' @return A list of 3 igraph objects:
#' \itemize{
#' \item "bap", the output bow-free acyclic path diagram,
#' \item "guu", the bidirected graph of significant covariances,
#' \item "gLV", the directed graph of latent variables (LV) underlying
#' significant covariances (i.e., the canonical graph, where bidirected
#' X <-> Y edges are substituted by directed edges X <- LV -> Y).
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom dagitty lavaanToGraph impliedConditionalIndependencies localTests
#' @importFrom stats na.omit var qchisq pchisq p.adjust
#' @importFrom corpcor is.positive.definite cor.shrink
#' @importFrom flip flip plot
#' @export
#'
#' @references
#'
#' Shipley B (2000). A new inferential test for path models based on DAGs. 
#' Struct. Equ. Modeling, 7(2): 206-218. 
#' https://doi.org/10.1207/S15328007SEM0702_4
#'
#' Brito C and Pearl J (2002). A New Identification Condition for 
#' Recursive Models With Correlated Errors. 
#' Structural Equation Modeling, 9(4): 459-474.
#'
#' Nowzohour C, Maathuis MH, Evans RJ, Buhlmann P (2017). Distributional 
#' equivalence and structure learning for bow-free acyclic path diagrams. 
#' Electron. J. Stat., 11(2): 5342-5374.
#' https://doi.org/10.1214/17-EJS1372
#'
#' @examples
#' 
#' # BAP estimation
#' G <- SEMbap(graph = sachs$graph, data = log(sachs$pkc))
#' 
#' # Model fitting
#' sem <- SEMrun(graph = G$bap, data = log(sachs$pkc), group = sachs$group)
#' 
#' # Graphs
#' gplot(G$bap)
#' plot(G$guu)
#' plot(G$gLV)
#'
SEMbap <- function(graph, data, gnet = NULL, d = 0, alpha = 0.05,
                   verbose = FALSE, ...)
{
	# Set SEM objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	dataY <- data[,nodes]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	
	# Shipley's d-separation local tests & overall C-test
	dag <- Shipley.test(graph = ig, data = dataY, verbose = verbose)
	dsep <- dag$dsep
	d_sep <- subset(dsep, p.adjust(dsep$pvalue, method = "bonferroni") < alpha)
	if (nrow(d_sep) == 0) {
		cat("\n", "ALL adjusted pvalues >", alpha, "!", "\n\n")
		guu <- make_empty_graph(n = 0)
		return(list(bap = guu, guu = guu, gLV = guu))
	} else {
		cat("\nN.significant local tests:", nrow(d_sep), "of", nrow(dsep),
		    "\n\n")
	}
	guu <- graph_from_data_frame(d_sep[, 1:2], directed = FALSE)
	
	# Covariance and latent variables graphs
	if (d > 0) {
		guu <- EXT_SET(graph = guu, gnet = gnet, d = d, dag = FALSE)
	}
	if(ecount(guu) == 0) return(list(bap = guu, guu = guu, gLV = guu))
	V(guu)$color <- "white"
	ftm <- as_edgelist(as.undirected(guu))
	ftmLV <- NULL
	for (i in 1:nrow(ftm)) {
		ftmLV <- rbind(ftmLV, cbind(rep(paste0("L", i), 2), ftm[i,]))
	}
	gLV <- graph_from_data_frame(ftmLV, directed = TRUE)
	V(gLV)$color <- ifelse(substr(V(gLV)$name, 1, 1) == "L", "yellow", "white")
	if(verbose) {
		plot(guu, main = "extended covariance graph (guu)")
		Sys.sleep(3)
		plot(gLV, main = "extended latent variables graph (gLV)")
		Sys.sleep(0)
	}
	
	# SEM fitting with interactome bow-free covariances
	guu <- as.directed(guu, mode = "mutual")
	Ug <- graph.union(g = list(dag$dag, guu)[!is.na(list(ig, guu))])
	E1 <- attr(E(Ug), "vnames")
	E0 <- attr(E(dag$dag), "vnames")
	E(Ug)$color <- ifelse(E1 %in% E0, "blue", "red")
	if(!is.null(E(Ug)$weight)) Ug <- delete_edge_attr(Ug, "weight")
	
	return(list(bap = Ug, guu = as.undirected(guu), gLV = gLV))
}


#' @title Missing edge testing implied by a graph
#'
#' @description Compute all the P-values of the d-separation tests 
#' implied by the missing edges of a given graph. The results of every 
#' test is then combined using the Fisher’s statistic in an overall test 
#' (Shipley's test) of the fitted model C = -2*sum(log(P-value(k))), where 
#' C is distributed as a chi-squared variate with df = 2k, as suggested 
#' by Shipley (2000).
#'
#' @param graph A directed graph as an igraph object.
#' @param data A data matrix with subjects as rows and variables as 
#' columns.
#' @param verbose If TRUE, Shipley's test results will be showed to 
#' screen (default = TRUE).
#' @param ... Currently ignored.
#' 
#' @import igraph
#' @export
#'
#' @return A list of two objects: (i) the list of all d-separation tests 
#' over missing edges in the input DAG or BAP, and (ii) the DAG or BAP 
#' used to perform the Shipley test.
#'
#' @examples
#' 
#' library(huge)
#' als.npn <- huge.npn(alsData$exprs)
#' 
#' sem <- SEMrun(alsData$graph, als.npn)
#' C.test0 <- Shipley.test(sem$graph, als.npn)
#' 
#' adjData <- diagonalizePsi(alsData$graph, als.npn, algo = "d-sep",
#'                           method = "BH",
#'                           alpha = 0.05)
#' C.test1 <- Shipley.test(sem$graph, adjData)
#'
Shipley.test <- function(graph, data, verbose = TRUE, ...)
{
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	graph <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	data <- as.matrix(data[, nodes])
	
	# Graph to DAG conversion
	if (!is_dag(graph)) {
		cat("\nWARNING: input graph is not acyclic!\n")
		cat(" Applying graph -> DAG conversion...\n")
		dag <- graph2dag(graph, data) # delete <-> & cycles
		df1 <- vcount(graph)*(vcount(graph) - 1)/2 - ecount(graph)
		df2 <- vcount(dag)*(vcount(dag) - 1)/2 - ecount(dag)
		cat(" \nDegrees of freedom:\n  Input graph =",
            df1, "\n          DAG =", df2, "\n\n")
	} else {
		dag <- graph
	}
	
	# Shipley's d-separation local tests & overall C-test
	dsep <- dsep.test(dag = dag, S = cov(data), n = nrow(data))
	if (verbose) {
		ctest <- -2*sum(log(dsep$pvalue))
		df <- 2*nrow(dsep)
		pv <- 1 - pchisq(q = ctest, df = df)
		print(c(C_test = ctest, df = df, pvalue = pv))
	}
	
	return(list(dsep = dsep, dag = dag))
}

#' @title Estimate the optimal DAG from an input graph
#'
#' @description Extract the optimal DAG from an input graph, using the 
#' LASSO-based algorithm, implememted in \code{\link[glmnet]{glmnet}}.
#'
#' @param graph An igraph object.
#' @param data A matrix whith rows corresponding to subjects, and
#' columns to graph nodes (variables).
#' @param gnet Reference "global" network as an igraph object. If given, 
#' new edges will be added to the final DAG only if present in the 
#' reference network.
#' @param d An integer value indicating the maximum length of indirect
#' interactions between pairs of nodes. If d = 1, direct interactions 
#' between nodes will be searched in the reference interactome (if given). 
#' If d > 1, indirect interactions of length d or shorter (i.e., with at 
#' most d - 1 connectors) between bow-free nodes will be searched. 
#' Setting d = 0, is equivalent to gnet = NULL.
#' @param beta Numeric value. Minimum absolute LASSO beta coefficient for 
#' a new interaction to be retained in the final model. By default, beta 
#' is set to 0.
#' @param lambdas A vector of regularization LASSO lambda values. 
#' Cross-validation (n > 100) or BIC-based (n <= 100) optimal lambdas 
#' for each response variable will be selected. If lambdas is NULL, the 
#' \code{\link[glmnet]{glmnet}} default is enabled. If lambdas is NA 
#' (default), the tuning-free scheme is enabled by fixing 
#' lambdas = sqrt(log(p)/n), as suggested by Janková and van de Geer (2015). 
#' This will both reduce computational time and provide the same result 
#' at each run.
#' @param verbose A logical value. If FALSE (default), the processed graphs
#' will not be plotted to screen.
#' @param ... Currently ignored.
#'
#' @details The optimal DAG is estimated using successive penalized (L1) 
#' regressions. If the input graph is not acyclic, a warning message will 
#' be raised, and a cycle-breaking algorithm will be applied (see 
#' \code{\link[SEMgraph]{graph2dag}} for details). Output DAG edges will 
#' be colored in blue, if they were present in the input graph, and in 
#' red, if they are generated after topological (or data-driven) sorting.
#'
#' @return A list of 3 igraph objects:
#' \enumerate{
#' \item "dag", the estimated DAG;
#' \item "dag.red", new estimated connections;
#' \item "dag.blue", connections preserved from the input graph.
#' }
#'
#' @import igraph
#' @import lavaan
#' @importFrom glmnet glmnet
#' @importFrom RcppEigen fastLm
#' @export
#'
#' @seealso \code{\link[SEMgraph]{modelSearch}}
#'
#' @references
#'
#' Shojaie A, Michailidis G (2010). Penalized likelihood methods for 
#' estimation of sparse high-dimensional directed acyclic graphs. 
#' Biometrika, 97(3): 519-538. https://doi.org/10.1093/biomet/asq038
#'
#' Tibshirani R, Bien J, Friedman J, Hastie T, Simon N, Taylor J, 
#' Tibshirani RJ (2012). Strong rules for discarding predictors in 
#' lasso‐type problems. Royal Statistical Society: Series B 
#' (Statistical Methodology), 74(2): 245-266. 
#' https://doi.org/10.1111/j.1467-9868.2011.01004.x
#'
#' @examples
#' 
#' # DAG estimation
#' G <- SEMdag(graph = sachs$graph, data = log(sachs$pkc), beta = 0.05)
#' 
#' # Model fitting
#' sem <- SEMrun(graph = G$dag, data = log(sachs$pkc), group = sachs$group)
#' 
#' # Graphs
#' par(mfrow=c(2,2), mar=rep(1,4))
#' plot(sachs$graph, layout=layout.circle, main="input graph")
#' plot(G$dag, layout=layout.circle, main = "Output DAG")
#' plot(G$dag.blue, layout=layout.circle, main = "Old inferred edges")
#' plot(G$dag.red, layout=layout.circle, main = "New inferred edges")
#'
SEMdag <- function(graph, data, gnet = NULL, d = 0, beta = 0, lambdas = NA,
                   verbose = FALSE, LO = "topo", ...)
{
	# Set SEM objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	if (!is_dag(ig)) {
		cat("\nWARNING: the input graph is not acyclic !\n")
		cat(" Applying graph -> DAG conversion ...\n")
		dag <- graph2dag(ig, data)
	} else {
		dag <- ig
	}
	X <- scale(data[, V(dag)$name])
	
	# Estimate DAG using top-down approach
	cat("Linear Ordering with", LO, "setting\n")
	x <- DAG_TD(graph = dag, X = X, beta = beta, LO = LO, lambdas = lambdas)
	colnames(x$adj) <- rownames(x$adj) <- colnames(X)
	ig1 <- graph_from_adjacency_matrix(x$adj, mode = "directed")
	ig2 <- quiet(properties(ig1)[[1]])
	
	# Mapping DAG edges on reference interactome
	if (d > 0) ig2 <- EXT_SET(graph = ig2, gnet = gnet, d = d, dag = TRUE)
	E1 <- attr(E(ig2), "vnames")
	E0 <- attr(E(ig), "vnames")
	E(ig2)$color <- ifelse(E1 %in% E0, "blue", "red")
	if (verbose & ecount(ig2) != 0) gplot(ig2)
	ig3 <- ig2 - E(ig2)[which(E(ig2)$color == "blue")]
	ig3 <- ig3 - vertices(V(ig3)$name[igraph::degree(ig3) == 0])
	ig4 <- ig2 - E(ig2)[which(E(ig2)$color == "red")]
	ig4 <- ig4 - vertices(V(ig4)$name[igraph::degree(ig4) == 0])
	
	return(list(dag = ig2, dag.red = ig3, dag.blue = ig4))
}

DAG_TD <- function(graph, X, beta, LO, lambdas, ...)
{
	n <- dim(X)[1]
	p <- dim(X)[2]
	if (LO == "topo") TO <- names(igraph::topo_sort(graph))
	rr <- rev(match(TO, colnames(X)))
	result <- matrix(0, p, p)
	A<- as_adj(graph, sparse = FALSE)
	l <- sqrt(log(p)/n)
	for (ii in 1:(p - 1)) {
		now <- rr[ii]
		this <- sort(rr[(ii + 1):p])
		pw <- 1 - A[this, now]
		if (sum(pw) <= 0) next
		if (length(this) > 1) {
			if (is.na(lambdas)) {
				lassom <- glmnet::glmnet(X[, this], X[, now], alpha = 1,
				                         lambda = l,
			                             penalty.factor = pw)
				bfit <- coefficients(lassom)[-1]
			} else {
				if (n > 100) {
					lassom <- glmnet::cv.glmnet(X[, this], X[, now],
					                            lambda = lambdas,
					                            penalty.factor = pw)
					bfit <- coefficients(lassom)[-1]
				} else {
					lassom <- glmnet::glmnet(X[, this], X[, now],
					                         lambda = lambdas,
					                         penalty.factor = pw)
					bic <- n * log(colSums((predict(lassom, X[, this]) - 
					                        X[, now])^2)/n) + 
					                        lassom$df*log(n) + 
					                        2*lassom$df*log(p - ii)
					bfit <- coefficients(lassom)[, which(bic == min(bic))[1]][-1]
				}
			}
			for (jj in 1:length(this)) {
				if (abs(bfit[jj]) > beta)
				result[this[jj], now] <- 1
			}
		} else {
			lmod <- summary(RcppEigen::fastLm(X[, now] ~ X[, this]))
			if (lmod$coef[2, 4] < 0.05) {
				result[this, now] <- 1
			}
		}
	}
	return(list(adj = result, TO = rev(rr)))
}

EXT_SET <- function(graph, gnet, d = 2, dag = TRUE, ...)
{
	SET1 <- as_edgelist(graph)
	if( nrow(SET1) == 0 ) {
		cat("\nWARNING: n.interactions = 0 !\n\n")
		return(guu = NULL)
	}
	
	ftm1 <- NULL
	for(j in 1:nrow(SET1)) {
		cat("\r", "edge set", j, "of", nrow(SET1))
		flush.console()
		a <- SET1[j, 1]
		b <- SET1[j, 2]
		ftm1 <- rbind(ftm1, c(a, b))
		v <- which(V(gnet)$name %in% c(a, b))
		
		if (length(v) == 2) {
			if (dag == FALSE) {
				sp <- distances(gnet, a, b, mode = "all", weights = NA)
			}
			if (dag == TRUE) {
				sp <- distances(gnet, a, b, mode = "out", weights = NA)
			}
			if (sp <= d) {
				ftm1[j,] <- c(a, b)
			} else {
				ftm1[j,] <- c(NA, NA)
			}
		} else {
			ftm1[j,] <- c(NA, NA)
		}
	}
	
	ftm1 <- rbind(na.omit(ftm1), SET1[SET1[, 1] == "group",])
	cat("\n\n", "N.selected interactions:", nrow(SET1),
	    " N.imported from interactome:", nrow(ftm1), "\n\n")
	guu <- graph_from_edgelist(ftm1, direct = dag)
	
	return(guu)
}

#' @title Compute the Average Causal Effect (ACE) for a given source-sink pair
#'
#' @description Compute total effects as ACEs of source variables X 
#' (i.e., incoming connectivity = 0) on sink variables Y (i.e., outgoing 
#' connectivity = 0), in a directed graph. The ACE will be estimated as 
#' the path coefficient of X (i.e., theta) in the linear equation 
#' Y ~ X + Z. Z is defined as the adjustment (or conditioning) set of 
#' Y over X, applying an "optimal" valid set (O-set), with the smallest 
#' asymptotic variance. Standard errors (SE), for each ACE, are computed 
#' following the \code{lavaan} standard procedure or a bootstrap-based 
#' procedure (see \code{\link[boot]{boot}} for details).
#'
#' @param graph An igraph object.
#' @param data A matrix or data.frame. Rows correspond to subjects, and
#' columns to graph nodes (variables).
#' @param group A binary vector. This vector must be as long as the
#' number of subjects. Each vector element must be 1 for cases and 0
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param method Multiple testing correction method. One of the values 
#' available in \code{\link[stats]{p.adjust}}. By default, method is set 
#' to "none" (i.e., no multiple test correction).
#' @param alpha Significance level for ACE selection (by default, alpha = 0.05).
#' @param boot The number of bootstrap samplings enabling bootstrap 
#' computation of ACE standard errors. If NULL (default), the bootstrap 
#' is disabled.
#' @param ... Currently ignored.
#'
#' @return A data.frame of ACE estimates between network sources and sinks.
#'
#' @import igraph
#' @import lavaan
#' @importFrom dagitty adjustmentSets
#' @export
#'
#' @references
#'
#' Witte J, Henckel L, Maathuis MH, Didelez V (2020). On efficient 
#' adjustment in causal graphs. arXiv:2002.06825 [math.ST].
#' URL: https://arxiv.org/abs/2002.06825
#' 
#' @examples
#' 
#' # ACE estimation, without group (default)
#' ace <- SEMace(graph = sachs$graph, data = log(sachs$pkc))
#' print(ace)
#' 
#' # ACE estimation, with group perturbation and multiple test correction
#' ace2 <- SEMace(graph = sachs$graph, data = log(sachs$pkc),
#'                group = sachs$group,
#'                method = "BH", alpha = 0.05)
#' print(ace2)
#'
SEMace <- function(graph, data, group = NULL, method = "none", alpha = 0.05,
                   boot = NULL, ...)
{
	# Set igraph and dagitty graph objects
	nodes <- colnames(data)[colnames(data) %in% V(graph)$name]
	ig <- induced_subgraph(graph, vids = which(V(graph)$name %in% nodes))
	
	if (!is_dag(ig)) {
		cat("\nWARNING: input graph is not acyclic!\n")
		cat(" Applying graph -> DAG conversion.\n")
		dag <- graph2dag(ig, data) # delete cycles & all <->
		if (!is_dag(dag)) return(theta = NULL)
	} else {
		dag <- ig
	}
	
	dagy <- graph2dagitty(dag, verbose = FALSE)
	
	# Set distance matrix and distance graph from source to target nodes
	D <- igraph::distances(dag, mode = "out", weights = NA)
	D <- ifelse(D == Inf, 0, D)
	din <- igraph::degree(dag, mode = "in")
	Vx <- V(dag)$name[din == 0]
	dout <- igraph::degree(dag, mode = "out")
	Vy <- V(dag)$name[dout == 0]
	Dxy <- D[c(Vx, Vy), c(Vx, Vy)]
	gD <- simplify(graph_from_adjacency_matrix(Dxy, mode = "directed",
	               weighted = TRUE))
	cat("\nFrequency distribution of path length from X to Y :")
	print(table(E(gD)$weight))
	cat("\n")
	
	# Compute total effect (ACE = theta) from DAG(CPDAG), MAG(PAG)
	
	theta <- NULL
	res <- NULL
	ftm <- as_edgelist(gD)
	
	for (i in 1:nrow(ftm)) {
		cat("\r", "ACE", i, "of", nrow(ftm))
		flush.console()
		x <- ftm[i, 1]
		y <- ftm[i, 2]
		
		# Adjustement SET Z using adjustmentSets() of dagitty package
		#Z <- dagitty::adjustmentSets(dag, x, y, type = "minimal",
		#                             effect = "total")
		#if (length(Z) == 0) next
		#X <- scale(data[, c(x, unlist(Z[[1]]))])
		#Y <- scale(data[, y])
		
		# Adjustement SET Z using OptAdjSet (Witte et al, 2020)
		#paths <- all_simple_paths(dag, from = x, to = y, mode = "out")
		#cn <- setdiff(unique(names(unlist(paths))), x)
		paths <- dagitty::paths(dagy, from = x, to = y, directed = TRUE)$paths
		cn <- setdiff(unique(unlist(strsplit(gsub("->", "", paths), "  "))), x)
		#pa_cn <- dagitty::parents(dag, cn)
		#forb <- c(dagitty::descendants(dag, cn),x)
		pa_cn <- V(dag)$name[parents(dag, cn)]
		forb <- c(V(dag)$name[descendants(dag, cn)], x)
		z <- setdiff(pa_cn, forb)
		Z <- scale(data[, c(y, x, z)])
		
		# LM fitting
		if (is.null(group)) {
			if (is.null(boot)) {
				try(est <- lmest(x, y, Z))
			} else {
				try(est <- boot.lmest(x, y, Z, R = boot))
			}
			try(res <- data.frame(lapply(est,
			    function(y) if(is.numeric(y)) round(y, 3) else y)))
		} else {
			try(est <- lmest2(x, y, Z, group, boot = boot))
			try(res <- data.frame(lapply(est,
			    function(y) if(is.numeric(y)) round(y, 3) else y)))
		}
		theta <- rbind(theta, res)
	}
	
	cat("\n")
	theta <- subset(theta, p.adjust(theta$pvalue, method = method) < alpha)
	
	return(theta)
}

lmest <- function(x, y, Z, ...)
{
	# LM fitting y ~ x + Z
	
	fit <- stats::lm.fit(as.matrix(Z[, -1]), Z[, 1])
	est <- as.numeric(fit$coefficients)[1]
	sigma <- sum(fit$residuals^2)/fit$df.residual
	
	if (is.na(sigma)) return(NULL)
	X <- as.matrix(Z[, -1])
	r <- fit$rank
	p <- ncol(X)
	
	if (r == p) {
		se <- sqrt(sigma*diag(solve(t(X)%*%X)))[1]
	} else {
		C <- t(X)%*%X
		E <- eigen(C)
		W <- E$vectors[1:p, 1:r]%*%diag(1/E$values[1:r])%*%t(E$vectors[1:p, 1:r])
		se <- sqrt(sigma*diag(W))[2]
	}
	
	z <- est/se
	res <- data.frame(sink = y, op = "<-", source = x, est = est, se = se,
	                  z = z,
	                  pvalue = 2*(1 - pnorm(abs(z))),
	                  ci.lower = (est - 1.96*se),
	                  ci.upper = (est + 1.96*se))
	return(res)
}

boot.lmest <- function(x, y, Z, R, ...)
{
	# LM fitting y ~ x + Z
	est <- function(Z, i) {
		stats::lm.fit(as.matrix(Z[i, -1]), Z[i, 1])$coefficients[1]
	}
	xboot <- boot::boot(Z, est, R = R)
	t0 <- xboot$t0
	se <- sd(xboot$t)
	z <- t0/se
	res <- data.frame(sink = y, op = "<-", source = x, est = t0, se = se,
	                  z = z,
	                  pvalue = 2*(1 - pnorm(abs(z))),
	                  ci.lower = (t0 - 1.96*se),
	                  ci.upper = (t0 + 1.96*se))
	return(res)
}

lmest2 <- function(x, y, Z, group, boot, ...)
{
	# LM fitting y ~ x + Z  w/n groups
	Z1 <- as.matrix(Z[group == 1,])
	Z0 <- as.matrix(Z[group == 0,])
	if (is.null(boot)) {
		est1 <- lmest(x, y, Z1)
		est0 <- lmest(x, y, Z0)
	} else {
		est1 <- boot.lmest(x, y, Z1, R = boot)
		est0 <- boot.lmest(x, y, Z0, R = boot)
	}
	d_est <- est1$est - est0$est
	d_se <- sqrt(est1$se^2 + est0$se^2)
	d_z <- d_est/d_se
	pvalue <- 2*(1 - pnorm(abs(d_z)))
	d_lower <- d_est - 1.96*d_se
	d_upper <- d_est + 1.96*d_se
	res <- cbind(est1[, 1:3], d_est, d_se, d_z, pvalue, d_lower, d_upper)
	return(res)
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
#' for control subjects. If NULL (default), group influence will not be
#' considered.
#' @param from Starting node name (i.e., source node).
#' @param to Ending node name (i.e., sink node).
#' @param path If path = "directed", all directed paths between the two 
#' nodes will be included in the fitted model. If path = "shortest", only 
#' shortest paths will be returned.
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
#' @import igraph
#' @import lavaan
#' @importFrom dagitty paths
#' @export
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
		cat("ValueError: infinite distance from", from, "to", to, ".\n\n")
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
		dag <- SEMgraph:::graph2dagitty(ig, verbose = FALSE)
		paths <- dagitty::paths(dag, from, to, directed = TRUE)$paths
		nodes <- unique(unlist(strsplit(gsub("->", "", paths), "  ")))
	}
	
	# Plot selected paths
	ig1 <- induced_subgraph(ig, vids = which(V(ig)$name %in% nodes))
	V(ig)$color <- "white"
	V(ig)$color[V(ig)$name %in% V(ig1)$name] <- "gold"
	E(ig)$color <- "gray40"
	E(ig)$color[attr(E(ig), "vnames") %in% attr(E(ig1), "vnames")] <- "darkorange3"
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

dsep.test <- function(dag, S, n, ...)
{
 	# d-sep (basis set) testing of a DAG
	ndf <- vcount(dag)*(vcount(dag) - 1)/2 - ecount(dag)
	idx <- as.numeric(topo_sort(dag, mode = "out"))
	A <- as_adj(dag, sparse = FALSE)[idx, idx]
	nodes <- colnames(A)
	SET <- NULL
	i <- 1
	for (r in 1:(length(nodes) - 1)) {
		for (s in (r + 1):length(nodes)) {
			cat("\r", "Basis set", i, "of", ndf)
			flush.console()
			if (A[r, s] != 0) next
			ed <- nodes[c(r, s)]
			pa.r <- nodes[A[, r] == 1]
			pa.s <- nodes[A[, s] == 1]
			dsep <- union(pa.r, pa.s)
			dsep <- setdiff(dsep, ed)
			B <- c(ed, dsep)
			if (length(B) > (n - 3)) next
			pvalue <- pcor.test(S, B, n, H0 = 0.05)
			set <- paste(B[-c(1:2)], collapse = ",")
			SET <- rbind(SET, data.frame(X = B[1], Y = B[2], SET = set,
			             pvalue))
			i <- i + 1
		}
	}
	cat("\n\n")
	return(SET)
}

pcor.test <- function(S, B, n, H0 = 0, ...)
{
	k <- solve(S[B, B])
	r <- -k[1, 2]/sqrt(k[1, 1]*k[2, 2])
	q <- length(B) - 2
	if (H0 == 0) {
		df <- n - 2 - q
		tval <- r*sqrt(df)/sqrt(1 - r^2)
		pval <- 2*pt(-abs(tval), df)
	} else {
		z <- atanh(r)
		se <- 1/sqrt(n - 3 - q)
		pval <- pchisq((z/se)^2, df = 1, ncp = (atanh(.05)/se)^2,
		               lower.tail = FALSE)
	}
	return(pval)
}

pb <- function(n_rep)
{
	pb <- progress::progress_bar$new(
		#format = "  |:bar| :percent elapsed = :elapsed  ~ :eta",
		format = "  got :current of :total rows in :elapsed",
		total = n_rep, width = 80, clear = FALSE, show_after = 0)
	progress <- function(x) pb$tick(tokens = list(trial = (1:n_rep)[x]))
 	return(progress)
}

quiet <- function(x) {
	sink(tempfile())
	on.exit(sink())
	invisible(force(x))
}
