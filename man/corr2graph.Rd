% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semUtils.R
\name{corr2graph}
\alias{corr2graph}
\title{Correlation matrix to graph}
\usage{
corr2graph(R, n, type = "marg", method = "none", alpha = 0.05, ...)
}
\arguments{
\item{R}{Correlation matrix.}

\item{n}{Sample size (i.e., the number of subjects).}

\item{type}{Graph building method. If \code{type} is either
\code{"marg"} or \code{"cond"}, marginal or conditional correlation
tests will be used, respectively.
If \code{type = "mst"}, input correlations are converted to distances
and a minimum spanning tree is generated from the distance matrix,
using Prim's algorithm (Prim, 1957).
If \code{type = "tmfg"}, a triangulate maximally graph is generated
from the given correlation matrix (Massara et al., 2016).}

\item{method}{Multiple testing correction method. One of the values
available in \code{\link[stats]{p.adjust}}. By default,
\code{method = "none"} (i.e., no multiple test correction).
See \code{\link[stats]{p.adjust}} for other correction methods.}

\item{alpha}{Significance level used to compute the correlation threshold.
By default, \code{alpha = 0.05}.}

\item{...}{Currently ignored.}
}
\value{
An igraph object.
}
\description{
Convert a correlation matrix to an igraph object.
}
\examples{

# Graphs creation
C1 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
                 type = "marg",
                 method = "BH")
C2 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
                 type = "cond",
                 method = "BH")
C3 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
                 type = "mst",
                 method = "BH")
C4 <- corr2graph(R = cor(log(sachs$pkc)), n = nrow(sachs$pkc),
                 type = "tmfg",
                 method = "BH")

# Graphs plots
old.par <- par(no.readonly = TRUE)
par(mfrow=c(2,2), mar= rep(2, 4))
plot(C1, layout=layout.circle, main= "marg"); box(col="gray")
plot(C2, layout=layout.circle, main= "cond"); box(col="gray")
plot(C3, layout=layout.circle, main= "mst"); box(col="gray")
plot(C4, layout=layout.circle, main= "tmfg"); box(col="gray")
par(old.par)

}
\references{
Palluzzi F, Grassi M (2021). SEMgraph: An R Package for Causal Network
Analysis of High-Throughput Data with Structural Equation Models.
<arXiv:2103.08332>

Massara GP, Di Matteo T and Aste T (2009). Network Filtering for Big
Data: Triangulated Maximally Filtered Graph.
Journal of complex Networks, 5(2): 161--178.
<https://doi.org/10.1093/comnet/cnw015>

Prim RC (1957). Shortest connection networks and some generalizations.
Bell System Technical Journal, 36(6):1389--1401.
<https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
}
\author{
Mario Grassi \email{mario.grassi@unipv.it}
}
