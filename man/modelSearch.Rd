% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semUtils.R
\name{modelSearch}
\alias{modelSearch}
\title{Optimal model search strategies}
\usage{
modelSearch(
  graph,
  data,
  gnet = NULL,
  d = 2,
  search = "basic",
  beta = 0,
  alpha = 0.05,
  pstop = TRUE,
  limit = 30000,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{graph}{Input network as an igraph object.}

\item{data}{A matrix or data.frame. Rows correspond to subjects, and
columns to graph nodes (variables).}

\item{gnet}{Reference directed network used to validate and import
nodes and interactions.}

\item{d}{Maximum allowed geodesic distance for directed or undirected
shortest path search. A distance \code{d = 0} disables shortest path
search (fixed in search = "basic"), while \code{d = 1} (fixed in
\code{search = "direct"}) only search for directed links (i.e., no
mediators are allowed).
A distance d > 1 (defaults to \code{d = 2} for "outer" and "inner"
strategies), will search for shortest paths with at most d - 1 mediators
between nodes sharing a significant estimated interaction.
Connectors are imported from the reference interactome, as specified
by the argument gnet. If the edges of the reference interactome are
weighted by P-value, as defined by the E(graph)$pv attribute, the
shortest path with the smallest sum of weights will be chosen (e.g.,
see \code{\link[SEMgraph]{weightGraph}} for graph weighting options).}

\item{search}{Search strategy. Four model search strategies are available:
\itemize{
\item "outer". The estimated DAG is extended using
\code{\link[SEMgraph]{extendGraph}} to find new indirect paths (i.e.,
inferred directed connections that may hide new mediators). New
interactions and mediators will be searched and imported from the
reference network (argument gnet, see above). Both DAG and extended
graph complexity can be controlled with beta > 0 and d > 1 arguments,
respectively (see below). The term "outer" means that new model mediator
variables are imported from an external resource (i.e., the reference
network).
\item "inner". This strategy is analogous to the "outer" one,
but disables external mediator search. In other words, new indirect
paths are generated by adding new interactions of the input model, so
that mediators will be nodes already present in the input graph. The
reference network is still used to validate new model paths. Also in
this case, beta > 0 and d > 1 are used.
\item "direct". The input graph structure is improved through direct
(i.e., adjacent) link search, followed by interaction validation and
import from the reference network, with no mediators
(i.e., \code{d = 1}).
\item "basic" (default). While the previous strategies rely on the
input graph and the reference network to integrate knowledge to the
final model, the "basic" strategy is data-driven. The input graph is
needed to define the topological order. The argument gnet is set
to NULL (i.e., no reference network is needed) and argument \code{d = 0}.
Model complexity can be still controlled by setting beta > 0.
}}

\item{beta}{Numeric value. Minimum absolute LASSO beta coefficient for
a new interaction to be retained in the estimated DAG backbone. Lower
beta values correspond to more complex DAGs. By default, beta
is set to 0 (i.e., maximum complexity).}

\item{alpha}{Significance level for false discovery rate (FDR) used
for either local d-separation tests (below \code{limit}) or conditional
independence (CI) test (above \code{limit}). This argument is used to
control data de-correlation. A higher \code{alpha} level includes more
hidden covariances, thus considering more sources of confounding.
If \code{alpha = 0}, data de-correlation is disabled.
By default, \code{alpha = 0.05}.}

\item{pstop}{A logical value. With the argument \code{pstop = TRUE}
(default), the algorithm can be halted when the Shipley's global model
test P-value > 0.05. If \code{pstop = FALSE}, the model search
algorithm stops when no additional edges can be added to the estimated
DAG.}

\item{limit}{An integer value corresponding to the number of missing
edges of the extracted acyclic graph. Beyond this limit, multicore
computation is enabled to reduce the computational burden.}

\item{verbose}{If TRUE, it shows intermediate graphs during the
execution (not recommended for large graphs).}

\item{...}{Currently ignored.}
}
\value{
The output model as well as the adjusted dataset are returned
as a list of 3 objects:
\itemize{
\item "fit", the fitted output model (lavaan object);
\item "graph", the output model as an igraph object;
\item "data", the adjusted dataset.
}
}
\description{
Four model search strategies are implemented combining
\code{SEMdag()}, \code{SEMbap()}, and \code{extendGraph()} functions.
All strategies estimate a DAG through the adjusted (de-correlate)
data matrix Z by iteratively update DAG and Z.
}
\details{
Search strategies can be ordered by decreasing conservativeness
respect to the input graph, as: "direct", "inner", "outer", and "basic".
The first three strategies are knowledge-based, since they require an
input graph and a reference network, together with data, for
knowledge-assisted model improvement. The last one does not require
any reference and the output model structure will be completely
determined by data.
Output model complexity can be limited using arguments d and beta.
While d is fixed to 0 or 1 in \code{"basic"} or \code{"direct"},
respectively; we suggest starting with \code{d = 2} (only one mediator)
for the other two strategies.
For knowledge-based strategies, we suggest to to start with
\code{beta = 0.1}. Then, beta can be relaxed (0 to < 0.1) to improve
model fitting, if needed. Since data-driven models can be complex,
we suggest to start from beta = 0.1 when using the "basic" strategy.
The beta value can be relaxed until a good model fit is obtained.
Argument alpha determines the extent of data adjustment: lower alpha
values for FDR correction correspond to a smaller number of significant
confounding factors, hence a weaker correction (default alpha = 0.05).
}
\examples{

\donttest{

# Comparison among different model estimation strategies

library(huge)
als.npn <- huge.npn(alsData$exprs)

# Models estimation
m1 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
      search = "direct", beta = 0, alpha = 0.05)
m2 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
      d = 2, search = "inner", beta = 0.05, alpha = 0.05)
m3 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = kegg,
      d = 2, search = "outer", beta = 0.05, alpha = 0.05)
m4 <- modelSearch(graph = alsData$graph, data = als.npn, gnet = NULL,
      search = "basic", beta = 0.1, alpha = 0.05)

# Graphs
#old.par <- par(no.readonly = TRUE)
#par(mfrow=c(2,2), mar= rep(1,4))
gplot(m1$graph, main = "direct graph")
gplot(m2$graph, main = "inner graph")
gplot(m3$graph, main = "outer graph")
gplot(m4$graph, main = "basic graph")
#par(old.par)

}

}
\author{
Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
}
