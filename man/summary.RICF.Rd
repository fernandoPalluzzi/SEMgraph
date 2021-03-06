% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semPaths.R
\name{summary.RICF}
\alias{summary.RICF}
\title{RICF model summary}
\usage{
\method{summary}{RICF}(object, ...)
}
\arguments{
\item{object}{A RICF fitted model object.}

\item{...}{Currently ignored.}
}
\value{
Shown the lavaan-formatted summary to console
}
\description{
Generate a summary for a RICF solver similar to
lavaan-formatted summary
}
\examples{
sem0 <- SEMrun(sachs$graph, log(sachs$pkc), algo = "ricf")
summary(sem0$fit)

}
\seealso{
\code{\link[SEMgraph]{SEMrun}}.
}
\author{
Mario Grassi \email{mario.grassi@unipv.it}
}
