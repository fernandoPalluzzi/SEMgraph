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
\description{
Generate a summary for a RICF model and show it to 
standard output.
}
\examples{
sem1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "ricf")
summary(sem1$fit)

}
\seealso{
\code{\link[SEMgraph]{SEMrun}}.
}
\author{
Mario Grassi \email{mario.grassi@unipv.it}
}