\name{moran-methods}
\docType{methods}
\alias{moran-methods}
\alias{moran,nbmatrix-method}
\title{Methods for Function \code{moran}}
\description{
 Methods for function \code{moran}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "nbmatrix",
    randomization = TRUE,
    alternative = "greater",
    verbose = FALSE)}}{
Calculates the Global Moran's I statistic for the neighborhood matrix
for a specified data column to be analyzed.
The parameters \code{randomization} and \code{alternative} are passed to \code{spdep::moran.test}}, 
which is internally used for the calculation.
}
}

\author{
Thomas Wieland
}

\keyword{methods}