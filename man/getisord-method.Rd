\name{getisord-methods}
\docType{methods}
\alias{getisord-methods}
\alias{getisord,nbmatrix-method}
\title{Methods for Function \code{getisord}}
\description{
 Methods for function \code{getisord}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "nbmatrix",
    alternative = "greater",
    verbose = FALSE)}}{
Calculates the Global Getis-Ord statistic for the neighborhood matrix
for a specified data column to be analyzed.
The parameters \code{randomization} and \code{alternative} are passed to \code{spdep::globalG.test}}, 
which is internally used for the calculation.
}
}

\author{
Thomas Wieland
}

\keyword{methods}