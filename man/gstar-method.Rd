\name{gstar-methods}
\docType{methods}
\alias{gstar-methods}
\alias{gstar,nbmatrix-method}
\title{Methods for Function \code{gstar}}

\description{
 Methods for function \code{gstar}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "nbmatrix",
    alternative = "greater",
    nsim = 999,
    verbose = FALSE)}}{
Calculates the Local Getis-Ord Gi* statistic for the neighborhood matrix
for a specified data column to be analyzed.
The parameters \code{alternative} and \code{nsim} are passed to \code{sfdep::local_gstar_perm}}, 
which is internally used for the calculation.
}
}

\author{
Thomas Wieland
}

\keyword{methods}