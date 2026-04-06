\name{growth_initial-methods}
\docType{methods}
\alias{growth_initial-methods}
\alias{growth_initial,infpan-method}
\title{Methods for Function \code{growth_initial}}
\description{
 Methods for function \code{growth_initial}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "infpan", time_units = 10, GI = 4, nls = TRUE, nls_start = list(a = 1, b = 0.1), add_constant = 1, verbose = FALSE)}}{
Estimation of \eqn{N} exponential growth models for the initial phase of an epidemic spread for \eqn{N} regions. 
Set argument \code{GI} for the calculation of the basic reproduction number, 
and control OLS/NLS estimation with arguments \code{nls}, \code{nls_start}, and \code{add_constant} (see \code{\link{exponential_growth}}).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}