\name{growth-methods}
\docType{methods}
\alias{growth-methods}
\alias{growth,sbm-method}
\title{Methods for Function \code{growth}}
\description{
 Methods for function \code{growth}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "sbm", S_iterations = 10, S_start_est_method = "bisect", seq_by = 10, nls = TRUE, verbose = FALSE)}}{
Estimation of \eqn{N} logistic growth models for \eqn{N} regions. Both OLS and NLS estimation are estimated by default (set \code{nls} = \code{FALSE} to skip NLS estimation). Parameters \code{S_iterations}, \code{S_start_est_method}, and \code{seq_by} are used to control the saturation value estimation (see \code{\link{logistic_growth}}).
}
}}
\author{
Thomas Wieland
}
\keyword{methods}