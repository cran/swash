\name{swash-methods}
\docType{methods}
\alias{swash-methods}
\alias{swash,sbm-method}
\title{Methods for Function \code{swash}}
\description{
 Methods for function \code{swash}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan", 
    verbose = FALSE)}}{
Performs the analysis of the input panel data with \eqn{N} regions and \eqn{T} time points
using the Swash-Backwash Model based on an \code{infpan} object.
The output is an object of class \code{"sbm"}. 
The results can be viewed using \code{summary(sbm)}.
See \code{\link{swash_backwash}} for further details of the model analysis.
}
}}

\author{
Thomas Wieland
}

\keyword{methods}