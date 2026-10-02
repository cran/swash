\name{plot_map-methods}
\docType{methods}
\alias{plot_map-methods}
\alias{plot_map,infpan-method}
\title{Methods for Function \code{plot_map}}
\description{
 Methods for function \code{plot_map}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan",
    attribute = "Cases",
    timepoint = NULL,
    verbose = FALSE,
    ...)}}{
Plots a map from an \code{infpan} object which displays the attribute specified with the \code{attribute} parameter
at the timepoint specified with \code{timepoint}.
To adjust the layout of the map, use the \code{...} parameters which are passed to the 
\code{plot} S3 method for the \code{sf} class (see the corresponding documentation).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}