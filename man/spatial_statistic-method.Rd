\name{spatial_statistic-methods}
\docType{methods}
\alias{spatial_statistic-methods}
\alias{spatial_statistic,infpan-method}
\title{Methods for Function \code{spatial_statistic}}
\description{
 Methods for function \code{spatial_statistic}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan",
    attribute = "Cases",
    timepoint = NULL,
    statistic = "nbstat",
    func = NULL,
    row.names = NULL,
    snap = NULL,
    queen = TRUE,
    weights_style = "W",
    weights_zeropolicy = NULL,
    verbose = FALSE,
    ...)}}{
Calculation of a specified spatial statistic for the data from an infpan object.
The user must choose a specific time point using the \code{timepoint} argument. 
If no time point is stated, the most current time point is chosen automatically.
}
}}
\author{
Thomas Wieland
}
\keyword{methods}