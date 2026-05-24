\name{calculate_rollsum-methods}
\docType{methods}
\alias{calculate_rollsum-methods}
\alias{calculate_rollsum,infpan-method}
\title{Methods for Function \code{calculate_rollsum}}

\description{
 Methods for function \code{calculate_rollsum}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan",
    k = 7,
    align = "center",
    fill = NA,
    col_name = NULL,
    verbose = FALSE)}}{
Calculates the rolling sum of the infections panel data for all time points for each region.
Set the rolling window by the parameter \code{k} (default: 7). 
Set the fill value for the observations left/within/right to the data range with parameter \code{fill} (default: NA).
Parameter \code{align} defines whether the index of the result should be left- or right-aligned or centered (default).
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_rs".
The method uses the function \code{rollsum} from the \code{zoo} package (see the corresponding documentation).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}