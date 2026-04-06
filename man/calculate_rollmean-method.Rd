\name{calculate_rollmean-methods}
\docType{methods}
\alias{calculate_rollmean-methods}
\alias{calculate_rollmean,infpan-method}
\title{Methods for Function \code{calculate_rollmean}}

\description{
 Methods for function \code{calculate_rollmean}
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
Calculates the rolling mean of the infections panel data for all time points for each region.
Set the rolling window by the parameter \code{k} (default: 7). 
Set the fill value for the observations left/within/right to the data range with parameter \code{fill} (default: NA).
Parameter \code{align} defines whether the index of the result should be left- or right-aligned or centered (default).
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_rm".
The method uses the function \code{rollmean} from the \code{zoo} package (see the corresponding documentation).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}