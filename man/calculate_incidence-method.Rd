\name{calculate_incidence-methods}
\docType{methods}
\alias{calculate_incidence-methods}
\alias{calculate_incidence,infpan-method}
\title{Methods for Function \code{calculate_incidence}}
\description{
 Methods for function \code{calculate_incidence}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object,
  use_column = "Cases",
  col_name = NULL,
  pop_factor = 100000,
  overwrite = FALSE,
  verbose = FALSE)}}{
Calculates the incidence of the infections panel data for all time points for each region.
Use \code{use_column} to specify which column should be used for the calculation of incidence. 
The following values are permitted: "Cases" (default, incremental cases), "Cum. cases" (cumulative cases), "Roll. mean" (rolling mean of cases), 
or "Roll. sum" (rolling sum of cases). If the specified column does not exist in the infections panel data of the \code{infpan} object,
the function raises an error. If in the \code{infpan} object, no "Population" column is defined, incidence calculation is not possible.
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_inc".
Set \code{overwrite} to TRUE, if an existing column should be overwritten.
}
}}
\author{
Thomas Wieland
}
\keyword{methods}