\name{calculate_cum-methods}
\docType{methods}
\alias{calculate_cum-methods}
\alias{calculate_cum,infpan-method}
\title{Methods for Function \code{calculate_cum}}

\description{
 Methods for function \code{calculate_cum}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan",
  col_name = NULL,
  overwrite = FALSE,
  verbose = FALSE)}}{
Calculates the cumulative values of the infections panel data for all time points for each region.
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_cum".
The method uses the function \code{cumsum} from the \code{base} package (see the corresponding documentation).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}