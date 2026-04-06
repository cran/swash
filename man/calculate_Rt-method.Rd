\name{calculate_Rt-methods}
\docType{methods}
\alias{calculate_Rt-methods}
\alias{calculate_Rt,infpan-method}
\title{Methods for Function \code{calculate_Rt}}
\description{
 Methods for function \code{calculate_Rt}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object = "infpan",
    GP = 4,
    correction = FALSE,
    col_name = NULL,
    verbose = FALSE)}}{
Calculates the effective reproduction number \eqn{R_t} for all time points for each region in the infections panel data.
Set the generation period by the parameter \code{GP} (default: 4). 
If \code{correction} is TRUE, values equal to zero are increased by one.
Set \code{overwrite} to TRUE, if an existing column should be overwritten.
The method uses the built-in function \code{R_t()}.
}
}}
\author{
Thomas Wieland
}
\keyword{methods}