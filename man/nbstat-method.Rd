\name{nbstat-methods}
\docType{methods}
\alias{nbstat-methods}
\alias{nbstat,infpan-method}
\title{Methods for Function \code{nbstat}}

\description{
 Methods for function \code{nbstat}
}

\section{Methods}{
\describe{

\item{\code{signature(
    object = "nbmatrix",
    link_data,
    ID_col,
    data_col,
    func = "sum",
    verbose = FALSE)}}{
Calculates a specified statistic (e.g., mean or sum) for the neighboring regions within a \code{nbmatrix} object.
The user must specify a (spatial) data frame with the parameter \code{link_data}.
The results are saved in the \code{nbmatrix} object.
If \code{func} is set to NULL, \code{nbstat} is used to specify the data column used for further analyses (e.g., \code{moran()}).
}
}}

\author{
Thomas Wieland
}

\keyword{methods}