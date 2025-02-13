\name{sbm-class}
\docType{class}
\alias{sbm-class}
\alias{confint,sbm-method}
\alias{plot,sbm-method}
\alias{print,sbm-method}
\alias{show,sbm-method}
\alias{summary,sbm-method}

\title{Class \code{"sbm"}}
\description{
The class \code{"sbm"} contains the results of the Swash-Backwash Model and the related input data as well as additional information. Use \code{summary(sbm)} and \code{plot(sbm)} for results summary and plotting, respectively.
}
\section{Objects from the Class}{
Objects can be created by the function \code{swash}.
}
\section{Slots}{
  \describe{
    \item{\code{R_0A}:}{Object of class \code{"numeric"} Model result: spatial reproduction number \eqn{R_{0A}}}
    \item{\code{integrals}:}{Object of class \code{"numeric"} Model result: integrals \eqn{S_A}, \eqn{I_A}, and \eqn{R_A} }
    \item{\code{velocity}:}{Object of class \code{"numeric"} Model result: velocity measures \eqn{t_{FE}} and \eqn{t_{LE}} }
    \item{\code{occ_regions}:}{Object of class \code{"data.frame"} Model result: Occurence at regional level }
    \item{\code{cases_by_date}:}{Object of class \code{"data.frame"} Total cases by date }
    \item{\code{input_data}:}{Object of class \code{"data.frame"} Input data }
    \item{\code{data_statistics}:}{Object of class \code{"numeric"} Diagnostics of input data }
    \item{\code{col_names}:}{Object of class \code{"character"} Column names in input data }
  }
}
\section{Methods}{
  \describe{
    \item{confint}{\code{signature(object = "sbm")}: Creates bootstrap confidence intervals for \code{sbm} objects. }
    \item{plot}{\code{signature(x = "sbm")}: Plots the results of the Swash-Backwash Model; two plots: edges over time, total infections per time unit }
    \item{print}{\code{signature(x = "sbm")}: Prints an \code{sbm} object; use \code{summary(sbm)} for results }
    \item{show}{\code{signature(object = "sbm")}: Prints an \code{sbm} object; use \code{summary(sbm)} for results }
    \item{summary}{\code{signature(object = "sbm")}: Prints a summary of \code{sbm} objects (results of the Swash-Backwash Model) }
	 }
}
\author{
Thomas Wieland
}
\examples{
showClass("sbm")
}
\keyword{classes}