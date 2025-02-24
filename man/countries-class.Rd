\name{countries-class}
\docType{class}
\alias{countries-class}
\alias{plot,countries-method}
\alias{show,countries-method}
\alias{summary,countries-method}

\title{Class \code{"countries"}}
\description{
The class \code{"countries"} contains the results of a two-country comparison analysis using the Swash-Backwash Model, including two \code{"sbm_ci"} classes for each country. Use \code{summary(countries)} and \code{plot(countries)} for results summary and plotting, respectively.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("countries", ...)}.
Objects can be created by the function \code{compare_countries(sbm1, sbm2)}. 
}
\section{Slots}{
  \describe{
    \item{\code{sbm_ci1}:}{Object of class \code{"sbm_ci"} Results of \code{"confint(sbm1)"} for country 1 }
    \item{\code{sbm_ci2}:}{Object of class \code{"sbm_ci"} Results of \code{"confint(sbm1)"} for country 2 }
    \item{\code{D}:}{Object of class \code{"numeric"} Results: Difference \eqn{D} between the samples with respect to the chosen indicator }
    \item{\code{D_ci}:}{Object of class \code{"numeric"} Results: \eqn{\alpha} confidence intervals of \eqn{D}}
    \item{\code{config}:}{Object of class \code{"list"} Configuration details for bootstrap sampling }
    \item{\code{country_names}:}{Object of class \code{"character"} User-stated country names }
    \item{\code{indicator}:}{Object of class \code{"character"} User-stated indicator to be tested }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "countries")}: Plots the results of a two-country comparison with the Swash-Backwash Model }
    \item{show}{\code{signature(object = "countries")}: Prints an \code{countries} object; use \code{summary(sbm_ci)} for results }
    \item{summary}{\code{signature(object = "countries")}: Prints a summary of a \code{countries} object (results of the two-country comparison) }
	 }
}
\author{
Thomas Wieland
}
\examples{
showClass("countries")
}
\keyword{classes}