\encoding{UTF-8}
\name{sbm_ci-class}
\docType{class}
\alias{sbm_ci-class}
\alias{plot,sbm_ci-method}
\alias{print,sbm_ci-method}
\alias{show,sbm_ci-method}
\alias{summary,sbm_ci-method}

\title{Class \code{"sbm_ci"}}
\description{
The class \code{"sbm_ci"} contains the results of the Swash-Backwash Model, confidence intervals for the model estimates, and the related input data as well as additional information. Use \code{summary(sbm_ci)} and \code{plot(sbm_ci)} for results summary and plotting, respectively. 
}
\section{Objects from the Class}{
Objects can be created by the function \code{confint(sbm)}.
}
\section{Slots}{
  \describe{
    \item{\code{R_0A}:}{Object of class \code{"numeric"} Model result: spatial reproduction number \eqn{R_{0A}}}
    \item{\code{integrals}:}{Object of class \code{"numeric"} Model result: integrals \eqn{S_A}, \eqn{I_A}, and \eqn{R_A} }
    \item{\code{velocity}:}{Object of class \code{"numeric"} Model result: velocity measures \eqn{t_{FE}} and \eqn{t_{LE}} }
    \item{\code{occ_regions}:}{Object of class \code{"data.frame"} Model result: Occurence at regional level }
    \item{\code{cases_by_date}:}{Object of class \code{"data.frame"} Total cases by date }
    \item{\code{cases_by_region}:}{Object of class \code{"data.frame"} Cumulative cases by region }
    \item{\code{input_data}:}{Object of class \code{"data.frame"} Input data }
    \item{\code{data_statistics}:}{Object of class \code{"numeric"} Diagnostics of input data }
    \item{\code{col_names}:}{Object of class \code{"character"} Column names in input data }
    \item{\code{integrals_ci}:}{Object of class \code{"list"} Confidence intervals for integrals \eqn{S_A}, \eqn{I_A}, and \eqn{R_A} }
    \item{\code{velocity_ci}:}{Object of class \code{"list"} Confidence intervals for velocity measures \eqn{t_{FE}} and \eqn{t_{LE}} }
    \item{\code{R_0A_ci}:}{Object of class \code{"numeric"} Confidence intervals for spatial reproduction number \eqn{R_{0A}} }
    \item{\code{iterations}:}{Object of class \code{"data.frame"} Results of bootstrap sampling iterations }
    \item{\code{ci}:}{Object of class \code{"numeric"} Lower and upper confidence intervals based on user input }
    \item{\code{config}:}{Object of class \code{"list"} Configuration details for bootstrap sampling }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "sbm_ci")}: Plots the results of bootstrap confidence intervals for the Swash-Backwash Model; one figure with six plots: \eqn{S_A}, \eqn{I_A}, \eqn{R_A}, \eqn{t_{FE}}, \eqn{t_{LE}}, and \eqn{R_{0A}} }
    \item{print}{\code{signature(x = "sbm_ci")}: Prints an \code{sbm_ci} object; use \code{summary(sbm_ci)} for results }
    \item{show}{\code{signature(object = "sbm_ci")}: Prints an \code{sbm_ci} object; use \code{summary(sbm_ci)} for results }
    \item{summary}{\code{signature(object = "sbm_ci")}: Prints a summary of \code{sbm_ci} objects (bootstrap confidence intervals for Swash-Backwash Model estimates) }
	 }
}
\author{
Thomas Wieland
}

\references{
Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap.

Ramachandran KM, Tsokos CP (2021) Mathematical Statistics with Applications in R (Third Edition). Ch. 13.3.1 (Bootstrap confidence intervals). \doi{https://doi.org/10.1016/B978-0-12-817815-7.00013-0}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.
}

\examples{
showClass("sbm_ci")
}
\keyword{classes}