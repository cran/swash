\encoding{UTF-8}
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
    \item{\code{SIR_regions}:}{Object of class \code{"data.frame"} Model result: Susceptible, infected and recovered regions over time }
    \item{\code{cases_by_date}:}{Object of class \code{"data.frame"} Total cases by date }
    \item{\code{cases_by_region}:}{Object of class \code{"data.frame"} Cumulative cases by region }
    \item{\code{input_data}:}{Object of class \code{"data.frame"} Input data }
    \item{\code{data_statistics}:}{Object of class \code{"numeric"} Diagnostics of input data }
    \item{\code{col_names}:}{Object of class \code{"character"} Original column names in input data }
  }
}
\section{Methods}{
  \describe{
    \item{confint}{\code{signature(object = "sbm")}: Creates bootstrap confidence intervals for \code{sbm} objects. }
    \item{plot}{\code{signature(x = "sbm")}: Plots the results of the Swash-Backwash Model; two plots: edges over time, total infections per time unit }
    \item{print}{\code{signature(x = "sbm")}: Prints an \code{sbm} object; use \code{summary(sbm)} for results }
    \item{show}{\code{signature(object = "sbm")}: Prints an \code{sbm} object; use \code{summary(sbm)} for results }
    \item{summary}{\code{signature(object = "sbm")}: Prints a summary of \code{sbm} objects (results of the Swash-Backwash Model) }
    \item{growth}{\code{signature(object = "sbm")}: Estimates logistic growth models from \code{sbm} objects }
    \item{growth_initial}{\code{signature(object = "sbm")}: Estimates exponential growth models from \code{sbm} objects for a given time period}
	 }
}
\author{
Thomas Wieland
}
\references{
Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. \emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}
\examples{
showClass("sbm")
}
\keyword{classes}