\encoding{UTF-8}
\name{loggrowth-class}
\docType{class}
\alias{loggrowth-class}
\alias{plot,loggrowth-method}
\alias{summary,loggrowth-method}
\alias{print,loggrowth-method}
\alias{show,loggrowth-method}

\title{Class \code{"loggrowth"}}
\description{
The class \code{"loggrowth"} contains the results of the logistic_growth() function. Use \code{summary(loggrowth)} and \code{plot(loggrowth)} for results summary and plotting, respectively.
}
\section{Objects from the Class}{
Objects can be created by the function \code{logistic_growth}.
}
\section{Slots}{
  \describe{
    \item{\code{LinModel}:}{Object of class \code{list} Results of the OLS helper model}
    \item{\code{GrowthModel_OLS}:}{Object of class \code{list} Results of the OLS fit (predicted, parameters, first derivative)}
    \item{\code{GrowthModel_NLS}:}{Object of class \code{list} Results of the NLS fit (predicted, parameters, first derivative)}
    \item{\code{t}:}{Object of class \code{numeric} Input time points data}
    \item{\code{y}:}{Object of class \code{numeric} Input infections data}
  \item{\code{config}:}{Object of class \code{list} Model fit configurations}

  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "loggrowth")}: Plots the results of the logistic growth model (observed, predicted, first derivative) }
    \item{summary}{\code{signature(object = "loggrowth")}: Prints a summary of \code{loggrowth} objects}
    \item{print}{\code{signature(x = "loggrowth")}: Prints an \code{loggrowth} object; use \code{summary(loggrowth)} for results }
	 }
}

\author{
Thomas Wieland
}

\references{
Chowell G, Simonsen L, Viboud C, Yang K (2014) Is West Africa Approaching a Catastrophic Phase or is the 2014 Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. \emph{PLoS currents 6}. \doi{https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting the 2015 ebola challenge. \emph{Epidemics 22}, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}

\examples{
showClass("loggrowth")
}

\keyword{classes}