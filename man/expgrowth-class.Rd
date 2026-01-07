\encoding{UTF-8}
\name{expgrowth-class}
\docType{class}
\alias{expgrowth-class}
\alias{summary,expgrowth-method}
\alias{print,expgrowth-method}
\alias{show,expgrowth-method}

\title{Class \code{"expgrowth"}}
\description{
The class \code{"expgrowth"} contains the results of the exponential_growth() function. Use \code{summary(expgrowth)} for results summary.
}

\section{Objects from the Class}{
Objects can be created by the function \code{exponential_growth}.
}

\section{Slots}{
  \describe{
  \item{\code{exp_gr}:}{Object of class \code{numeric} Exponential growth rate}
    \item{\code{y_0}:}{Object of class \code{numeric} Baseline}
    \item{\code{R0}:}{Object of class \code{numeric} Basic reproduction number, derived from exponential growth model}
    \item{\code{doubling}:}{Object of class \code{numeric} Doubling rate, derived from exponential growth model}
    \item{\code{model_data}:}{Object of class \code{lm} OLS model}
  \item{\code{config}:}{Object of class \code{list} Model fit configurations}
  }
}

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "expgrowth")}: Prints a summary of \code{expgrowth} objects}
    \item{print}{\code{signature(x = "expgrowth")}: Prints an \code{expgrowth} object; use \code{summary(expgrowth)} for results }
	 }
}

\author{
Thomas Wieland
}

\references{
Bonifazi G et al. (2021) A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. \emph{The European Physical Journal Plus} 136, 386. \doi{https://doi.org/10.1140/epjp/s13360-021-01339-6}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. \emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting the 2015 ebola challenge. \emph{Epidemics 22}, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
}

\examples{
showClass("expgrowth")
}

\keyword{classes}