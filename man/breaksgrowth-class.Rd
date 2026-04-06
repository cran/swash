\encoding{UTF-8}
\name{breaksgrowth-class}
\docType{class}
\alias{breaksgrowth-class}
\alias{summary,breaksgrowth-method}
\alias{print,breaksgrowth-method}
\alias{show,breaksgrowth-method}
\alias{plot,breaksgrowth-method}

\title{Class \code{"breaksgrowth"}}
\description{
The class \code{"breaksgrowth"} contains the results of the breaks_growth() function. 
Use \code{summary(breaksgrowth)} for results summary.
}

\section{Objects from the Class}{
Objects can be created by the function \code{breaks_growth}.
}

\section{Slots}{
  \describe{
  \item{\code{GrowthModel_OLS}:}{Object of class \code{list} Results of the OLS fit (predicted, parameters)}
  \item{\code{t}:}{Object of class \code{numeric} Input time points data}
  \item{\code{y}:}{Object of class \code{numeric} Input infections data}
  \item{\code{config}:}{Object of class \code{list} Model fit configurations}
  }
}

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "breaksgrowth")}: Prints a summary of \code{breaksgrowth} objects}
    \item{plot}{\code{signature(x = "breaksgrowth")}: Plots the results of the breakpoint analysis }
    \item{print}{\code{signature(x = "breaksgrowth")}: Prints an \code{breaksgrowth} object; use \code{summary(breaksgrowth)} for results }
	 }
}

\author{
Thomas Wieland
}

\references{
Bai J, Perron P (2003) Computation and analysis of multiple structural change models. \emph{Journal of Applied Econometrics} 18(1), 1-22. \doi{https://doi.org/10.1002/jae.659}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Zeileis C, Kleiber W, Krämer K, Hornik, K (2003) Testing and dating of structural changes in practice. \emph{Computational Statistics & Data Analysis} 44(1-2), 109-123. \doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
}

\examples{
showClass("breaksgrowth")
}

\keyword{classes}