\encoding{UTF-8}
\name{hawkes-class}
\docType{class}
\alias{hawkes-class}
\alias{summary,hawkes-method}
\alias{print,hawkes-method}
\alias{show,hawkes-method}
\alias{plot,hawkes-method}

\title{Class \code{"hawkes"}}
\description{
The class \code{"hawkes"} contains the results of the hawkes_growth() function. 
Use \code{summary(hawkes)} for results summary.
}

\section{Objects from the Class}{
Objects can be created by the function \code{hawkes_growth}.
}

\section{Slots}{
  \describe{
  \item{\code{t}:}{Object of class \code{numeric} Input time points data}
  \item{\code{y}:}{Object of class \code{numeric} Input infections data}
  \item{\code{mu}:}{Object of class \code{numeric} Estimated \eqn{\mu} parameter}
  \item{\code{alpha}:}{Object of class \code{numeric} Estimated \eqn{\alpha} parameter}
  \item{\code{beta}:}{Object of class \code{numeric} Estimated \eqn{\beta} parameter}
  \item{\code{br}:}{Object of class \code{numeric} Estimated breaking ratio (\eqn{\alpha}/\eqn{\beta})}
  \item{\code{y_pred}:}{Object of class \code{numeric} Predicted values of y}
  \item{\code{fit_metrics}:}{Object of class \code{list} Fit metrics for model, output from built-in function \code{fit_metrics}}
  \item{\code{config}:}{Object of class \code{list} Model fit configurations}
  }
}

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "hawkes")}: Prints a summary of \code{hawkes} objects}
    \item{print}{\code{signature(x = "hawkes")}: Prints an \code{hawkes} object; use \code{summary(hawkes)} for results }
    \item{plot}{\code{signature(x = "hawkes")}: Plots the results of the Hawkes model (observed, predicted) }
	 }
}

\author{
Thomas Wieland
}

\references{
Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L. (2018) SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. In: \emph{Proceedings of the 2018 World Wide Web Conference. 
WWW’18}. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. \doi{https://doi.org/10.1145/3178876.3186108}
}

\examples{
showClass("hawkes")
}

\keyword{classes}