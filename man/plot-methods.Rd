\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\title{Methods for Function \code{plot}}
\description{
Methods for function \code{plot}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "sbm")}}{
Plots the results of the Swash-Backwash Model; two plots: edges over time, total infections per time unit
}

\item{\code{signature(x = "sbm_ci")}}{
Plots the results of bootstrap confidence intervals for the Swash-Backwash Model; one figure with six plots: \eqn{S_A}, \eqn{I_A}, \eqn{R_A}, \eqn{t_{FE}}, \eqn{t_{LE}}, and \eqn{R_{0A}}
}

\item{\code{signature(x = "loggrowth")}}{
Plots the results of the logistic growth model (observed, predicted, first derivative)
}

}}
\author{
Thomas Wieland
}
\keyword{methods}