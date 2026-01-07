\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\title{Methods for Function \code{summary}}
\description{
Methods for function \code{summary}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "sbm")}}{
Prints a summary of \code{sbm} objects (results of the Swash-Backwash Model)
}

\item{\code{signature(object = "sbm_ci")}}{
Prints a summary of \code{sbm_ci} objects (bootstrap confidence intervals for Swash-Backwash Model estimates)
}

\item{\code{signature(object = "countries")}}{
Prints a summary of a \code{countries} object built with the function \code{\link{compare_countries}}
}

\item{\code{signature(object = "loggrowth")}}{
Prints a summary of a \code{loggrowth} object built with the function \code{\link{logistic_growth}}
}

\item{\code{signature(object = "expgrowth")}}{
Prints a summary of a \code{expgrowth} object built with the function \code{\link{exponential_growth}}
}

}}
\keyword{methods}