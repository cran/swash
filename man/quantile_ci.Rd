\name{quantile_ci}
\alias{quantile_ci}

\title{
Computing Quantiles for a given Numeric Vector
}
\description{
Computes quantiles for a given vector \code{x} and the related confidence intervals (\code{lower}, \code{upper}).
}
\usage{
quantile_ci(
  x, 
  alpha = 0.05
  )
}
\arguments{
  \item{x}{
     A numeric vector
}
  \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
}
\details{
Helper function for \code{plot(sbm_ci)}, but may be used separately.
}
\value{
A numeric vector with lower and upper quantile
}
\author{
Thomas Wieland
}

\examples{
numeric_vector <- c(1,9,5,6,3,10,20,6,9,14,3,5,8,6,11)
# any numeric vector

quantile_ci(numeric_vector)
}