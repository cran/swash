\name{hist_ci}
\alias{hist_ci}
\title{
Creating Histograms with Confidence Intervals
}
\description{
Plot of a histogram of a given vector \code{x} and the related confidence intervals (\code{lower}, \code{upper}).
}
\usage{
hist_ci(
  x, 
  alpha = 0.05,
  col_bars = "grey", 
  col_ci = "red",
  ...
  )
}
\arguments{
  \item{x}{
     A numeric vector
}
  \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
  \item{col_bars}{
     Color of bars in histogram
}
  \item{col_ci}{
     Color of lines for confidence interval
}
  \item{\dots}{
     Additional arguments passed to \code{barplot()}
}
}
\details{
Helper function for \code{plot(sbm_ci)}, but may be used separately.
}
\value{
Histogram plot, no returned value
}
\author{
Thomas Wieland
}
\examples{
numeric_vector <- c(1,9,5,6,3,10,20,6,9,14,3,5,8,6,11)
# any numeric vector

hist_ci(numeric_vector)
}