\encoding{UTF-8}
\name{binary_metrics}
\alias{binary_metrics}

\title{
Fit metrics of observed and expected binary variables
}
\description{
Calculation of fit metrices for binary variables (Sensitivity, specificity, accuracy) 
}
\usage{
binary_metrics(
  observed, 
  expected,
  no_information_rate = "negative"
  )
}

\arguments{
  \item{observed}{
     Numeric vector: Y observed
}
  \item{expected}{
     Numeric vector: Y expected
}
  \item{no_information_rate}{
     \code{bool} argument which indicates whether the no-information rate is calculated based on negatives or positives
}
}

\details{
The function computes model performance metrices for binary outcomes. Observed and expected data must be stated by the user. The function returns sensitivity, specificity, accurracy, and no-information rate.
}

\value{
\code{list} with two entries:
    \item{\code{fit_metrics}:}{\code{list} with fit metrics (\code{sens}, \code{spec}, ...)}
    \item{\code{observed_expected}:}{\code{data.frame} with observed, expected and hit (1/0)}
}

\references{
Altman DG, Bland JM (1994) Diagnostic tests. 1: Sensitivity and specificity. \emph{British Medical Journal} 308, 1552. \doi{https://doi.org/10.1136/bmj.308.6943.1552}.

Boehmke B, Greenwell B (2020) Hands-On Machine Learning with R (1 ed.). Taylor & Francis, New York, NY.
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{metrics}, \link{binary_metrics_glm}}
}

\examples{
obs <- c(1,1,0,0,0,0,1,0,1)
exp <- c(0,1,0,0,0,0,1,0,0)

binary_metrics(
  obs,
  exp
)
}