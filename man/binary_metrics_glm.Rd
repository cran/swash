\encoding{UTF-8}
\name{binary_metrics_glm}
\alias{binary_metrics_glm}

\title{
Fit metrics for binary logit model
}
\description{
Calculation of fit metrices for binary variables (Sensitivity, specificity, accuracy) out of binary logit models (glm object) 
}
\usage{
binary_metrics_glm(
  logit_model, 
  threshold = 0.5
  )
}

\arguments{
  \item{logit_model}{
     \code{glm} object with binary logit model
}
  \item{threshold}{
     Threshold for destinction of probability with respect to TRUE or FALSE
}
}

\details{
The function computes model performance metrices for binary outcomes. A binary logit model (glm) must be stated by the user. The function returns sensitivity, specificity, accurracy, and no-information rate.
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
\code{\link{metrics}, \link{binary_metrics}}
}

\examples{
dep <- c(1,1,0,0,0,0,1,0,1, 1)
x <- c(2,3,1,1,0,1,3,2,1,3)

testmodel <-
  glm(
    dep~x,
    family=binomial()
  )
  
summary(testmodel)

binary_metrics_glm(testmodel)
}