\encoding{UTF-8}
\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\title{Methods for Function \code{confint}}
\description{
Methods for function \code{confint}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "sbm", iterations = 100, samples_ratio = 0.8, alpha = 0.05, replace = TRUE)}}{
Creates bootstrap confidence intervals for \code{sbm} objects. The argument \code{iterations} indicates the number of bootstrap samples which are drawn. Since the initial data in the Swash-Backwash Model should be balanced, \emph{entity-based bootstrap sampling} is carried out. This means that not, for example, 80\% of all observations are included in each sample at a sample ratio equal to \eqn{p} = 0.8 (\code{samples_ratio = 0.8}), but rather all observations for 80\% of the regions. The significance level for the confidence intervals \eqn{\alpha} is set by the argument \code{alpha} (default: 0.05, which corresponds to a 95\% confidence level). 
}
}}
\references{

\emph{Swash-Backwash Model:}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

\emph{Bootstrapping und bootstrap confidence intervals:}

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap.

Ramachandran KM, Tsokos CP (2021) Mathematical Statistics with Applications in R (Third Edition). Ch. 13.3.1 (Bootstrap confidence intervals). \doi{https://doi.org/10.1016/B978-0-12-817815-7.00013-0}
}
\author{
Thomas Wieland
}
\seealso{
\code{\link{sbm_ci-class}}
}

\keyword{methods}