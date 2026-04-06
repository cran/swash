\encoding{UTF-8}
\name{hawkes_growth}
\alias{hawkes_growth}

\title{
Hawkes Process Model for Epidemic Data
}

\description{
Estimation of Hawkes Process models from incremental infections data
}

\usage{
hawkes_growth(
  y,
  optim_method = "L-BFGS-B",
  verbose = FALSE
  )
}

\arguments{
  \item{y}{
     \code{numeric vector} with incremental infections data over time (e.g., daily infections)
}
  \item{optim_method}{
     \code{character} specifying the optimization algorithm, passed to \code{stats::optim}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
This function allows the estimation of a Hawkes Process model, with the time decay being expressed as exponential function,
which results in three estimated parameters (\eqn{\mu}, \eqn{\alpha}, and \eqn{\beta}).
The user must specify the dependent variable (incremental infections). 
The estimation is performed using nonlinear estimation via \code{stats::optim}.
See the corresponding documentation for available optimization methods (default: "L-BFGS-B").
}

\value{
object of class \code{\link{hawkes-class}}
}

\references{
Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L. (2018) SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. In: \emph{Proceedings of the 2018 World Wide Web Conference. 
WWW’18}. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. \doi{https://doi.org/10.1145/3178876.3186108}
}

\author{
Thomas Wieland
}

\seealso{
\link{logistic_growth}, \link{exponential_growth}, \link{breaks_growth}
}

\examples{
data(Infections)
# Confirmed SARS-CoV-2 cases in Germany

hawkes_BS <- hawkes_growth(
  y = Infections$infections_daily
)
# Hawkes Process model

summary(hawkes_BS)
# Summary of Hawkes model estimates

plot(hawkes_BS)
# Plot of Hawkes Process model
}