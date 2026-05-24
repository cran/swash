\encoding{UTF-8}
\name{breaks_growth}
\alias{breaks_growth}

\title{
Time Series Model with Breakpoints
}

\description{
Estimation of breakpoints in linear regression models from incremental (e.g., daily) infections data
}

\usage{
breaks_growth(
  y, 
  t,
  ln = FALSE,
  add_constant = 1,
  alpha = 0.05,
  ...,
  verbose = FALSE
  )
}

\arguments{
  \item{y}{
     \code{numeric vector} with cumulative infections data over time
}
  \item{t}{
     \code{vector} of class \code{numeric} or \code{Date} with time points or dates
}
  \item{ln}{
     \code{bool} argument which indicates whether dependent variable should be transformed by natural logarithm
}
  \item{add_constant}{
     Numeric constant to be added to y if zero values occur
}
  \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
  \item{...}{
     Other parameters passed to \code{strucchange::breakpoints} (see the corresponding documentation)
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
This function allows detecting breakpoints in a linear regression time series model. 
The user must specify the dependent variable (e.g, daily infections) and the time variable (time counter or date values). 
The estimation is performed using OLS.
The function internally uses the function \code{breakpoints} from the \code{strucchange} package (Zeileis et al. 2003),
where breakpoints are identified using the Bai-Perron algorithm (Bai & Perron 2003). 
}

\value{
object of class \code{\link{breaksgrowth-class}}
}

\references{
Bai J, Perron P (2003) Computation and analysis of multiple structural change models. \emph{Journal of Applied Econometrics} 18(1), 1-22. \doi{https://doi.org/10.1002/jae.659}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Zeileis C, Kleiber W, Krämer K, Hornik K (2003) Testing and dating of structural changes in practice. \emph{Computational Statistics & Data Analysis} 44(1-2), 109-123. \doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{exponential_growth}, \link{logistic_growth}}, \link{hawkes_growth}
}

\examples{
data(Infections)
# Confirmed SARS-CoV-2 cases in Germany

breakpoints_infections <- breaks_growth(
  y = Infections$infections_daily,
  t = Infections$day,
  ln = TRUE,
  verbose = TRUE
)
# Breakpoints for time series of infections

summary(breakpoints_infections)
# Summary of breakpoints

plot(breakpoints_infections)
# Plot breakpoints
}