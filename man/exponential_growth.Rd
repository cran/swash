\encoding{UTF-8}
\name{exponential_growth}
\alias{exponential_growth}

\title{
Exponential Growth Model for Epidemic Data
}
\description{
Estimation of exponential growth models from daily infections data
}
\usage{
exponential_growth(
  y, 
  t, 
  GI = 4,
  nls = TRUE,
  nls_start = list(a = 1, b = 0.1),
  add_constant = 1,
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
  \item{GI}{
     Generation interval for computing \eqn{R_0}
}
 \item{nls}{
     Nonlinear estimation? \code{TRUE} or \code{FALSE}
}
 \item{nls_start}{
     A \code{list} with start values for the two parameters to be estimated
}
  \item{add_constant}{
     Numeric constant to be added to y if zero values occur (only relevant for OLS estimation)
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
This function allows the estimation of an exponential growth model. 
The user must specify the dependent variable (daily infections) and the time variable (time counter or date values). 
The estimation is performed using a linearized model as an OLS estimator, and, if \code{nls=TRUE}, also by NLS. 
The results are the exponential growth rate \eqn{r}, basic reproduction number \eqn{R_0}, and the doubling rate.
}

\value{
object of class \code{\link{expgrowth-class}}
}


\references{
Bonifazi G et al. (2021) A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. \emph{The European Physical Journal Plus} 136, 386. \doi{https://doi.org/10.1140/epjp/s13360-021-01339-6}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. \emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting the 2015 ebola challenge. \emph{Epidemics 22}, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{growth_initial}, \link{logistic_growth}}
}

\examples{
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_ZH <-
  COVID19Cases_geoRegion[
  (COVID19Cases_geoRegion$geoRegion == "ZH")
  & (COVID19Cases_geoRegion$sumTotal > 0)
  ,]
# COVID cases for Zurich

expgrowth_ZH <- exponential_growth(
  y = COVID19Cases_ZH$sumTotal[1:28], 
  t = COVID19Cases_ZH$datum[1:28] 
)
# Exponential growth model for the first 4 weeks

summary(expgrowth_ZH)
# Summary of exponential growth model

plot(expgrowth_ZH)
# Plot of exponential growth model

expgrowth_ZH@GrowthModel_OLS$fit_metrics
expgrowth_ZH@GrowthModel_OLS$fit_metrics
# Fit metrics for OLS and NLS models
}