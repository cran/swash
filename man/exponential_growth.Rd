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
  GI = 4
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
}

\details{
This function allows the estimation of an exponential growth model. The user must specify the dependent variable (daily infections) and the time variable (time counter or date values). The estimation is performed using a linearized model as an OLS estimator. The results are the exponential growth rate \eqn{r}, basic reproduction number \eqn{R_0}, and the doubling rate.
}

\value{
\code{list} with four entries:
    \item{\code{exp_gr}:}{\code{numeric} value of exponential growth rate}
    \item{\code{R0}:}{\code{numeric} value of basic reproduction number}
    \item{\code{doubling}:}{\code{numeric} value of doubling rate}
    \item{\code{model_data}:}{\code{lm} object of regression model}
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
  COVID19Cases_geoRegion[(COVID19Cases_geoRegion$geoRegion == "ZH")
                         & (COVID19Cases_geoRegion$sumTotal > 0),]
# COVID cases for Zurich

expgrowth_BS <- exponential_growth (
  y = COVID19Cases_ZH$sumTotal[1:28], 
  t = COVID19Cases_ZH$datum[1:28] 
)
# Exponential growth model for the first 4 weeks

expgrowth_BS$doubling
# Doubling rate
}