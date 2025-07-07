\encoding{UTF-8}
\name{growth_initial}
\alias{growth_initial}

\title{
Exponential Growth Models for Regional Infections
}
\description{
Estimates \eqn{N} exponential growth models for a given time period in \eqn{N} regions.
}
\usage{
growth_initial(
  object, 
  time_units = 10,
  GI = 4
  )
}

\arguments{
  \item{object}{
     object of class \code{sbm}
}
  \item{time_units}{
     \code{numeric} value for the analysis time (time units from start)
}
  \item{GI}{
     Generation interval for computing \eqn{R_0}
}
}

\details{
The function estimates exponential growth models for regional infections based on a \code{sbm} object. Such models are design for the analysis of the \emph{initial} phase of an epidemic spread. The user must state how much time units (from start) are included. See \code{\link{exponential_growth}} for further details of the estimation.
}

\value{
\code{list} with two entries:
    \item{\code{results}:}{Object of class \code{"data.frame"} Results of the exponential growth models (growth rate, basic reproduction number, doubling rate) for each region}
    \item{\code{exponential_growth_models}:}{Object of class \code{"list"} List with \eqn{N} entries for \eqn{N} exponential models}
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

\examples{
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion \%in\% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

CH_covidwave1 <- 
  swash (
    data = COVID19Cases_geoRegion, 
    col_cases = "entries", 
    col_date = "datum", 
    col_region = "geoRegion"
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

CH_covidwave1_initialgrowth_3weeks <- 
  growth_initial(
    CH_covidwave1,
    time_units = 21
    )
CH_covidwave1_initialgrowth_3weeks$results
# Exponential models for sbm object CH_covidwave1 
# initial growth in the first 3 weeks
}