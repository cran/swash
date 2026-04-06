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
  GI = 4,
  nls = TRUE,
  nls_start = list(a = 1, b = 0.1),
  add_constant = 1,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{time_units}{
     \code{numeric} value for the analysis time (time units from start)
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
The method estimates exponential growth models for regional infections based on an \code{infpan} object. 
Such models are design for the analysis of the \emph{initial} phase of an epidemic spread. 
The user must state how much time units (from start) are included. 
See \code{\link{exponential_growth}} for further details of the estimation.
}

\value{
object of class \code{\link{growthmodels-class}}
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

infpan_CH <- load_infections_paneldata(
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion",
    other_cols = c("Population" = "pop"), 
    verbose = TRUE
  )
# Import as infections panel data set (class infpan)

CH_covidwave1_initialgrowth_3weeks <- 
  growth_initial(
    infpan_CH,
    time_units = 21
  )
summary(CH_covidwave1_initialgrowth_3weeks)
# Exponential models for infpan object CH_covidwave1 
# initial growth in the first 3 weeks
}