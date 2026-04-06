\encoding{UTF-8}
\name{growth}
\alias{growth}

\title{
Logistic Growth Models for Regional Infections
}

\description{
Estimates \eqn{N} logistic growth models for \eqn{N} regions.
}

\usage{
growth(
  object, 
  S_iterations = 10, 
  S_start_est_method = "bisect", 
  seq_by = 10, 
  nls = TRUE,
  add_constant = 1,
  overwrite = FALSE,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{S_iterations}{
     Number of iterations for saturation value search
}
  \item{S_start_est_method}{
     Method for saturation value search, either "bisect" or "trial_and_error"
}
  \item{seq_by}{
     No of segments for the "trial_and_error" estimation of the saturation value
}
  \item{nls}{
     Nonlinear estimation? \code{TRUE} or \code{FALSE}
}
  \item{add_constant}{
     Numeric constant to be added to y if zero values occur (only relevant for OLS estimation)
}
  \item{overwrite}{
     \code{bool} argument which indicates whether the column containing cumulative cases should be overwritten if already existing
} 
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The function estimates logistic growth models for regional infections based on a \code{infpan} object. 
See \code{\link{logistic_growth}} for further details.
}

\value{
object of class \code{\link{growthmodels-class}}
}

\references{
Chowell G, Simonsen L, Viboud C, Yang K (2014) Is West Africa Approaching a Catastrophic Phase or is the 2014 Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. \emph{PLoS currents 6}. \doi{https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting the 2015 ebola challenge. \emph{Epidemics 22}, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{logistic_growth}, \link{exponential_growth}}
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

CH_covidwave1_growth <- 
  growth(infpan_CH)
summary(CH_covidwave1_growth)
# Logistic growth models for infpan object infpan_CH
}