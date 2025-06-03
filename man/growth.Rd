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
  nls = TRUE
  )
}

\arguments{
  \item{object}{
     object of class \code{sbm}
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
}
\details{
The function estimates logistic growth models for regional infections based on a \code{sbm} object. See \code{\link{logistic_growth}} for further details.
}
\value{
\code{list} with two entries:
    \item{\code{results}:}{Object of class \code{"data.frame"} Results of the logistic growth models (coefficients and derivates)}
    \item{\code{logistic_growth_models}:}{Object of class \code{"list"} List with \eqn{N} entries for \eqn{N} growth models resp. \code{loggrowth} objects }
}
\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
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

CH_covidwave1_growth <- growth(CH_covidwave1)
CH_covidwave1_growth
# Logistic growth models for sbm object CH_covidwave1
}