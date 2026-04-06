\encoding{UTF-8}
\name{growth_hawkes}
\alias{growth_hawkes}

\title{
Hawkes Processs models for Regional Infections
}

\description{
Estimates \eqn{N} Hawkes process models for \eqn{N} regions.
}

\usage{
growth_hawkes(
  object, 
  optim_method = "L-BFGS-B",
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{optim_method}{
     \code{character} value for the optimization method. Passed to argument \code{method} in \code{stats} function \code{optim()}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The function estimates Hawkes process models for regional infections based on an \code{infpan} object. 
See \code{\link{hawkes_growth}} for further details of the estimation.
}

\value{
object of class \code{\link{growthmodels-class}}
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

CH_covidwave1_Hawkes <- 
  growth_hawkes(infpan_CH)
summary(CH_covidwave1_Hawkes)
# Hawkes process models for infpan object infpan_CH 
}