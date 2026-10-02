\encoding{UTF-8}
\name{moran}
\alias{moran}

\title{
Calculate Global Moran's I statistic
}

\description{
Calculating the Global Moran's I statistic for spatial autocorrelation 
from an \code{nbmatrix} object
}

\usage{
moran(
  object,
  randomization = TRUE,
  alternative = "greater",
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     Object of class \code{nbmatrix}
}
  \item{randomization}{
     \code{bool} argument that indicates which test is conducted (under randomization or normality); passed to \code{spdep::moran.test}
}
  \item{alternative}{
     a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided; passed to \code{spdep::moran.test}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Moran's I (Moran 1950) is a statistical measure of spatial autocorrelation. 
Its calculation requires a spatial weights matrix that defines the spatial relationships between observations 
and assigns weights to neighboring observations (e.g., binary weights based on contiguity).
A typical application is, for example, investigating whether spatial autocorrelation exists 
in the incidence of infectious diseases, that is, whether areas with high (or low) incidence are surrounded by 
neighboring areas that also exhibit high (or low) incidence (Wieland 2020, 2022).
This function internally uses \code{spdep::moran.test} for calculation (see the corresponding documentation).
The results are written in the \code{nbmatrix} instance and displayed in the summary.
}

\value{
object of class \code{\link{nbmatrix-class}}
}

\references{
Cliff AD, Ord JK (1981) Spatial Processes: Models & Applications. 
\url{https://archive.org/details/spatialprocesses0000clif}

GeoDa Center (2026) GeoDa. An Introduction to Spatial Data Science. 
\url{https://geodacenter.github.io/}

Moran PAP (1950) Notes on Continuous Stochastic Phenomena. 
\emph{Biometrika} 37(1/2) 17–23. 
\doi{https://doi.org/10.2307/2332142}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
\emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. 
\emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{nbmatrix}}, \code{\link{nbstat}}, \code{\link{getisord}}, \code{\link{gstar}}
}

\examples{
data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

Corona_nbmat <- 
  nbmatrix (
    RKI_Corona_counties, 
    ID_col="AGS",
    verbose = TRUE
  )
# Creating neighborhood matrix

Corona_nbstat2 <- nbstat(
  Corona_nbmat,
  link_data = RKI_Corona_counties, 
  ID_col = "AGS", 
  data_col = "cases_per_", 
  func = NULL,
  verbose = TRUE
)
# New nbmatrix instance
# Defining "cases_per_" (COVID-19 cases per 100,000 inhabitants) as analysis column

Corona_SpatialStatistics <- moran(
  Corona_nbstat2,
  verbose = TRUE
  )
# Calculating Global Moran's I

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object
}