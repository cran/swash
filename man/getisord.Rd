\encoding{UTF-8}
\name{getisord}
\alias{getisord}

\title{
Calculate Global Getis-Ord statistic
}

\description{
Calculating the Global Getis-Ord statistic for hotspot analysis 
from an \code{nbmatrix} object
}

\usage{
getisord(
  object,
  alternative = "greater",
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     Object of class \code{nbmatrix}
}
  \item{alternative}{
     a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided"; passed to \code{spdep::globalG.test}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The Global Getis-Ord statistic (Getis & Ord 1992) is a statistical measure of spatial association that assesses 
the spatial concentration of high or low values. 
Its calculation requires a spatial weights matrix that defines the spatial relationships between observations 
and assigns weights to neighboring observations (e.g., binary weights based on contiguity).
A typical application is, for example, investigating whether high or low incidence of infectious diseases 
is spatially clustered, that is, whether areas with high (or low) incidence tend to be surrounded by neighboring areas 
that also exhibit high (or low) incidence.
This function internally uses \code{spdep::globalG.test} for calculation (see the corresponding documentation).
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

Getis A, Ord JK (1992). The analysis of spatial association by use of distance statistics. 
\emph{Geographical Analysis} 24 (3), 189–206. \doi{https://doi.org/10.1111/j.1538-4632.1992.tb00261.x}

Ord JK, Getis A (1995) Local Spatial Autocorrelation Statistics: Distributional Issues and an Application. 
\emph{Geographical Analysis} 27, 286-306. \doi{https://doi.org/10.1111/j.1538-4632.1995.tb00912.x}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{nbmatrix}}, \code{\link{nbstat}}, \code{\link{moran}}, \code{\link{gstar}}
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

Corona_SpatialStatistics <- getisord(
  Corona_nbstat2,
  verbose = TRUE
)
# Calculating Global Getis-Ord statistic

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object
}