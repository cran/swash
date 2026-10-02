\encoding{UTF-8}
\name{gstar}
\alias{gstar}

\title{
Calculate Local Getis-Ord Gi* statistic
}

\description{
Calculating the Local Getis-Ord Gi* for detecting significant hotspots 
from an \code{nbmatrix} object
}

\usage{
gstar(
  object,
  alternative = "greater",
  nsim = 999,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     Object of class \code{nbmatrix}
}
  \item{alternative}{
     a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided"; passed to \code{sfdep::local_gstar_perm}
}
  \item{nsim}{
     \code{numeric} value with number of simulations to run; passed to \code{sfdep::local_gstar_perm}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Unlike the global Getis-Ord statistic, the Getis-Ord \eqn{G_i^*} statistic (Getis & Ord 1992) is calculated for each region in the study. 
Based on the calculated \eqn{z}-scores and their \eqn{p}-values, the analysis identifies where 
regions with high or low values tend to form spatial clusters.
A typical application is, for example, investigating whether high or low incidence of infectious diseases 
is spatially clustered, that is, identifying hotspots and/or coldspots.
This function internally uses \code{sfdep::local_gstar_perm.test} for calculation (see the corresponding documentation).
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
\code{\link{nbmatrix}}, \code{\link{getisord}}, \code{\link{nbstat}}, \code{\link{moran}}
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

Corona_SpatialStatistics <- gstar(
  Corona_nbstat2,
  verbose = TRUE
  )
# Calculating Local Getis-Ord Gi*

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object
}