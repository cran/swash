\encoding{UTF-8}
\name{Oesterreich_Faelle}
\docType{data}
\alias{Oesterreich_Faelle}
\title{Austria Daily COVID-19 cases by region 2020-02-26 to 2020-05-31}
\description{
  A dataset containing COVID-19 cases by region (NUTS 3) and time periods (days) for Austria (Source: BMSGPK).
}
\usage{
data(Oesterreich_Faelle)
}
\format{
  A data.frame with multiple columns:
  \describe{
    \item{NUTS3}{(character) Region for which the data was collected.}
    \item{Datum}{(Date) Date of record.}
    \item{Faelle}{(integer) Number of reported cases on this date.}
  }
}
\details{
The original data was originally published by BMSGPK at a smaller spatial scale level (political districts, "Politische Bezirke"). The data was linked to a corresponding shapefile from Statistik Austria (2022), joined to the NUTS3 level via a spatial join, and summed over the Austrian NUTS3 regions. The spatial join is based on polygon centroids of the political districts level; in cases where the centroid was outside the polygon, it was placed inside the polygon manually.
}
\source{
BMSGPK, Oesterreichisches COVID-19 Open Data Informationsportal (2022) COVID-19: Zeitliche Darstellung von Daten zu Covid19-Faellen je Bezirk. \url{https://www.data.gv.at/katalog/dataset/4b71eb3d-7d55-4967-b80d-91a3f220b60c} (retrieved 2022-06-23)

Statistik Austria (2022) Politische Bezirke. \url{https://www.data.gv.at/katalog/dataset/stat_gliederung-osterreichs-in-politische-bezirke131e2} (retrieved 2022-06-27)

Wieland T (2022) C19dNUTS: Dataset of Regional COVID-19 Deaths per 100,000 Pop (NUTS). R package v1.0.1. \doi{https://doi.org/10.32614/CRAN.package.C19dNUTS}
}
\examples{
data(Oesterreich_Faelle)
# Get Austrian COVID19 cases at NUTS 3 level
# (first wave, same final date as in Swiss data: 2020-05-31)

AT_covidwave1 <- 
  swash (
    data = Oesterreich_Faelle,
    col_cases = "Faelle",
    col_date = "Datum",
    col_region = "NUTS3"
  )
# Swash-Backwash Model for Austrian COVID19 cases
# Spatial aggregate: NUTS 3

summary(AT_covidwave1)
# Summary of model results
}