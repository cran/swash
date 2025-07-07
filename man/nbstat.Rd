\encoding{UTF-8}
\name{nbstat}
\alias{nbstat}

\title{
Calculate Neighbourhood Statistics from Polygons
}
\description{
Calculating descriptive neighbourhood statistics based on regions (polygons) with contiguous boundaries and resulting a data frame
}
\usage{
nbstat(
  polygon_sf, 
    ID_col, 
    link_data, 
    data_ID_col, 
    data_col, 
    func = "sum",
    row.names = NULL
  )
}

\arguments{
  \item{polygon_sf}{
     \code{sf} object with polygons
}
  \item{ID_col}{
     Column of \code{polygon_sf} with unique ID of each polygon
}
  \item{link_data}{
     \code{data.frame} to merge with
}
  \item{data_ID_col}{
     Column with unique ID of each polygon in \code{data.frame}
}
  \item{data_col}{
     Column with regarded numeric values in \code{data.frame}
}
  \item{func}{
     Descriptive statistic (\code{FUN}) to be computed for \code{data_col} of the neighbouring regions 
}
  \item{row.names}{
     row.names for the \code{sf} object
}
}

\details{
The function is based on \code{spdep::poly2nb}, which creates neighbours lists. The input is a \code{sf} object (spatial data frame) and the results are 1) a \code{nb} list (\code{poly2nb} result) and 2) a \code{data.frame}.
}

\value{
\code{list} with three entries:
    \item{\code{nbmat}:}{Object of class \code{"data.frame"} Dataset neighbouring regions }
    \item{\code{nbmat_data}:}{Object of class \code{"data.frame"} Dataset neighbouring regions and linked data}
    \item{\code{nbmat_data_aggreagte}:}{Object of class \code{"data.frame"} Dataset with statistic by region}
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. \emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{nbmatrix}}
}

\examples{
data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

Corona_nbstat <- 
  nbstat (
    RKI_Corona_counties, 
    ID_col="AGS",
    link_data = RKI_Corona_counties, 
    data_ID_col = "AGS", 
    data_col = "EWZ", 
    func = "sum"
  )
Corona_nbstat$nbmat_data_aggregate
# Sum of population (EWZ) of neighbouring counties
}