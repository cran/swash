\encoding{UTF-8}
\name{nbmatrix}
\alias{nbmatrix}

\title{
Construct Neighbourhood Matrix from Polygons
}
\description{
Building a neighbourhood matrix based on regions (polygons) with contiguous boundaries and resulting a data frame 
}
\usage{
nbmatrix(
  polygon_sf, 
  ID_col,
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
  \item{row.names}{
     row.names for the \code{sf} object
}
}

\details{
The function is based on \code{spdep::poly2nb}, which creates neighbours lists. The input is a \code{sf} object (spatial data frame) and the results are 1) a \code{nb} list (\code{poly2nb} result) and 2) a \code{data.frame}.
}

\value{
\code{list} with two entries:
    \item{\code{nb}:}{Object of class \code{"sb"} Neighbours list; see the \code{spdep:poly2nb} documentation }
    \item{\code{nbmat}:}{Object of class \code{"data.frame"} Dataset neighbouring regions }
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. \emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{nbstat}}
}

\examples{
data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

Corona_nbmat <- 
  nbmatrix (
    RKI_Corona_counties, 
    ID_col="AGS"
  )
# Creating neighborhood matrix
}