\encoding{UTF-8}
\name{nbmatrix}
\alias{nbmatrix}

\title{
Construct Neighbourhood Matrix from Polygons
}

\description{
Building a neighbourhood matrix based on regions (polygons) with contiguous boundaries 
}

\usage{
nbmatrix(
  polygon_sf, 
  ID_col,
  row.names = NULL,
  snap = NULL,
  queen = TRUE,
  weights_style = "W",
  weights_zeropolicy = NULL,
  verbose = FALSE
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
     row.names for the \code{sf} object, passed to \code{spdep::poly2nb}
}
  \item{snap}{
     Snapping distance, passed to \code{spdep::poly2nb}
}
  \item{queen}{
     If TRUE, a single shared boundary point meets the contiguity condition; passed to \code{spdep::poly2nb}
}
  \item{weights_style}{
     weights_style can take values “W”, “B”, “C”, “U”, “minmax” and “S”; passed to the \code{style} argument of \code{spdep::nb2listw}
}
  \item{weights_zeropolicy}{
      if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors; 
      passed to the \code{zero.policy} argument of \code{spdep::nb2listw}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The function is based on \code{spdep::poly2nb}, which creates neighbours lists,
and \code{spdep::nb2listw}, which computes spatial weights for neighbours lists.
Please consult the corresponding documentation for a more detailed description of the parameters.
The input is a \code{sf} object (spatial data frame) and the results 
(contained in a \code{{nbmatrix-class}} object) include  
1) a \code{nb} list (\code{poly2nb} result),
2) a \code{listw} object (spatial weights),
3) a \code{data.frame} with a linear neighbor matrix,
4) a \code{data.frame} with the number of neighbors,
and other elements (statistics etc.)
}

\value{
object of class \code{\link{nbmatrix-class}}
}

\references{
Cliff AD, Ord JK (1981) Spatial Processes: Models & Applications. \url{https://archive.org/details/spatialprocesses0000clif}

GeoDa Center (2026) GeoDa. An Introduction to Spatial Data Science. \url{https://geodacenter.github.io/}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. \emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{nbstat}}, \code{\link{moran}}, \code{\link{getisord}}, \code{\link{gstar}}
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

Corona_nbstat <- nbstat(
  Corona_nbmat,
  link_data = RKI_Corona_counties, 
  ID_col = "AGS", 
  data_col = "EWZ", 
  func = "sum",
  verbose = TRUE
  )
# Sum of population (EWZ) of neighboring counties

plot(
  Corona_nbstat,
  main = "No. of inhabitants"
)
# Plot simple map of "EWZ"
}