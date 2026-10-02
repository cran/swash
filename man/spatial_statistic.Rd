\encoding{UTF-8}
\name{spatial_statistic}
\alias{spatial_statistic}

\title{
Calculating Spatial Statistics
}

\description{
Calculation of a specified spatial statistic for the data from an infpan object.
}

\usage{
spatial_statistic(
    object,
    attribute = "Cases",
    timepoint = NULL,
    statistic = "nbstat",
    func = NULL,
    row.names = NULL,
    snap = NULL,
    queen = TRUE,
    weights_style = "W",
    weights_zeropolicy = NULL,
    verbose = FALSE,
    ...
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
\item{attribute}{
     Attribute for which the statistic is calculated, must be defined in \code{infpan} object (e.g., "Cases", "Incidence")
}
  \item{timepoint}{
     Time point for which the statistic is calculated
}
  \item{statistic}{
     Statistic to be computed
}
  \item{func}{
     If \code{statistic == "nbstat"}: Descriptive statistic (\code{FUN}) to be computed for \code{data_col} of the neighbouring regions
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
  \item{...}{
     Further arguments passed to the spatial statistics functions (see details)
} 
}

\details{
The user must choose a specific time point using the \code{timepoint} argument. 
If no time point is stated, the most current time point is chosen automatically.

The \code{...} arguments are passed to the specific spatial statistics function
which depends on what is specified by argument \code{statistic}.
See \code{\link{moran}}, \code{\link{getisord}}, \code{\link{gstar}}.
}

\value{
\code{nbmatrix} object including the spatial statistic to be calculated
}

\references{
Cliff AD, Ord JK (1981) Spatial Processes: Models & Applications. 
\url{https://archive.org/details/spatialprocesses0000clif}

GeoDa Center (2026) GeoDa. An Introduction to Spatial Data Science.
\url{https://geodacenter.github.io/}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
\emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. 
\emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}, \code{\link{nbstat}}, \code{\link{nbmatrix}}, \code{\link{moran}}, \code{\link{getisord}}, \code{\link{gstar}}
}

\examples{
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

data(K4kant20260101gf_ch2007Poly)
# Cantonal boundaries of Switzerland

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
  other_cols = c(
    "Population" = "pop"
  ), 
  verbose = TRUE
)
# Import as infections panel data set (class infpan)

infpan_CH <- add_geodata(
  infpan_CH,
  K4kant20260101gf_ch2007Poly,
  unit_col = "cant_abbrev"
)
# Adding geodata (sf) to infpan object

summary(infpan_CH)
# Summary of infpan object

timestamps(infpan_CH)
# Time stamps of infpan object

moran_incidence <-
  spatial_statistic(
    infpan_CH,
    statistic = "moran",
    randomization = FALSE
  )
# Calculating Moran's I test under normality
# Result = nbmatrix object
}