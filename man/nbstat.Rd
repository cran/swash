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
  object,
  link_data,
  ID_col,
  data_col,
  func = "sum",
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     \code{nbmatrix} object
}
  \item{link_data}{
     \code{data.frame} to merge with
}
  \item{ID_col}{
     Column with unique ID of each polygon feature in the link_data \code{data.frame}
}
  \item{data_col}{
     Column with regarded numeric values in \code{data.frame}
}
  \item{func}{
     Descriptive statistic (\code{FUN}) to be computed for \code{data_col} of the neighbouring regions 
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Calculates a specified statistic (e.g., mean or sum) for the neighboring regions within a \code{nbmatrix} object.
The user must specify a (spatial) data frame with the parameter \code{link_data}.
The results are saved in the \code{nbmatrix} object.
If \code{func} is set to NULL, \code{nbstat} is used to specify the data column used for further analyses (e.g., \code{moran()}).
}

\value{
object of class \code{\link{nbmatrix-class}}
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
\code{\link{nbmatrix}}, \code{\link{moran}}, \code{\link{getisord}}, \code{\link{gstar}}
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