\name{plot_regions}
\alias{plot_regions}

\title{
Plots of Regional Infections Over Time
}

\description{
Plots regional infection curves in \eqn{N} plots for \eqn{N} regions.
}

\usage{
plot_regions(
  object, 
  col = "red",
  lty = "solid",
  scale = FALSE, 
  normalize_by_col = NULL, 
  normalize_factor = 1,
  plot_rollmean = FALSE,
  rollmean_col = "blue",
  rollmean_lty = "solid",
  rollmean_k = 7,
  rollmean_align = "center",
  rollmean_fill = NA,
  growth_col = "orange",
  growth_lty = "solid", 
  growth_per_time_unit = 1
  )
}

\arguments{
  \item{object}{
     object of class \code{sbm}
}
  \item{col}{
     Color of line plot
}
\item{lty}{
     Type of line plot
}
  \item{scale}{
     Set y axis of the plots uniformly based on the maximum value across all regions? (boolean, default: \code{FALSE})
}
  \item{normalize_by_col}{
     Normalize infection numbers by stating a column in the input data frame (e.g., regional population)
}
  \item{normalize_factor}{
     Multiply density/incidence with a factor (say, 100,000 inhabitants)
}
 \item{plot_rollmean}{
     Plot rolling mean?
}
 \item{rollmean_col}{
     Color of rolling mean line
}
 \item{rollmean_lty}{
     Type of rolling mean line
}
\item{rollmean_k}{
     integer width of the rolling window (see \code{zoo::rollmean})
}
\item{rollmean_align}{
     specifying whether the index of the result should be left- or right-aligned or centered (default) compared to the rolling window of observations (see \code{zoo::rollmean})
}
\item{rollmean_fill}{
      filling values at the left/within/to the right of the data range (see \code{zoo::rollmean})
}
\item{growth_col}{
      Color of growth line
}
\item{growth_lty}{
      Line type of growth line
}
\item{growth_per_time_unit}{
      Time unit for growth rates
}
}

\details{
Plots regional infection curves in \eqn{N} plots for \eqn{N} regions, with the number of columns equals 4 and the number of rows is calculated based the size of \eqn{N}.
}

\value{
Plot only, no returned value
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

CH_covidwave1 <- 
  swash (
    data = COVID19Cases_geoRegion, 
    col_cases = "entries", 
    col_date = "datum", 
    col_region = "geoRegion"
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

plot_regions(CH_covidwave1)
# Plot of regional infections
}