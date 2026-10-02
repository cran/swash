\encoding{UTF-8}
\name{add_geodata}
\alias{add_geodata}

\title{
Adding Geodata to Infpan Object
}

\description{
Assigns geodata (\code{sf} object) to an \code{infpan} object.
}

\usage{
add_geodata(
    object,
    sf_geodata,
    unit_col,
    verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{sf_geodata}{
     \code{sf} object containing features (polygons) representing the study area (e.g., counties)
}
  \item{unit_col}{
     Unique ID column of the regions in the \code{sf} object
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Adds geodata (\code{sf} object) to the \code{infpan} object.
Included geodata is mandatory for plotting indicators via \code{plot_map}.
}

\value{
Input \code{infpan} object including geodata
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
\emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}, 
\code{\link{load_infections_paneldata}}, 
\code{\link{calculate_incidence}}, 
\code{\link{calculate_cum}}, 
\code{\link{calculate_rollmean}},
\code{\link{calculate_rollsum}},
\code{\link{calculate_Rt}}
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
    other_cols = c("Population" = "pop"), 
    verbose = TRUE
  )
# Import as infections panel data set (class infpan)

infpan_CH <- calculate_incidence(
  infpan_CH, 
  col_name = "incidence",
  verbose = TRUE
)
# Calculate incidence of cases as "incidence"

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

plot_map(
  infpan_CH,
  attribute = "Incidence",
  main = "COVID-19 Indicence 2020-05-31",
  breaks = c(0, 0.05, 0.1, 0.2, 1),
  verbose = TRUE
)
# Simple map of incidence for the most current date (timepoint=NULL)
}