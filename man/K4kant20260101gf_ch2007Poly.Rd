\name{K4kant20260101gf_ch2007Poly}
\alias{K4kant20260101gf_ch2007Poly}
\docType{data}
\title{
Swiss Canton Boundaries
}

\description{
Geodata (sf) of Swiss canton boundaries
}


\usage{data("K4kant20260101gf_ch2007Poly")}

\format{
  A data frame with 26 observations (each one represents a Swiss canton).
  \describe{
    \item{\code{id}}{Swiss canton ID}
    \item{\code{name}}{Swiss canton name}
    \item{\code{geometry}}{Geometry}
    \item{\code{cant_abbrev}}{Swiss canton abbreviation}
)}
}

\details{
The data is included as published by the Federal Statistical Office of Switzerland (BfS)
but extended by the Geometry column and the canton abbreviation.
}

\source{
Bundesamt fuer Statistik (2026). Kartengeometrien ThemaKart - Set 2026.
\url{https://dam-api.bfs.admin.ch/hub/api/dam/assets/36390072/thumbnail?width=1980&height=1080} (accessed 2026-08-30)
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

\keyword{datasets}