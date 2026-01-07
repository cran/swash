\encoding{UTF-8}
\name{swash}
\alias{swash}

\title{
Swash-Backwash Model for the Single Epidemic Wave
}
\description{
Analysis of regional infection/surveillance data using the Swash-Backwash Model for the single epidemic wave by Cliff and Haggett (2006)
}
\usage{
swash(
  data, 
  col_cases, 
  col_date, 
  col_region
  )
}

\arguments{
  \item{data}{
     \code{data.frame} with regional infection data
}
  \item{col_cases}{
     Column containing the cases (numeric)
}
  \item{col_date}{
     Column containing the time points (e.g., days)
}
  \item{col_region}{
     Column containing the unique identifier of the regions (e.g., name, NUTS 3 code)
}
}
\details{
The function performs the analysis of the input panel data using the Swash-Backwash Model. The output is an object of class \code{"sbm"}. The results can be viewed using \code{summary(sbm)}. The user must state panel data with daily infections.
}

\value{
object of class \code{\link{sbm-class}}
}

\references{
Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{sbm-class}}
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

summary(CH_covidwave1)
# Summary of Swash-Backwash Model

plot(CH_covidwave1)
# Plot of Swash-Backwash Model edges and total epidemic curve
}