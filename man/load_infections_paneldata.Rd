\encoding{UTF-8}
\name{load_infections_paneldata}
\alias{load_infections_paneldata}

\title{
Import of infections panel data
}

\description{
Loading infections panel data (data.frame) and creating an object of class infpan
}

\usage{
load_infections_paneldata(
  data,
  col_cases, 
  col_date, 
  col_region,
  other_cols = NULL,
  time_format = "\%Y-\%m-\%d",
  time_unit = "days",
  verbose = FALSE
  )
}

\arguments{
  \item{data}{
     \code{data.frame} with regional infection data
}
  \item{col_cases}{
     \code{character}, Column containing the cases (numeric)
}
  \item{col_date}{
     \code{character}, Column containing the time points (e.g., days)
}
  \item{col_region}{
    \code{character}, Column containing the unique identifier of the regions (e.g., name, NUTS 3 code)
}
  \item{other_cols}{
     \code{list}, Further columns in the input data
}
  \item{time_format}{
     \code{character}, Time format of the values in \code{col_date}
}
\item{time_unit}{
     \code{character}, Time unit of the values in \code{col_date}, e.g., "days"
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The function import user-given infections panel data. 
The input data is checked in several ways (e.g., whether data is balanced or not).
Other relevant columns from the input data may be defined in the \code{character} vector \code{other_cols}: 
"R_t" (Effective reproduction number), "Cum. cases" (Cumulative cases), "Incidence" Incidence (per xxx pop),
"Population" (Population size of the region), "Roll. mean" (Rolling mean of cases), and "Roll. sum" (Rolling sum of cases).

The output is an object of class \code{infpan}. The results can be viewed using \code{summary(infpan)}.
From an instance of class \code{infpan}, all built-in analyses for infections panel data may be conducted,
e.g., the Swash-Backwash Model (\code{swash(infpan)}) or logistic growth models (\code{growth(infpan)}).
}

\value{
object of class \code{\link{infpan-class}}
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}
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

infpan_CH <- load_infections_paneldata(
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion",
    other_cols = c("Population" = "pop"), 
    verbose = TRUE
  )
# Import as infections panel data set (class infpan)

CH_covidwave1 <-
  swash(
    infpan_CH,
    verbose = TRUE
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model
}