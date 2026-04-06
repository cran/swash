\encoding{UTF-8}
\name{calculate_incidence}
\alias{calculate_incidence}

\title{
Incidence from Infection Numbers
}

\description{
Calculation of the incidence from infection numbers and population for infections panel data.
}

\usage{
calculate_incidence(
  object,
  use_column = NULL,
  col_name = NULL,
  pop_factor = 100000,
  overwrite = FALSE,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{use_column}{
     \code{character} value specifying which column should be used for incidence calculation
}
  \item{col_name}{
     \code{character} value specifying the column name of the computed incidence
}
  \item{pop_factor}{
     \code{numeric} value specifying the factor with which the incidence should be multiplied (e.g., cases/pop*100000)
}  
  \item{overwrite}{
     \code{bool} argument which indicates whether the column should be overwritten if already existing
}  
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Calculates the incidence of the infections panel data for all time points for each region.
Use \code{use_column} to specify which column should be used for the calculation of incidence. 
The following values are permitted: "Cases" (default, incremental cases), "Cum. cases" (cumulative cases), "Roll. mean" (rolling mean of cases), 
or "Roll. sum" (rolling sum of cases). If the specified column does not exist in the infections panel data of the \code{infpan} object,
the function raises an error.
If in the \code{infpan} object, no "Population" column is defined, incidence calculation is not possible.
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_inc".
Set \code{overwrite} to TRUE, if an existing column should be overwritten.
}

\value{
\code{infpan} object including column with incidence values in the infections panel data
}

\references{
an der Heiden M, Hamouda O (2020) Schätzung der aktuellen Entwicklung der SARS-CoV-2-Epidemie in Deutschland - Nowcasting. \emph{Epidemiologisches Bulletin} 17, 10-15. \doi{https://doi.org/10.25646/6692}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}, \code{\link{load_infections_paneldata}}, \code{\link{calculate_Rt}}, \code{\link{calculate_cum}}, \code{\link{calculate_rollmean}}
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

infpan_CH <- calculate_incidence(
  infpan_CH, 
  col_name = "incidence",
  verbose = TRUE
)
# Calculate incidence of cases as "incidence"

summary(infpan_CH)
# Summary of infpan object
}