\encoding{UTF-8}
\name{calculate_cum}
\alias{calculate_cum}

\title{
Cumulative Infection Numbers
}

\description{
Calculation of the cumulative values of infection numbers for infections panel data.
}

\usage{
calculate_cum(
  object,
  col_name = NULL,
  overwrite = FALSE,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{col_name}{
     \code{character} value specifying the column name of the computed cumulative values
}
  \item{overwrite}{
     \code{bool} argument which indicates whether the column should be overwritten if already existing
}  
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Calculates the cumulative values of the infections panel data for all time points for each region.
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_cum".
Set \code{overwrite} to TRUE, if an existing column should be overwritten.
The method uses the function \code{cumsum} from the \code{base} package (see the corresponding documentation).
}

\value{
\code{infpan} object including column with cumulative values in the infections panel data
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}, \code{\link{load_infections_paneldata}}, \code{\link{calculate_Rt}}, \code{\link{calculate_incidence}}, \code{\link{calculate_rollmean}}
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

infpan_CH <- calculate_cum(
  infpan_CH, 
  col_name = "cumulatives",
  verbose = TRUE
)
# Calculate rolling mean of cases as "cumulatives"

summary(infpan_CH)
# Summary of infpan object
}