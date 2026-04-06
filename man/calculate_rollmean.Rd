\encoding{UTF-8}
\name{calculate_rollmean}
\alias{calculate_rollmean}

\title{
Rolling Means of Infection Numbers
}

\description{
Calculation of the rolling means of infection numbers for infections panel data.
}

\usage{
calculate_rollmean(
  object,
  k = 7,
  align = "center",
  fill = NA,
  col_name = NULL,
  overwrite = FALSE,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{k}{
     \code{integer} width of the rolling window (default: 7)
}
  \item{align}{
     \code{character} specifying whether the rolling mean should be left- or right-aligned or centered (default) compared to the rolling window (default: center)
}
  \item{fill}{
     \code{numeric} value or \code{NA} for the filling value at the left/within/right end of the data range
}
  \item{col_name}{
     \code{character} value specifying the column name of the computed rolling means
}
  \item{overwrite}{
     \code{bool} argument which indicates whether the column should be overwritten if already existing
}  
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
Calculates the rolling mean of the infections panel data for all time points for each region.
Set the rolling window by the parameter \code{k} (default: 7). 
Set the fill value for the observations left/within/right to the data range with parameter \code{fill} (default: NA).
Parameter \code{align} defines whether the index of the result should be left- or right-aligned or centered (default).
If \code{col_name} is \code{NULL}, the column is defined as "<Column name of cases>_rm".
Set \code{overwrite} to TRUE, if an existing column should be overwritten.
The method uses the function \code{rollmean} from the \code{zoo} package (see the corresponding documentation).
}

\value{
\code{infpan} object including column with rolling means in the infections panel data
}

\references{
Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{infpan-class}}, \code{\link{load_infections_paneldata}}, \code{\link{calculate_rollsum}},
\code{\link{calculate_Rt}}, \code{\link{calculate_cum}}, \code{\link{calculate_incidence}}
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

infpan_CH <- calculate_rollmean(
  infpan_CH, 
  col_name = "RollingMean",
  verbose = TRUE
)
# Calculate rolling mean of cases as "RollingMean"

summary(infpan_CH)
# Summary of infpan object
}