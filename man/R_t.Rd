\encoding{UTF-8}
\name{R_t}
\alias{R_t}

\title{
Effective Reproduction Number for Epidemic Data
}
\description{
Calculation of the effective reproduction number for infection/surveillance data
}
\usage{
R_t(
  infections, 
  GP = 4,
  correction = FALSE
  )
}

\arguments{
  \item{infections}{
     \code{numeric vector} with infection data
}
  \item{GP}{
     Generation period, in days
}
  \item{correction}{
     Correction of values equal to zero? (Recommended)
}
}
\details{
The function calculates the effective reproduction number (=growth factor), \eqn{R_t}, of an infections time series.
}
\value{
\code{list} with two entries:
    \item{\code{R_t}:}{Object of class \code{"numeric"} \eqn{R_t} values}
    \item{\code{infections_data}:}{Object of class \code{"data.frame"} Dataset with infections data and \eqn{R_t} }
}
\references{
Bonifazi G et al. (2021) A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. \emph{The European Physical Journal Plus} 136, 386. \doi{https://doi.org/10.1140/epjp/s13360-021-01339-6}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
}
\author{
Thomas Wieland
}
\seealso{
\code{\link{logistic_growth}}
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

COVID19Cases_BS <-
  COVID19Cases_geoRegion[(COVID19Cases_geoRegion$geoRegion == "ZH")
                         & (COVID19Cases_geoRegion$sumTotal > 0),]
# COVID cases for Zurich

Rt_BS <- R_t(infections = COVID19Cases_BS$entries)
# Effective reproduction number

Rt_BS
}