\encoding{UTF-8}
\name{swash}
\alias{swash}

\title{
Swash-Backwash Model for the Single Epidemic Wave
}

\description{
Analysis of regional infection/surveillance data stored in \code{infpan} object 
using the Swash-Backwash Model for the single epidemic wave by Cliff and Haggett (2006).
}

\usage{
swash(
  object, 
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The method performs the analysis of the input panel data with \eqn{N} regions and \eqn{T} time points
using the Swash-Backwash Model based on an \code{infpan} object.
The output is an object of class \code{sbm}. 
The results can be viewed using \code{summary(sbm)}.
The built-in function \code{swash_backwash} is used for the analysis. 
See \code{\link{swash_backwash}} for further details.
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
\code{\link{sbm-class}}, \code{\link{swash_backwash}}
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