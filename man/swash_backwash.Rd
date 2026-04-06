\encoding{UTF-8}
\name{swash_backwash}
\alias{swash_backwash}

\title{
Swash-Backwash Model for the Single Epidemic Wave
}

\description{
Analysis of regional infection/surveillance data using the Swash-Backwash Model for the Single Epidemic Wave by Cliff and Haggett (2006).
}

\usage{
swash_backwash(
  infpan = NULL,
  data = NULL,
  col_cases = NULL, 
  col_date = NULL, 
  col_region = NULL,
  time_format = "\%Y-\%m-\%d",
  verbose = FALSE
  )
}

\arguments{
  \item{infpan}{
     \code{infpan} object containing regional infection data
}  \item{data}{
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
  \item{time_format}{
     \code{character}, Time format of the values in \code{col_date}
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}
\details{
The function performs the analysis of the input panel data with \eqn{N} regions and \eqn{T} time points using the Swash-Backwash Model. 
The user must state panel data with daily infections.

The Swash-Backwash Model (SBM) for the Single Epidemic Wave is the spatial equivalent of the classic epidemiological SIR (Susceptible-Infected-Recovered) model. 
It was developed by Cliff and Haggett (2006) to model the velocity of spread of infectious diseases across space. 
Current applications can be found, for example, in Smallman-Raynor et al. (2022a,b). 
The function  \code{swash_backwash()} enables the calculation of the Swash-Backwash Model for user-supplied panel data on regional infections. 
It calculates the model and creates a model object of the \code{sbm} class defined in this package. 
This class can be used to visualize results (\code{summary()}, \code{plot()}) and calculate bootstrap confidence intervals for the model estimates (\code{confint(sbm)}); 
the latter returns an object of class \code{sbm_ci} as defined in this package. Two \code{sbm_ci} objects for different countries may be compared with \code{compare_countries()}, 
which allows the estimation of mean differences of a user-specified model parameter (e.g., spatial reproduction number \eqn{R_{OA}}) between two countries. 
This makes it possible to check whether the spatial spread velocity of a communicable disease is significantly different in one country than in another country; 
the result is an object of class \code{countries}. 

To calculate the SBM model based on an \code{infpan} object, use the corresponding method \code{swash(infpan)}.
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
\code{\link{sbm-class}, \link{swash}}
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
  swash_backwash(
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