\encoding{UTF-8}
\name{compare_countries}
\alias{compare_countries}

\title{
Two-country Comparison of Swash-Backwash Model Parameters
}
\description{
This function enables bootstrap estimates for the mean difference of Swash-Backwash Model parameters of two countries to be compared.
}
\usage{
compare_countries(
  sbm1, 
  sbm2, 
  country_names = c("Country 1", "Country 2"), 
  indicator = "R_0A", 
  iterations = 20, 
  samples_ratio = 0.8, 
  alpha = 0.05, 
  replace = TRUE
  )
}
\arguments{
  \item{sbm1}{
     A \code{sbm} object for country 1
}
  \item{sbm2}{
     A \code{sbm} object for country 2
}
  \item{country_names}{
     \code{list} with user-given country names (two entries)
}
  \item{indicator}{
     \code{character}, indicator to be analyzed (\code{"S_A"}, \code{"I_A"}, \code{"R_A"}, \code{"t_LE"}, \code{"t_LE"}, or \code{"R_0A"} (default and recommended: \code{"R_0A"}))
}
  \item{iterations}{
     Number of iterations for resampling (default: 100)
}
  \item{samples_ratio}{
     Proportion of regions included in each sample (default: 0.8)
}
  \item{alpha}{
     Significance level \eqn{\alpha} for the confidence intervals (default: 0.05)
}
  \item{replace}{
     Resampling with replacement (\code{TRUE} or \code{FALSE}, default: \code{TRUE} = bootstrap resampling)
}
}
\details{
The combination of the Swash-Backwash Model and bootstrap resampling allows the estimation of mean differences of a user-specified model parameter (e.g., spatial reproduction number \eqn{R_{OA}}) between two countries. This makes it possible to check whether the spatial spread velocity of a communicable disease is significantly different in one country than in another country. Since the initial data in the Swash-Backwash Model should be balanced, \emph{entity-based bootstrap sampling} is carried out in the \code{compare_countries()} function. This means that not, for example, 80\% of all observations are included in each sample at a sample ratio equal to \eqn{p} = 0.8, but rather all observations for 80\% of the regions. For both countries, \eqn{B} bootstrap samples (default: 100) are drawn for which the Swash-Backwash Model is calculated. Based on the distribution of indicators, confidence intervals are calculated at the user-specified significance level \eqn{\alpha}. The \code{compare_countries()} function calculates the differences of the user's desired indicator between the two samples, \eqn{D}, and also calculates \eqn{\alpha} confidence intervals for this.
}
\value{
object of class \code{countries}, see \code{\link{countries-class}}
}
\references{

\emph{Swash-Backwash Model:}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

\emph{Bootstrapping und bootstrap confidence intervals:}

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap.

Ramachandran KM, Tsokos CP (2021) Mathematical Statistics with Applications in R (Third Edition). Ch. 13.3.1 (Bootstrap confidence intervals). \doi{https://doi.org/10.1016/B978-0-12-817815-7.00013-0}
}
\author{
Thomas Wieland
}
\seealso{
\code{\link{swash}}, \code{\link{countries-class}}
}
\examples{
data(COVID19Cases_geoRegion)
# Get Swiss COVID19 cases at NUTS 3 level

data(Oesterreich_Faelle)
# Get Austrian COVID19 cases at NUTS 3 level
# (first wave, same final date as in Swiss data: 2020-05-31)

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

AT_covidwave1 <- 
  swash (
    data = Oesterreich_Faelle,
    col_cases = "Faelle",
    col_date = "Datum",
    col_region = "NUTS3"
  )
# Swash-Backwash Model for Austrian COVID19 cases
# Spatial aggregate: NUTS 3

AT_vs_CH <- 
  compare_countries(
    CH_covidwave1, 
    AT_covidwave1,
    country_names = c("Switzerland", "Austria"))
# Country comparison Switzerland vs. Austria
# default config: 20 iterations, alpha = 0.05, sample ratio = 80\%,
# indicator: R_0A

summary(AT_vs_CH)
# Summary of country comparison

plot(AT_vs_CH)
# Plot of country comparison
}