\encoding{UTF-8}
\name{swash-package}
\alias{swash-package}
\docType{package}
\title{
Implementation of the Swash-Backwash Model for the Single Epidemic Wave and additional functions in R 
}
\description{
Swash-Backwash Model for the single epidemic wave (Cliff and Haggett 2006) with additional functions for bootstrap confidence intervals and data management
}
\details{
The Swash-Backwash Model for the Single Epidemic Wave is the spatial equivalent of the classic epidemiological SIR (Susceptible-Infected-Recovered) model. It was developed by Cliff and Haggett (2006) to model the velocity of spread of infectious diseases across space. Current applications can be found, for example, in Smallman-Raynor et al. (2022a,b). This package enables the calculation of the Swash-Backwash Model for user-supplied panel data on regional infections. The core of this is the \code{swash()} function, which calculates the model and creates a model object of the \code{sbm} class defined in this package. This class can be used to visualize results (\code{summary()}, \code{plot()}) and calculate bootstrap confidence intervals for the model estimates (\code{confint(sbm)}). The package also contains additional helper functions.
}
\author{
Thomas Wieland
}
\references{

\emph{Swash-Backwash Model:}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022a) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022b) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

\emph{Basics of epidemiological modeling:}

Li, MY (2018) An Introduction to Mathematical Modeling of Infectious Diseases. \doi{https://doi.org/10.1007/978-3-319-72122-4}

Nishiura H, Chowell G (2009) The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. In Chowell G, Hyman JM, Bettencourt LMA (eds.) Mathematical and statistical estimation approaches in epidemiology, 103–121. \doi{https://doi.org/10.1007/978-90-481-2313-1_5}

\emph{Spatio-temporal analysis and modeling of infectious diseases:}

Bourdin S, Jeanne L, Nadou F, Noiret G (2021) Does lockdown work? A spatial analysis of the spread and concentration of Covid-19 in Italy. Regional Studies, 55, 1182–1193. \doi{https://doi.org/10.1080/00343404.2021.1887471}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. \emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Viboud C, Bjørnstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT (2006) Synchrony, Waves, and Spatial Hierarchies in the Spread of Influenza. \emph{Science} 312,447-451. \doi{https://doi.org/10.1126/science.1125237}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

\emph{Panel data:}

Greene, WH (2012) \emph{Econometric Analysis}. Ch. 11.

Wooldridge, JM (2012) \emph{Introductory Econometrics. A Modern Approach}. Ch. 13.

\emph{Bootstrapping und bootstrap confidence intervals:}

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap.

Ramachandran KM, Tsokos CP (2021) Mathematical Statistics with Applications in R (Third Edition). Ch. 13.3.1 (Bootstrap confidence intervals). \doi{https://doi.org/10.1016/B978-0-12-817815-7.00013-0}
}

\keyword{package}

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