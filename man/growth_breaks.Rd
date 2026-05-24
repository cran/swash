\encoding{UTF-8}
\name{growth_breaks}
\alias{growth_breaks}

\title{
Time Series Model with Breakpoints for Regional Infections
}

\description{
Conducts \eqn{N} breakpoints analyses for infection time series in \eqn{N} regions.
}

\usage{
growth_breaks(
  object,
  ln = FALSE,
  add_constant = 1,
  alpha = 0.05,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{ln}{
     \code{bool} argument which indicates whether the dependent variable should be transformed by natural logarithm
}
  \item{add_constant}{
     Numeric constant to be added to y if zero values occur
}
  \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
  \item{verbose}{
     \code{bool} argument which indicates whether progress messages are displayed
}
}

\details{
The method detects breakpoints in regional infections time series based on an \code{infpan} object. 
The function internally uses the function \code{breakpoints} from the \code{strucchange} package (Zeileis et al. 2003),
where breakpoints are identified using the Bai-Perron algorithm (Bai & Perron 2003).
See \code{\link{breaks_growth}} for further details of the estimation.
}

\value{
object of class \code{\link{growthmodels-class}}
}

\references{
Bai J, Perron P (2003) Computation and analysis of multiple structural change models. \emph{Journal of Applied Econometrics} 18(1), 1-22. \doi{https://doi.org/10.1002/jae.659}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Zeileis C, Kleiber W, Krämer K, Hornik, K (2003) Testing and dating of structural changes in practice. \emph{Computational Statistics & Data Analysis} 44(1-2), 109-123. \doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
}

\author{
Thomas Wieland
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

CH_covidwave1_breaks <- 
  growth_breaks(infpan_CH)
summary(CH_covidwave1_breaks)
# Breakpoints for infpan object infpan_CH 
}