\encoding{UTF-8}
\name{as_balanced}
\alias{as_balanced}
\title{
Correction of Non-balanced Panel Dataset with Regional Infection Data
}
\description{
This function corrects non-balanced input panel data by replacing missing entries with a user-given constant (e.g., 0).
}
\usage{
as_balanced(
  data, 
  col_cases, 
  col_date, 
  col_region, 
  fill_missing = 0
  )
}
\arguments{
  \item{data}{
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
  \item{fill_missing}{
     Constant to fill missing values (default and recommended: 0)
}
}
\details{
The Swash-Backwash Model for the Single Epidemic Wave does not necessarily require balanced panel data in order for the calculations to be carried out. However, for a correct estimation it is implicitly assumed that the input data is balanced. The function corrects non-balanced panel data. It is executed automatically whithin the \code{swash()} function (when using the function \code{is_balanced()}), but can also be used separately.
}
\value{
  \item{data }{Corrected input dataset (\code{data.frame})}
}
\references{

\emph{Swash-Backwash Model:}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

\emph{Panel data:}

Greene, WH (2012) \emph{Econometric Analysis}. Ch. 11.

Wooldridge, JM (2012) \emph{Introductory Econometrics. A Modern Approach}. Ch. 13.
}
\author{
Thomas Wieland
}
\seealso{
\code{\link{is_balanced}}
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

COVID19Cases_geoRegion_balanced <- 
  is_balanced(
  data = COVID19Cases_geoRegion,
  col_cases = "entries",
  col_date = "datum",
  col_region = "geoRegion"
)
# Test whether "COVID19Cases_geoRegion" is balanced panel data 

COVID19Cases_geoRegion_balanced$data_balanced
# Balanced? TRUE or FALSE

if (COVID19Cases_geoRegion_balanced$data_balanced == FALSE) {
  COVID19Cases_geoRegion <- 
    as_balanced(
    COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion"
  )
}
# Correction of dataset "COVID19Cases_geoRegion"
# not necessary as parameter balance of is_balanced is set TRUE by default
}