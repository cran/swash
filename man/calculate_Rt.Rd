\encoding{UTF-8}
\name{calculate_Rt}
\alias{calculate_Rt}

\title{
Effective Reproduction Number
}

\description{
Calculation of the effective reproduction number for infections panel data.
}

\usage{
calculate_Rt(
  object,
  GP = 4,
  correction = FALSE,
  col_name = NULL,
  overwrite = FALSE,
  verbose = FALSE
  )
}

\arguments{
  \item{object}{
     object of class \code{infpan}
}
  \item{GP}{
     Generation period, in time units (typically days)
}
  \item{correction}{
     Correction of values equal to zero? (Recommended)
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
Calculates the effective reproduction number \eqn{R_t} for all time points for each region in the infections panel data.
Set the generation period by the parameter \code{GP} (default: 4). 
If \code{correction} is TRUE, values equal to zero are increased by one.
The method uses the built-in function \code{R_t()}.
}

\value{
\code{infpan} object including \eqn{R_t} column in the infections panel data
}

\references{
an der Heiden M, Hamouda O (2020) Schätzung der aktuellen Entwicklung der SARS-CoV-2-Epidemie in Deutschland - Nowcasting. \emph{Epidemiologisches Bulletin} 17, 10-15. \doi{https://doi.org/10.25646/6692}

Bonifazi G, Lista L, Menasce D, Mezzetto M, Pedrini D, Spighi R, Zoccoli A (2021) A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. \emph{The European Physical Journal Plus} 136, 386. \doi{https://doi.org/10.1140/epjp/s13360-021-01339-6}

Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
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

infpan_CH <- calculate_Rt(
  infpan_CH,
  verbose = TRUE
  )
# Calculate effective reproduction number

summary(infpan_CH)
# Summary of infpan object
}