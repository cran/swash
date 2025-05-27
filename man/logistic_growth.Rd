\encoding{UTF-8}
\name{logistic_growth}
\alias{logistic_growth}

\title{
Logistic Growth Model for Epidemic Data
}
\description{
Estimation of logistic growth models from cumulative infections data, linearized OLS and/or NLS
}
\usage{
logistic_growth(
  y, 
  t, 
  S = NULL,
  S_start = NULL, 
  S_end = NULL, 
  S_iterations = 10, 
  S_start_est_method = "bisect", 
  seq_by = 10,
  nls = TRUE
  )
}

\arguments{
  \item{y}{
     \code{numeric vector} with cumulative infections data over time
}
  \item{t}{
     \code{vector} of class \code{numeric} or \code{Date} with time points or dates
}
  \item{S}{
     Saturation value for the model
}
  \item{S_start}{
     Start value of the saturation value for estimation
}
  \item{S_end}{
     End value of the saturation value for estimation
}
  \item{S_iterations}{
     Number of iterations for saturation value search
}
  \item{S_start_est_method}{
     Method for saturation value search, either "bisect" or "trial_and_error"
}
 \item{seq_by}{
     No of segments for the "trial_and_error" estimation of the saturation value
}
 \item{nls}{
     Nonlinear estimation? \code{TRUE} or \code{FALSE}
}
}
\details{
This function allows the estimation of a logistic growth model. The user must specify the dependent variable (cumulative infections) and the time variable (time counter or date values). The estimation is performed using a linearized model as an OLS estimator and as an NLS estimator. For the former, the saturation value can either be specified by the user or found using a search algorithm. The parameters from the OLS fit are used as starting values for the NLS estimation.
}
\value{
object of class \code{\link{loggrowth-class}}
}
\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}
}
\author{
Thomas Wieland
}
\seealso{
\code{\link{loggrowth-class}}
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

loggrowth_BS <- logistic_growth (
  y = as.numeric(COVID19Cases_BS$sumTotal), 
  t = COVID19Cases_BS$datum, 
  S = 5557,
  S_start = NULL, 
  S_end = NULL, 
  S_iterations = 10, 
  S_start_est_method = "bisect", 
  seq_by = 10,
  nls = TRUE
)
# Logistic growth model with stated saturation value

summary(loggrowth_BS)
# Summary of logistic growth model

plot(loggrowth_BS)
# Plot of logistic growth model
}