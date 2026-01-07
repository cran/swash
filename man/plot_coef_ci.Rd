\encoding{UTF-8}
\name{plot_coef_ci}
\alias{plot_coef_ci}

\title{
Plot Point Estimates With Confidence Intervals
}
\description{
Plotting point estimates with confidence intervals from regression results
}
\usage{
plot_coef_ci(
  point_estimates,
  confint_lower,
  confint_upper,
  coef_names,
  p = NULL,
  estimate_colors = NULL,
  confint_colors = NULL,
  auto_color = FALSE,
  alpha = 0.05,
  set_estimate_colors = c("red", "grey", "green"),
  set_confint_colors = c("#ffcccb", "lightgray", "#CCFFCC"),
  skipvars = NULL,
  plot.xlab = "Independent variables",
  plot.main = "Point estimates with CI",
  axis.at = seq(-30, 40, by = 5),
  pch = 15,
  cex = 2,
  lwd = 5,
  y.cex = 0.8
  )
}

\arguments{
  \item{point_estimates}{
     \code{numeric} vector containing point estimates
}
  \item{confint_lower}{
     \code{numeric} vector containing lower confidence intervals
}
   \item{confint_upper}{
     \code{numeric} vector containing upper confidence intervals
}
   \item{coef_names}{
     \code{character} vector containing coefficient names
}
   \item{p}{
     \code{numeric} vector containing \eqn{p} values of the coefficients (optional)
}
   \item{estimate_colors}{
     vector containing colors for the point estimates (optional)
}
   \item{confint_colors}{
     vector containing colors for the confidence intervals (optional)
}
   \item{auto_color}{
     \code{bool} value which indicates whether the colors are found automatically based on coef and CI values
}
   \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
   \item{set_estimate_colors}{
     Colors for point estimates (significant negative, not significant, significant positive)
}
   \item{set_confint_colors}{
     Colors for confidence intervals (significant negative, not significant, significant positive)
}
   \item{skipvars}{
     List with coefficients to be dropped
}
   \item{plot.xlab}{
     Label of x axis
}
   \item{plot.main}{
     Plot title
}
   \item{axis.at}{
     Position of y axis
}
   \item{pch}{
     Point type
}
   \item{cex}{
     Point size
}
   \item{lwd}{
     Line width (confidence intervals)
}
   \item{y.cex}{
     Font size of y axis
}
}

\details{
The function checks whether the input vectors have the same length. If \code{auto_color} is \code{TRUE}, the colors from \code{set_estimate_colors} and \code{set_confint_colors} are used, and the significance level is determined based on the coefficient and confidence interval values (all three below 0 = significant negative, all three above 0 = significant positive).
}

\value{
Coefficients plot, no returned value
}

\references{
Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
}

\author{
Thomas Wieland
}

\examples{
data(did_fatalities_splm_coef)
# Results of a difference-in-differences model

plot_coef_ci(
  point_estimates = did_fatalities_splm_coef$Estimate,
  confint_lower = did_fatalities_splm_coef$CI_lower_Bonferroni,
  confint_upper = did_fatalities_splm_coef$CI_upper_Bonferroni,
  coef_names = did_fatalities_splm_coef$Var,
  skipvars = c(
    "Alpha_share", 
    "lambda",
    "rho",
    "log(D_Infections_daily_7dsum_per100000_lag2weeks)",
    "vacc_cum_per100000_lag2weeks"
    ),
  lwd = 13,
  pch = 19,
  auto_color = TRUE
)
# Plot with point estimates and confidence intervals
}