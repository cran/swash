\encoding{UTF-8}
\name{did_fatalities_splm_coef}
\docType{data}
\alias{did_fatalities_splm_coef}

\title{
Results from a Difference-in-Differences Model
}

\description{
Example data frame with results from a difference-in-differences model
}

\usage{
data(did_fatalities_splm_coef)
}

\format{
A \code{data.frame} with multiple columns:
\describe{
  \item{Var}{Coefficient name}
  \item{Estimate}{Coef. estimate}
  \item{Std_Error_Bonferroni}{Coef. standard error}
  \item{t_value_Bonferroni}{Coef. t value}
  \item{Pr_t_Bonferroni}{Coef. p value}
  \item{CI_lower_Bonferroni}{Coef. lower confidence interval}
  \item{CI_upper_Bonferroni}{Coef. upper confidence interval}
}
}

\details{
Data frame with results from a difference-in-differences model (SPLM model), example data
}

\source{
Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
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