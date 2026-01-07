\encoding{UTF-8}
\name{Infections}
\docType{data}
\alias{Infections}

\title{
Infections
}

\description{
Daily SARS-CoV-2 infection in Germany spring 2020
}

\usage{
data(Infections)
}

\format{
A data.frame with multiple columns:
\describe{
  \item{infectedtest_CW2}{Calendar week 2020 of conducted test}
  \item{infection_date}{Estimated date of infection}
  \item{infections_daily}{Daily infections}
  \item{infections_daily_lwr}{Daily infections lower confidence interval}
  \item{infections_daily_upr}{Daily infections upper confidence interval}
  \item{infections_cum}{Cumulative infections}
  \item{infections_cum_lwr}{Cumulative infections lower confidence interval}
  \item{infections_cum_upr}{Cumulative infections upper confidence interval}
  \item{R4}{Estimated effective reproduction number \eqn{R_t} with generation interval = 4}
  \item{R4_lwr}{Estimated effective reproduction number \eqn{R_t} with generation interval = 4 lower confidence interval}
  \item{R4_upr}{Estimated effective reproduction number \eqn{R_t} with generation interval = 4 upper confidence interval}
  \item{R7}{Estimated effective reproduction number \eqn{R_t} with generation interval = 7}
  \item{R7_lwr}{Estimated effective reproduction number \eqn{R_t} with generation interval = 7 lower confidence interval}
  \item{R7_upr}{Estimated effective reproduction number \eqn{R_t} with generation interval = 7 upper confidence interval}
  \item{onsets_of_symptoms}{Daily onsets of symptoms}
  \item{onsets_of_symptoms_lwr}{Daily onsets of symptoms lower confidence interval}
  \item{onsets_of_symptoms_upr}{Daily onsets of symptoms upper confidence interval}
  \item{reported_cases}{Daily reported cases}
  \item{day}{Time counter (day)}
  \item{ln_inf_cum}{Nat. log. of cumulative infections}
  \item{ln_inf_daily}{Nat. log. of daily infections}
  \item{ln_R4}{Nat. log. of estimated effective reproduction number \eqn{R_t} with generation interval = 4}
  \item{ln_R7}{Nat. log. of estimated effective reproduction number \eqn{R_t} with generation interval = 7}
  \item{infection_date_CW}{Calendar week of infection data (numeric)}
  \item{infection_date_CW2}{Calendar week of infection data (categorical)}
  \item{infectedtest_CW}{Calendar week of conducted test}
  \item{conducted_tests}{No. of conducted tests}
  \item{negative_tests}{No. of negative tests}
  \item{positive_tests}{No. of positive tests}
  \item{positive_tests_share}{Share of positive tests (average per day)}
  \item{conducted_tests_index}{No. of conducted tests (average per day), index (CW 14 = 100)}
  \item{conducted_tests_dailyaverage}{No. of conducted tests, average per day}
  \item{positive_tests_dailyaverage}{Positive tests, average per day}
  \item{infections_daily_testweighted}{Daily infections weighted by test volume}
  \item{ln_inf_daily_tw}{Nat. log. of daily infections weighted by test volume}
}
}

\details{
Example data with daily SARS-CoV-2 infections in Germany. See Wieland (2020) for data sources and method of backdating infections.
}

\source{
Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}
}

\examples{
data(Infections)
}