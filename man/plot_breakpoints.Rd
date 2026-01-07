\encoding{UTF-8}
\name{plot_breakpoints}
\alias{plot_breakpoints}

\title{
Calculating and plotting breakpoints in time series
}
\description{
The function calculates breakpoints for time series and provides plots of the analysis results
}
\usage{
plot_breakpoints(
  dataset, 
    formula,
    alpha = 0.05,
    line.col = "chocolate",
    ci.col = c(
      rgb(0, 0, 255, maxColorValue = 255, alpha = 0.5),
      rgb(0, 255, 0, maxColorValue = 255, alpha = 0.5),
      rgb(255, 0, 0, maxColorValue = 255, alpha = 0.5)
    ),
    legend.show = TRUE, 
    legend.pos = c(
      "topright", 
      "topright", 
      "topright", 
      "topright"
      ),
    xlab = "Time", 
    ylab = "Y", 
    ylim = NULL,
    plot.main = c(
      "Breakpoints", 
      "Segment model fits", 
      "BIC and Residual Sum of Squares", 
      "Model without breaks (one segment)"
    ),
    output.full = TRUE,
    separate_plots = FALSE,
    ...
  )
}

\arguments{
  \item{dataset}{
     \code{data.frame} object with (infections) time series
}
  \item{formula}{
     \code{formula} object for linear regression formula
}
  \item{alpha}{
     Significance level \eqn{\alpha} for 1-\eqn{\alpha}*100 confidence intervals
}
  \item{line.col}{
     Line color in time series plot
}
  \item{ci.col}{
     Colors for confidence intervals in 3 plots
}
  \item{legend.show}{
     \code{bool}: Show legend?
}
  \item{legend.pos}{
     Position of legend
}
  \item{xlab}{
     Label of x axis
}
  \item{ylab}{
     Label of y axis
}
  \item{ylim}{
     Limits of y axis
}
  \item{plot.main}{
     Titles for 4 plots
}
  \item{output.full}{
     \code{bool}: All 4 plots or just one (breakpoints)?
}
  \item{separate_plots}{
     \code{bool}: All 4 plots in one or separate?
}
  \item{...}{
     Further arguments passed to breakpoints()
}
}

\details{
The function uses the breakpoints() function from the \emph{strucchange} package (Zeileis et al. 2003) to estimate breaks in time series given by the user. These breakpoints are visualized, confidence intervals and model diagnostics are calculated, and for \eqn{M} breakpoints, \eqn{M+1} models are estimated and returned.
}

\value{
\code{data.frame} with one row for each model segment with coefficients, \eqn{R^2}, etc.
}

\references{
Wieland T (2020) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Zeileis A, Kleiber C, Krämer W, Hornik K (2003) Testing and dating of structural changes in practice. \emph{Computational Statistics & Data Analysis} 44(1-2) 109–123. \doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{logistic_growth}, \link{exponential_growth}}
}

\examples{
data(Infections)
# Confirmed SARS-CoV-2 cases in Germany

plot_breakpoints(
  Infections, 
  log(infections_daily) ~ day,
  output.full = TRUE
  )
# Breakpoints for time series
}