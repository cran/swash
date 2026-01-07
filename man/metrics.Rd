\encoding{UTF-8}
\name{metrics}
\alias{metrics}

\title{
Fit metrics of observed and expected numeric variables
}
\description{
Calculation of fit metrics for observed and expected numeric variables 
(e.g. \eqn{R^2}, \eqn{RMSE}, \eqn{MAE}, \eqn{MAPE}).
}
\usage{
metrics(
  observed,
  expected,
  plot = TRUE,
  plot.main = "Observed vs. expected",
  xlab = "Observed",
  ylab = "Expected",
  point.col = "blue",
  point.pch = 19,
  line.col = "red",
  plot_residuals.main = "Residuals",
  legend.cex = 0.7
)
}

\arguments{
  \item{observed}{
     Numeric vector of observed values.
  }
  \item{expected}{
     Numeric vector of expected or predicted values.
  }
  \item{plot}{
     Logical. If \code{TRUE}, diagnostic plots for observed vs. expected values
     and relative residual distributions are created.
  }
  \item{plot.main}{
     Character string. Title of the observed vs. expected plot.
  }
  \item{xlab}{
     Character string. Label of the x-axis.
  }
  \item{ylab}{
     Character string. Label of the y-axis.
  }
  \item{point.col}{
     Color of points in the observed vs. expected plot.
  }
  \item{point.pch}{
     Plotting character used for points.
  }
  \item{line.col}{
     Color of the identity line (\eqn{y = x}).
  }
  \item{plot_residuals.main}{
     Character string. Title of the residuals bar plot.
  }
  \item{legend.cex}{
     Numeric. Character expansion factor for legends.
  }
}

\details{
The function computes several goodness-of-fit metrics comparing observed and
expected numeric values. In addition to classical error measures such as
mean squared error (MSE), root mean squared error (RMSE), mean absolute error (MAE),
and mean absolute percentage error (MAPE), the coefficient of determination
(\eqn{R^2}) is calculated.

If \code{plot = TRUE}, the function produces:
\itemize{
  \item a scatter plot of observed versus expected values including the identity line,
  \item a bar plot of relative residual frequencies.
}
}

\value{
A \code{list} with two elements:
\item{\code{fit_metrics}}{
  A \code{list} containing the computed fit metrics:
  \code{SQR}, \code{SAR}, \code{SQT}, \code{R2}, \code{MSE}, \code{RMSE},
  \code{MAE}, and \code{MAPE}.
}
\item{\code{observed_expected}}{
  A \code{data.frame} containing observed values, expected values,
  residuals, and derived residual measures.
}
}

\references{
Boehmke B, Greenwell B (2020). \emph{Hands-On Machine Learning with R} (1st ed.). Taylor & Francis, New York, NY.
}

\author{
Thomas Wieland
}

\seealso{
\code{\link{metrics}}
}

\examples{
obs <- c(10, 12, 15, 18, 20)
exp <- c(11, 13, 14, 17, 21)

metrics(
  observed = obs,
  expected = exp
)
}