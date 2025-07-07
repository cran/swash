\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot,sbm,ANY-method}
\alias{plot,sbm_ci,ANY-method}
\alias{plot,loggrowth,ANY-method}
\alias{plot,countries,ANY-method}
\title{Methods for Function \code{plot}}

\description{
Methods for function \code{plot} for different S4 classes: \code{sbm}, \code{sbm_ci}, and \code{loggrowth}.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "sbm")}}{
  \code{plot.sbm(x, y, ...)}:
  Plots the results of the Swash-Backwash Model. This generates two plots:
  \itemize{
    \item Edges over time.
    \item Total infections per time unit.
  }
  
  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{sbm} representing the results of the Swash-Backwash Model.
    \item \code{y}: Optional argument for additional customization, such as plot style or axis labels.
    \item \code{...}: Additional graphical parameters that can be passed to control plot appearance.
  }

  \strong{Details}: This method is used to visualize the output of the Swash-Backwash Model, providing insight into the dynamics of the modeled epidemic.
}

\item{\code{signature(x = "sbm_ci")}}{
  \code{plot.sbm_ci(x, y, ...)}:
  Plots the results of bootstrap confidence intervals for the Swash-Backwash Model. This generates a single figure with six subplots:
  \itemize{
    \item \eqn{S_A} (susceptible population),
    \item \eqn{I_A} (infected population),
    \item \eqn{R_A} (recovered population),
    \item \eqn{t_{FE}} (final epidemic time),
    \item \eqn{t_{LE}} (last epidemic time),
    \item \eqn{R_{0A}} (basic reproduction number).
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{sbm_ci} containing the bootstrap confidence intervals for the Swash-Backwash Model.
    \item \code{y}: Optional argument for additional customization, such as plot style or axis labels.
    \item \code{...}: Additional graphical parameters for fine-tuning the plots.
  }

  \strong{Details}: This method is used to visualize the bootstrap confidence intervals for various parameters of the Swash-Backwash Model.
}

\item{\code{signature(x = "countries")}}{
  \code{plot.sbm(x, y = NULL, col_bars = "grey", col_ci = "red")}:
  Plots the results of the between-countries analysis via Swash-Backwash Model. This generates four plots:
  \itemize{
    \item Indicator for country 1
    \item Indicator for country 2
    \item Boxplots of the distribution of the indicator in country 1 and 2
    \item Distribution of the difference between the indicators of country 1 and 2
  }
  
  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{countries} representing the results of the Swash-Backwash Model country analysis.
    \item \code{y}: Not relevant
    \item \code{col_bars}: Color of bars
    \item \code{col_ci}: Color of confidence intervals
  }

  \strong{Details}: This method is used to visualize the output of the Swash-Backwash Model, providing insight into the dynamics of the modeled epidemic.
}

\item{\code{signature(x = "loggrowth")}}{
  \code{plot.loggrowth(x, y, ...)}:
  Plots the results of the logistic growth model, including:
  \itemize{
    \item Observed values,
    \item Predicted values,
    \item First derivative (growth rate).
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{loggrowth} containing the data for the logistic growth model.
    \item \code{y}: Optional argument for additional customization of the plot (e.g., color, labels).
    \item \code{...}: Additional arguments for graphical parameters.
  }

  \strong{Details}: This method is useful for visualizing the observed and predicted growth patterns in an epidemic or similar phenomena modeled by logistic growth.
}

}}
\author{
Thomas Wieland
}
\keyword{methods}