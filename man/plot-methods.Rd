\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot,infpan,ANY-method}
\alias{plot,sbm,ANY-method}
\alias{plot,sbm_ci,ANY-method}
\alias{plot,loggrowth,ANY-method}
\alias{plot,expgrowth,ANY-method}
\alias{plot,countries,ANY-method}
\alias{plot,breaksgrowth,ANY-method}
\alias{plot,hawkes,ANY-method}
\title{Methods for Function \code{plot}}

\description{
Methods for function \code{plot} for different S4 classes: 
\code{infpan}, \code{sbm}, \code{sbm_ci}, \code{loggrowth}, \code{expgrowth}, 
\code{hawkes}, and \code{breaksgrowth}.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "infpan")}}{
  \code{plot.infpan(x, y, ...)}:
  Plots regional infections against time
  
  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{infpan} including infections panel data.
    \item \code{y}: Optional argument for additional customization, such as plot style or axis labels.
    \item \code{...}: Additional graphical parameters that can be passed to control plot appearance.
  }

  \strong{Details}: This method is used to visualize case data by region for \eqn{N} regions and \eqn{T} time points.
}

\item{\code{signature(x = "sbm")}}{
  \code{plot.sbm(x, 
    y = NULL, 
    col_edges = "blue",
    xlab_edges = "Time", 
    ylab_edges = "Regions",
    main_edges = "Edges",
    col_SIR = c("blue", "red", "green"),
    lty_SIR = c("solid", "solid", "solid"),
    lwd_SIR = c(1,1,1),
    xlab_SIR = "Time", 
    ylab_SIR = "Regions",
    main_SIR = "SIR integrals",
    col_cases = "red",
    lty_cases = "solid",
    lwd_cases = 1,
    xlab_cases = "Time", 
    ylab_cases = "Infections",
    main_cases = "Daily infections",
    xlab_cum = "Cases",
    ylab_cum = "Regions",
    main_cum = "Cumulative infections per region",
    horiz_cum = TRUE,
    separate_plots = FALSE)}:
    
  Plots the results of the Swash-Backwash Model. This generates four plots:
  \itemize{
    \item Edges over time.
    \item SIR integrals.
    \item All-over infections.
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
    \item Observed values
    \item Predicted values
    \item First derivative
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{loggrowth} containing the data for the logistic growth model.
    \item \code{y}: Optional argument for additional customization of the plot (e.g., color, labels).
    \item \code{...}: Additional arguments for graphical parameters.
  }

  \strong{Details}: This method is useful for visualizing the observed and predicted growth patterns in an epidemic or similar phenomena modeled by logistic growth.
}

\item{\code{signature(x = "expgrowth")}}{
  \code{plot.expgrowth(x, y, ...)}:
  Plots the results of the exponentai growth model, including:
  \itemize{
    \item Observed values
    \item Predicted values
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{expgrowth} containing the data for the exponential growth model.
    \item \code{y}: Optional argument for additional customization of the plot (e.g., color, labels).
    \item \code{...}: Additional arguments for graphical parameters.
  }

  \strong{Details}: This method is useful for visualizing the observed and predicted growth patterns in the initial phase of an epidemic or similar phenomena modeled by exponential growth.
}

\item{\code{signature(x = "hawkes")}}{
  \code{plot.hawkes(x, y, ...)}:
  Plots the results of the Hawkes process model, including:
  \itemize{
    \item Observed values
    \item Predicted values
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{hawkes} containing the data for the Hawkes model.
    \item \code{y}: Optional argument for additional customization of the plot (e.g., color, labels).
    \item \code{...}: Additional arguments for graphical parameters.
  }

  \strong{Details}: This method is useful for visualizing the observed 
  and predicted growth patterns of an epidemic or similar 
  phenomena modeled as Hawkes processes.
}

\item{\code{signature(x = "breaksgrowth")}}{
  \code{plot.breaksgrowth(x, y, ...)}:
  Plots the results of a breakpoint analysis, including:
  \itemize{
    \item Time series data
    \item Breakpoints
  }

  \strong{Arguments}:
  \itemize{
    \item \code{x}: An object of class \code{breaksgrowth} containing the data for the breakspoints model.
    \item \code{y}: Optional argument for additional customization of the plot (e.g., color, labels).
    \item \code{...}: Additional arguments for graphical parameters.
  }

  \strong{Details}: This method is useful for visualizing the derived breakpoints.
}

}}

\author{
Thomas Wieland
}

\keyword{methods}