\name{plot_regions-methods}
\docType{methods}
\alias{plot_regions-methods}
\alias{plot_regions,sbm-method}
\title{Methods for Function \code{plot_regions}}
\description{
 Methods for function \code{plot_regions}
}
\section{Methods}{
\describe{

\item{\code{signature(
    object = "sbm", 
    col = "red",
    lty = "solid",
    scale = FALSE, 
    normalize_by_col = NULL, 
    normalize_factor = 1, 
    plot_rollmean = FALSE,
    rollmean_col = "blue",
    rollmean_lty = "solid",
    rollmean_k = 7,
    rollmean_align = "center",
    rollmean_fill = NA,
    plot_growth = FALSE,
    growth_rollmean = FALSE,
    growth_col = "orange",
    growth_lty = "solid",
    growth_per_time_unit = 1)}}{
Plots regional infection curves in \eqn{N} plots for \eqn{N} regions, with the number of columns equals 4 and the number of rows is calculated based the size of \eqn{N}. Set the color by the argument \code{col}. If \code{scale} is \code{TRUE}, the y axis of the plots is set uniformly based on the maximum value across all regions. If the input data contains a column to normalize the infection numbers (such as regional population), the user may use this data to normalize the infection numbers by setting \code{normalize_by_col}. If this density value should be multiplied by a factor (e.g. regional infections per 100,000 inhabitants), this can be set with the argument \code{normalize_by_col} (default: 1).
}
}}
\author{
Thomas Wieland
}
\keyword{methods}