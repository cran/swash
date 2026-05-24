\encoding{UTF-8}
\name{infpan-class}
\docType{class}
\alias{infpan-class}
\alias{plot,infpan-method}
\alias{print,infpan-method}
\alias{show,infpan-method}
\alias{summary,infpan-method}
\alias{swash,infpan-method}

\title{Class \code{"infpan"}}

\description{
The class \code{"infpan"} contains infections panel data for \eqn{N} regions and \eqn{T} time points as well as additional information. 
Use \code{summary(infpan)} and \code{plot(infpan)} for results summary and plotting, respectively.
}

\section{Objects from the Class}{
Objects can be created by importing infections panel data using the function \code{load_infections_paneldata}.
}

\section{Slots}{
  \describe{
    \item{\code{input_data}:}{Object of class \code{"data.frame"} Model result: Input infections panel data}
    \item{\code{data_statistics}:}{Object of class \code{"numeric"} Data statistics (\eqn{N} regions, \eqn{T} time points, test whether data is balanced, etc.)}
    \item{\code{index_col_names}:}{Object of class \code{"character"} Column names of regions and time points }
    \item{\code{cases_col_name}:}{Object of class \code{"character"} Column name of incremental cases }
    \item{\code{other_cols}:}{Object of class \code{"character"} Names of other relevant columns derived from incremental case data, e.g. effective reproduction number \eqn{R_t} }
    \item{\code{time_format}:}{Object of class \code{"character"} Format of time points in time column }
    \item{\code{time_unit}:}{Object of class \code{"character"} Time unit, default: "days" }
    \item{\code{timestamp}:}{Object of class \code{"list"} Time stamps of any update of the instance }
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "infpan")}: Plots case data by region for \eqn{N} regions and \eqn{T} time points }
    \item{calculate_Rt}{\code{signature(x = "infpan")}: Calculates the effective reproduction number \eqn{R_t} from \code{infpan} objects. Returns updated \code{infpan} instance.}
    \item{calculate_cum}{\code{signature(x = "infpan")}: Calculates cumulative cases from \code{infpan} objects. Returns updated \code{infpan} instance.}
    \item{calculate_rollmean}{\code{signature(x = "infpan")}: Calculates rolling means of cases from \code{infpan} objects. Returns updated \code{infpan} instance.}
    \item{calculate_rollsum}{\code{signature(x = "infpan")}: Calculates rolling sums of cases from \code{infpan} objects. Returns updated \code{infpan} instance.}
    \item{calculate_incidence}{\code{signature(x = "infpan")}: Calculates incidences of cases from \code{infpan} objects. Returns updated \code{infpan} instance.}
    \item{print}{\code{signature(x = "infpan")}: Prints an \code{infpan} object; use \code{summary(infpan)} for results }
    \item{show}{\code{signature(object = "infpan")}: Shows an \code{infpan} object; use \code{summary(infpan)} for results }
    \item{summary}{\code{signature(object = "infpan")}: Prints a summary of \code{infpan} objects }
    \item{swash}{\code{signature(object = "infpan")}: Performs a Swash-Backwash Model analysis from \code{infpan} objects. Returns \code{sbm} instance. }
    \item{growth}{\code{signature(object = "infpan")}: Estimates logistic growth models from \code{infpan} objects. Returns \code{growthmodels} instance. }
    \item{growth_initial}{\code{signature(object = "infpan")}: Estimates exponential growth models from \code{infpan} objects for a given time period. Returns \code{growthmodels} instance. }
    \item{growth_hawkes}{\code{signature(object = "infpan")}: Estimates Hawkes process models from \code{infpan} objects. Returns \code{growthmodels} instance. }
    \item{breaks_growth}{\code{signature(object = "infpan")}: Estimates breakpoints in linear regression models from \code{infpan} objects. Returns \code{growthmodels} instance. }
	 }
}

\author{
Thomas Wieland
}

\references{
Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
}

\examples{
showClass("infpan")
}

\keyword{classes}