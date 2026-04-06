\encoding{UTF-8}
\name{growthmodels-class}
\docType{class}
\alias{growthmodels-class}
\alias{print,growthmodels-method}
\alias{show,growthmodels-method}
\alias{summary,growthmodels-method}

\title{Class \code{"growthmodels"}}

\description{
The class \code{"growthmodels"} contains the results of growth model analyses and the related input data as well as additional information.
The \code{swash} package includes the following model analyses under the heading "growth models": Exponential growth models,
logistic growth models, Hansen Process models, and time series models with breakpoints.
Use \code{summary(growthmodels)} for results summary.
See the corresponding functions for details: \link{exponential_growth}, \link{logistic_growth}, \link{hawkes_growth}, \link{breaks_growth}.
}

\section{Objects from the Class}{
Objects can be created by the functions \code{exponential_growth}, \code{logistic_growth}, \code{breaks_growth}, or \code{hawkes_growth}.
}

\section{Slots}{
  \describe{
    \item{\code{results}:}{Object of class \code{"data.frame"} Model results as a table with coefficents, fit metrics, etc.}
    \item{\code{growth_models}:}{Object of class \code{"list"} containing all models }
    \item{\code{model_type}:}{Object of class \code{"character"} describing the type of model }
    \item{\code{results_cols}:}{Object of class \code{"character"} Vector with column names containing results }
    \item{\code{results_cols_names}:}{Object of class \code{"character"} Vector with descriptions of the column names }
    \item{\code{data_statistics}:}{Object of class \code{"numeric"} Diagnostics of input data }
    \item{\code{time_format}:}{Object of class \code{"character"} Format of time points in time column }
    \item{\code{timestamp}:}{Object of class \code{"list"} Time stamps of any update of the instance }
  }
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "growthmodels")}: Prints an \code{growthmodels} object; use \code{summary(growthmodels)} for results }
    \item{show}{\code{signature(object = "growthmodels")}: Prints an \code{growthmodels} object; use \code{summary(growthmodels)} for results }
    \item{summary}{\code{signature(object = "growthmodels")}: Prints a summary of \code{growthmodels} objects (model results) }
	 }
}

\author{
Thomas Wieland
}

\references{
Bai J, Perron P (2003) Computation and analysis of multiple structural change models. 
\emph{Journal of Applied Econometrics} 18(1), 1-22. 
\doi{https://doi.org/10.1002/jae.659}

Bonifazi G et al. (2021) A simplified estimate of the effective reproduction number Rt using its relation 
with the doubling time and application to Italian COVID-19 data. 
\emph{The European Physical Journal Plus} 136, 386. 
\doi{https://doi.org/10.1140/epjp/s13360-021-01339-6}

Chowell G, Simonsen L, Viboud C, Yang K (2014) Is West Africa Approaching a Catastrophic Phase or is the 2014 
Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. 
\emph{PLoS currents} 6. 
\doi{https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic 
exhibits both global exponential and local polynomial growth rates. 
\emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. 
\doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting 
the 2015 ebola challenge. 
\emph{Epidemics 22}, 62–70. 
\doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L. (2018) SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. In: \emph{Proceedings of the 2018 World Wide Web Conference. 
WWW’18}. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. \doi{https://doi.org/10.1145/3178876.3186108}

Wieland T (2020) A phenomenological approach to assessing the effectiveness 
of COVID-19 related nonpharmaceutical interventions in Germany. 
\emph{Safety Science} 131, 104924. 
\doi{https://doi.org/10.1016/j.ssci.2020.104924}

Wieland T (2020) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
\emph{REGION} 7(2), 43–83. 
\doi{https://doi.org/10.18335/region.v7i2.324}

Zeileis C, Kleiber W, Krämer K, Hornik, K (2003) Testing and dating of structural changes in practice.
\emph{Computational Statistics & Data Analysis} 44(1-2), 109-123. 
\doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
}

\examples{
showClass("growthmodels")
}

\keyword{classes}