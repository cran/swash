\encoding{UTF-8}
\name{swash-package}
\alias{swash-package}
\docType{package}

\title{
swash: Health Geography Toolbox for Model-Based Analysis of Infections Panel Data 
}

\description{

The R library \emph{swash} provides a toolbox for quantitative analyses in health geography with a focus on the spatial spread of infectious diseases.
It bundles functions developed by the author between 2020 and 2023 for the quantitative analysis of (panel) infection data during the COVID-19 pandemic.
The aim was to consolidate these methods and analytical tools into a unified and coherent framework.
The target audience of this R package consists of researchers and practitioners in the fields of health geography, spatial epidemiology, and statistics.

In order to use all functionalities, the user should import her/his infections panel data using the function \code{load_infections_paneldata()}, 
which returns an instance of class \code{infpan}. The panel data is checked whether it is balanced and whether it includes missing values. 
From an \code{infpan} object, the user may utilize the following built-in analysis models and visualization functions:

\itemize{
    \item{\emph{Swash-Backwash Model for the Single Epidemic Wave}, including further analysis towards bootstrap-based inference and country comparison as well as visualization}
    \item{\emph{Growth Analysis} with \emph{logistic} and \emph{exponential growth models} for cumulative or incremental infections, 
    whereby the former is intended for the entire infection wave, and the latter for the initial phase of the infection wave; including visualization}
    \item{\emph{Hawkes process} models for incremental infections; including visualization}
    \item{\emph{Breakpoints analysis} using the \emph{Bai-Perron algorithm} implemented in \code{strucchange::breakpoints}; including visualization}
    \item{Calculation of further epidemic indicators from the infections panel data such as the \emph{effective reproduction number}}
    \item{\emph{Plots} of infection curves by region}
}

\code{infpan} objects and objects resulting from the functions mentioned above have summary() and plot() methods. 
All mentioned functions may be used stand-alone as well.
}

\details{
Based on an \code{infpan} object, several indicators may be calculated from incremental infections values, such as incidence or effective reprodction number \eqn{R_t}.
Infection curves may be plotted by \code{plot(infpan)}. All built-in model analyses may be conducted based on an instance of class \code{infpan}.

The Swash-Backwash Model (SBM) for the Single Epidemic Wave is the spatial equivalent of the classic epidemiological SIR (Susceptible-Infected-Recovered) model. 
It was developed by Cliff and Haggett (2006) to model the velocity of spread of infectious diseases across space. 
Current applications can be found, for example, in Smallman-Raynor et al. (2022a,b). 
This package enables the calculation of the Swash-Backwash Model for user-supplied panel data on regional infections. 
The core of this is the \code{swash_backwash()} function, which calculates the model and creates a model object of the \code{sbm} class defined in this package. 
This class can be used to visualize results (\code{summary()}, \code{plot()}) and calculate bootstrap confidence intervals for the model estimates (\code{confint(sbm)}); 
the latter returns an object of the built-in class \code{sbm_ci}. Two \code{sbm_ci} objects for different countries may be compared with \code{compare_countries()}, 
which allows the estimation of mean differences of a user-specified model parameter (e.g., spatial reproduction number \eqn{R_{OA}}) between two countries. 
This makes it possible to check whether the spatial spread velocity of a communicable disease is significantly different in one country than in another country; 
the result is an object of class \code{countries}. To calculate the SBM model based on an \code{infpan} object, use the corresponding method \code{swash(infpan)}.

The library allows for estimating growth models based on time series of infections. 
Logistic and exponential growth models (see, e.g., Chowell et al. 2014, 2015, Pell et al. 2018, Wieland 2020a, 2020b) 
as well as Hawkes process models (see, e.g. Rizoiu et al. 2018) are provided. 
Additionally, breakpoints in time series may be detected (see, e.g., Wieland 2020b).
A model for a single time series may be estimated with the built-in functions \code{logistic_growth()}, \code{exponential_growth()}, \code{hawkes_growth()}, 
or \code{breaks_growth()}, respectively. These function return objects of class \code{loggrowth}, \code{expgrowth}, \code{hawkes},
and \code{breaksgrowth}, respectively, all of them defined in this package. Plotting is available via the \code{plot} method.
Estimating such a model based on an \code{infpan} object is provided by the \code{infpan} methods \code{growth()}, \code{growth_initial()}, 
\code{growth_hawkes()}, and \code{growth_breaks()}, respectively, all of them resulting in an object of class \code{growthmodels}.

The package also contains other functions for spatio-temporal analysis, including spatial statistics (\code{nbstat()} for neighborhood statistics)
and fit metrics (\code{metrics()}, \code{binary_metrics()}, \code{binary_metrics_glm()}).
The package includes example data from the COVID-19 pandemic.
}

\author{
Thomas Wieland
}

\references{
Chowell G, Simonsen L, Viboud C, Yang K (2014) Is West Africa Approaching a Catastrophic Phase or is the 2014 Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. \emph{PLoS currents 6}. \doi{https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. \emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Cliff AD, Haggett P (2006) A swash-backwash model of the single epidemic wave. \emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Li, MY (2018) An Introduction to Mathematical Modeling of Infectious Diseases. \doi{https://doi.org/10.1007/978-3-319-72122-4}

Nishiura H, Chowell G (2009) The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. In Chowell G, Hyman JM, Bettencourt LMA (eds.) Mathematical and statistical estimation approaches in epidemiology, 103–121. \doi{https://doi.org/10.1007/978-90-481-2313-1_5}

Pell B, Kuang Y, Viboud C, Chowell G (2018) Using phenomenological models for forecasting the 2015 ebola challenge. \emph{Epidemics 22}, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L. (2018) SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. In: \emph{Proceedings of the 2018 World Wide Web Conference. 
WWW’18}. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. \doi{https://doi.org/10.1145/3178876.3186108}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022a) Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. \emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022b) Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. \emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

Viboud C, Bjørnstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT (2006) Synchrony, Waves, and Spatial Hierarchies in the Spread of Influenza. \emph{Science} 312, 447-451. \doi{https://doi.org/10.1126/science.1125237}

Wieland T (2020a) Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. \emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2020b) A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. \emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Wieland T (2022) Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. \emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}

Wieland T (2025) Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. \emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}
}

\keyword{package}

\examples{
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_geoRegion <-
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion \%in\% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

infpan_CH <- load_infections_paneldata(
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion",
    other_cols = c(
      "Population" = "pop"
        ), 
    verbose = TRUE
  )
# Import as infections panel data set (class infpan)

is(infpan_CH)
# "infpan"

plot(
  infpan_CH,
  plot_rollmean = TRUE
  )
# Plot cases

infpan_CH <- calculate_Rt(
  infpan_CH,
  verbose = TRUE
  )
# Calculate effective reproduction number

summary(infpan_CH)
# Summary of infpan object

timestamps(infpan_CH)
# Time stamps of infpan object

CH_covidwave1 <-
  swash(
    infpan_CH,
    verbose = TRUE
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model
}