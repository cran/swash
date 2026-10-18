\encoding{UTF-8}
\name{swash-package}
\alias{swash-package}
\docType{package}

\title{
swash: Health Geography Toolbox for Model-Based Analysis of Infections Panel Data 
}

\description{
 
The R library \code{swash} provides a toolbox for quantitative analyses in health geography with respect to the 
spatial spread of infectious diseases. The focus of the package is the analysis of infections panel data,
for example, daily or weekly infections for \eqn{T} time points and \eqn{N} regions within a country.
Examples of research questions include: 
How quickly does an infectious disease spread across a country? 
Are there differences in the spread velocity between countries?
Are there regional differences in the growth rates of the infection wave? 
When was the peak (inflection point) of infections reached in each region?
Are there significant hotspots or clusters of infections?

The target audience of this \code{R} package consists of researchers and practitioners in the fields of health geography, 
spatial epidemiology, and statistics. However, many applications are also relevant to other fields of research, 
such as econometrics or geoinformation science. The package bundles functions and analysis workflows developed by the author 
between 2020 and 2023 for the quantitative analysis of (panel) infection data during the COVID-19 pandemic. 
The aim was to consolidate these methods and analysis tools into a unified and coherent framework.
Users simply need to import their infection panel data once; 
based on this, they can calculate additional indicators and run the built-in analyses successively.

In order to use all functionalities, the user should import her/his infections panel data using the function \code{load_infections_paneldata()}, 
which returns an instance of class \code{infpan}. The panel data is checked whether it is balanced and whether it includes missing values. 
From an \code{infpan} object, the user may utilize the following built-in analysis and visualization functions:

\itemize{
    \item{\emph{Swash-Backwash Model for the Single Epidemic Wave}, including further analysis towards bootstrap-based inference and country comparison as well as visualization}
    \item{\emph{Growth Analysis} with \emph{logistic} and \emph{exponential growth models} for cumulative or new (e.g., daily, weekly) infections, 
    including visualization}
    \item{\emph{Hawkes process} models for new infection cases (e.g., daily, weekly), including visualization}
    \item{\emph{Breakpoints analysis} using the \emph{Bai-Perron algorithm} implemented in \code{strucchange::breakpoints}, 
    including visualization}
    \item{Spatial statistics for the identification of clusters/hotspots, such as \emph{Moran's I} and \emph{Getis-Ord}, 
    implemented in \code{spdep} and \code{sfdep}}
    \item{Calculation of further epidemic indicators from the infections panel data such as the \emph{effective reproduction number}}
    \item{\emph{Plots} of infection curves by region}
    \item{Map plots of infection indicators and (local) spatial statistics}
}

\code{infpan} objects and objects resulting from the functions mentioned above have \code{summary()} and \code{plot()} methods. 
All mentioned functions may be used stand-alone as well. 
Every change (e.g., the calculation of a new indicator) within an object is tracked with a timestamp.
}

\details{
Based on an \code{infpan} object, several indicators may be calculated from new infections values, such as incidence or effective reproduction number \eqn{R_t}.
Infection curves may be plotted by \code{plot(infpan)}. It is also possible to import spatial features (\code{sf} spatial data frame) into an \code{infpan} object
in order to plot maps of specified indicator at a given time point. All built-in model analyses mentioned below may be conducted based on an instance of class \code{infpan},
but may also be used independently. 

The Swash-Backwash Model for the Single Epidemic Wave (SBM) is the spatial equivalent of the classic epidemiological SIR (Susceptible-Infected-Recovered) model. 
It was developed by Cliff and Haggett (2006) to model the velocity of spread of infectious diseases across space. 
Current applications may be found, for example, in Smallman-Raynor et al. (2022a,b). 
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
Logistic and exponential growth models are estimated via OLS (Ordinary Least Squares) and/or NLS (Nonlinear least squares).
All model estimations include the calculation of goodness-of-fit metrics (e.g., \eqn{R^2}, \eqn{RMSE}, \eqn{MAE}, \eqn{MAPE}).
Additionally, breakpoints in time series (see, e.g., Wieland 2020b) may be detected based on the \emph{Bai-Perron algorithm} (Bai and Perron 2003) 
which is included in the \code{breakpoints} function of the \code{strucchange} package (Zeileis et al. 2003) 
accessed in the functions of this package.
A model for a single time series may be estimated with the built-in functions \code{logistic_growth()}, \code{exponential_growth()}, \code{hawkes_growth()}, 
or \code{breaks_growth()}, respectively. These function return objects of class \code{loggrowth}, \code{expgrowth}, \code{hawkes},
and \code{breaksgrowth}, respectively, all of them defined in this package. Plotting is available via the \code{plot()} method.
Estimating such a model based on an \code{infpan} object is provided by the \code{infpan} methods \code{growth()}, \code{growth_initial()}, 
\code{growth_hawkes()}, and \code{growth_breaks()}, respectively, all of them resulting in an object of class \code{growthmodels}.

Spatial statistics regarding autocorrelation or clustering can be calculated either based on an \code{infpan} object or independently of it,
with the class \code{nbmatrix} existing for the latter.
Neighboring spatial features can be identified, and descriptive statistics (e.g., sum, mean) can be calculated for the neighboring regions. 
Within the \code{nbmatrix} class, statistics such as Moran's I and Getis-Ord may be calculated.
The methods \code{moran()}, \code{getisord()}, and \code{gstar()} are available for this purpose.
The \code{swash} package internally uses the packages \code{spdep} (Bivand and Wong 2018) and \code{sfdep} (Parry and Locke 2024).

The library also includes various auxiliary statistical functions that are integrated into the main functions and methods 
but can also be used standalone. 
For example, \code{metrics()} and \code{binary_metrics()} calculate goodness-of-fit metrics for observed and predicted values 
for continuous and binary outcome variables, respectively.
The \code{plot_coef_ci()} function allows regression model results (estimates and confidence intervals) to be visualized in an appealing way.
The package includes example data from the COVID-19 pandemic.
}

\author{
Thomas Wieland
}

\references{
Bivand RS, Wong DWS (2018) 
Comparing implementations of global and local indicators of spatial association. 
\emph{TEST} 27, 716–748. \doi{https://doi.org/10.1007/s11749-018-0599-x}

Chowell G, Simonsen L, Viboud C, Yang K (2014) 
Is West Africa Approaching a Catastrophic Phase or is the 2014 Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. 
\emph{PLoS currents 6}. \doi{https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81}

Chowell G, Viboud C, Hyman JM, Simonsen L (2015) 
The Western Africa ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. 
\emph{PLOS Currents Outbreaks}, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. \doi{https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261}

Cliff AD, Haggett P (2006) 
A swash-backwash model of the single epidemic wave. 
\emph{Journal of Geographical Systems} 8(3), 227-252. \doi{https://doi.org/10.1007/s10109-006-0027-8}

Li, MY (2018) 
An Introduction to Mathematical Modeling of Infectious Diseases. 
\doi{https://doi.org/10.1007/978-3-319-72122-4}

Nishiura H, Chowell G (2009) 
The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. 
In Chowell G, Hyman JM, Bettencourt LMA (eds.) Mathematical and statistical estimation approaches in epidemiology, 103–121. 
\doi{https://doi.org/10.1007/978-90-481-2313-1_5}

Parry J, Locke D (2024) 
sfdep: Spatial Dependence for Simple Features. R package version 0.2.5.
\doi{https://doi.org/10.32614/CRAN.package.sfdep}
  
Pell B, Kuang Y, Viboud C, Chowell G (2018) 
Using phenomenological models for forecasting the 2015 ebola challenge. 
\emph{Epidemics} 22, 62–70. \doi{https://doi.org/10.1016/j.epidem.2016.11.002}

Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L. (2018) SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. In: \emph{Proceedings of the 2018 World Wide Web Conference. 
WWW’18}. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. \doi{https://doi.org/10.1145/3178876.3186108}

Rotejanaprasert C, Chinpong K, Lawson AB, Chienwichai P, Maude RJ (2024)
Evaluation and comparison of spatial cluster detection methods for improved decision making of disease surveillance: 
a case study of national dengue surveillance in Thailand. 
\emph{BMC Medical Research Methodology} 24, 14. \doi{https://doi.org/10.1186/s12874-023-02135-9}

Smallman-Raynor MR, Cliff AD, Stickler PJ (2022a)
Meningococcal Meningitis and Coal Mining in Provincial England: Geographical Perspectives on a Major Epidemic, 1929–33. 
\emph{Geographical Analysis} 54, 197–216. \doi{https://doi.org/10.1111/gean.12272}

Smallman-Raynor MR, Cliff AD, The COVID-19 Genomics UK (COG-UK) Consortium (2022b) 
Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021.
\emph{Epidemiology and Infection} 150, e145. \doi{https://doi.org/10.1017/S0950268822001285}.

Viboud C, Bjørnstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT (2006) 
Synchrony, Waves, and Spatial Hierarchies in the Spread of Influenza. 
\emph{Science} 312, 447-451. \doi{https://doi.org/10.1126/science.1125237}

Wieland T (2020a) 
Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
\emph{REGION} 7(2), 43–83. \doi{https://doi.org/10.18335/region.v7i2.324}

Wieland T (2020b) 
A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. 
\emph{Safety Science} 131, 104924. \doi{https://doi.org/10.1016/j.ssci.2020.104924}

Wieland T (2022) 
Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. 
\emph{European Journal of Geography} 13(4), 18-33. \doi{https://doi.org/10.48088/ejg.t.wie.13.4.018.033}

Wieland T (2025) 
Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. 
\emph{Journal of Public Health} 33(11), 2497-2511. \doi{https://doi.org/10.1007/s10389-024-02218-x}

Zeileis C, Kleiber W, Krämer K, Hornik K (2003) Testing and dating of structural changes in practice. 
\emph{Computational Statistics & Data Analysis} 44(1-2), 109-123. 
\doi{https://doi.org/10.1016/S0167-9473(03)00030-6}
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