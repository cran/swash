#-------------------------------------------------------------------------------
# Name:        swash_test
# Purpose:     Tests and examples for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     2.2.0
# Last update: 2026-10-17 10:12
# Copyright (c) 2025-2026 Thomas Wieland
#-------------------------------------------------------------------------------

library(swash)


# Switzerland:

data(COVID19Cases_geoRegion)
data(K4kant20260101gf_ch2007Poly)

table(COVID19Cases_geoRegion$geoRegion)
table(COVID19Cases_geoRegion$datum)

COVID19Cases_geoRegion <-
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion %in% c("CH", "CHFL"),]
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

infpan_CH <- calculate_rollmean(
  infpan_CH, 
  col_name = "RollingMean",
  verbose = TRUE
)
# Calculate rolling mean of cases as "RollingMean"

infpan_CH <- calculate_cum(
  infpan_CH, 
  col_name = "cumulatives",
  verbose = TRUE
)
# Calculate cumulative values of cases as "cumulatives"

infpan_CH <- calculate_incidence(
  infpan_CH, 
  col_name = "incidence",
  verbose = TRUE
)
# Calculate incidence of cases as "incidence"

infpan_CH <- add_geodata(
  infpan_CH,
  K4kant20260101gf_ch2007Poly,
  unit_col = "cant_abbrev"
)
# Adding geodata (sf) to infpan object

summary(infpan_CH)
# Summary of infpan object

timestamps(infpan_CH)
# Time stamps of infpan object

plot_map(
  infpan_CH,
  attribute = "Incidence",
  main = "COVID-19 Incidence 2020-05-31",
  breaks = c(0, 0.05, 0.1, 0.2, 1),
  verbose = TRUE
)
# Simple map of incidence for the most current date (timepoint=NULL)

moran_incidence <-
  spatial_statistic(
    infpan_CH,
    statistic = "moran",
    randomization = FALSE
  )
# Calculating Moran's I test under normality
# Result = nbmatrix object


CH_covidwave1_growth <- 
  growth(infpan_CH)
CH_covidwave1_growth
summary(CH_covidwave1_growth)
# Logistic growth models for infpan object infpan_CH

CH_covidwave1_initialgrowth_3weeks <- 
  growth_initial(
    infpan_CH,
    time_units = 21
  )
CH_covidwave1_initialgrowth_3weeks
summary(CH_covidwave1_initialgrowth_3weeks)
# Exponential models for infpan object CH_covidwave1 
# initial growth in the first 3 weeks


CH_covidwave1_Hawkes <- 
  growth_hawkes(infpan_CH)
CH_covidwave1_Hawkes
summary(CH_covidwave1_Hawkes)
# Hawkes process models for infpan object infpan_CH 


CH_covidwave1_breaks <- 
  growth_breaks(infpan_CH)
CH_covidwave1_breaks
summary(CH_covidwave1_breaks)
# Breakpoints for infpan object infpan_CH 


CH_covidwave1 <-
  swash(
    infpan_CH,
    verbose = TRUE
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model

# infpan_CH@timestamp

COVID19Cases_geoRegion_balanced <- 
  is_balanced(
  data = COVID19Cases_geoRegion,
  col_cases = "entries",
  col_date = "datum",
  col_region = "geoRegion"
)
# Test whether "COVID19Cases_geoRegion" is balanced panel data 

COVID19Cases_geoRegion_balanced$data_balanced
# Balanced? TRUE or FALSE

CH_covidwave1 <- 
  swash_backwash(
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion",
    verbose = TRUE
  )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model

CH_covidwave1 <- 
  swash_backwash(
    infpan=infpan_CH,
    verbose = TRUE
  )
# Same Swash-Backwash Model analysis
# with infpan object

plot(CH_covidwave1)
# Plot of Swash-Backwash Model edges and total epidemic curve

plot(
  infpan_CH,
  normalize_by_col = "pop",
  plot_rollmean = TRUE
  )

CH_covidwave1_confint <- 
  confint(
    CH_covidwave1, 
    iterations = 100
    )
# Bootstrap confidence intervals with 100 iterations

summary(CH_covidwave1_confint)
# Summary of confidence intervals

plot(CH_covidwave1_confint)
# Plot of confidence intervals


# Austria:

data(Oesterreich_Faelle)

table(Oesterreich_Faelle$NUTS3)
table(Oesterreich_Faelle$Datum)

AT_covidwave1 <- 
  swash_backwash(
    data = Oesterreich_Faelle,
    col_cases = "Faelle",
    col_date = "Datum",
    col_region = "NUTS3"
  )
# Swash-Backwash Model for Austrian COVID19 cases
# Spatial aggregate: NUTS 3

summary(AT_covidwave1)

plot(AT_covidwave1)


AT_vs_CH <- 
  compare_countries(
    CH_covidwave1, 
    AT_covidwave1,
    country_names = c("Switzerland", "Austria"),
    iterations = 10
    )

AT_vs_CH

plot(AT_vs_CH)


COVID19Cases_ZH <-
  COVID19Cases_geoRegion[
    (COVID19Cases_geoRegion$geoRegion == "ZH")
    & (COVID19Cases_geoRegion$sumTotal > 0),]
# COVID cases for Zurich


loggrowth_ZH <- logistic_growth(
  y = COVID19Cases_ZH$sumTotal, 
  t = COVID19Cases_ZH$datum, 
  S = 3600,
  S_start = NULL, 
  S_end = NULL, 
  S_iterations = 10, 
  S_start_est_method = "bisect", 
  seq_by = 10,
  nls = TRUE
)
# Logistic growth model with stated saturation value

summary(loggrowth_ZH)
# Summary of logistic growth model estimates

plot(loggrowth_ZH)
# Plot of logistic growth model


Rt_BS <- R_t(infections = COVID19Cases_ZH$entries)
Rt_BS
# Effective reproduction number


expgrowth_ZH <- exponential_growth (
  y = COVID19Cases_ZH$sumTotal[1:28], 
  t = COVID19Cases_ZH$datum[1:28] 
)
# Exponential growth model for the first 4 weeks

summary(expgrowth_ZH)
# Summary of exponential growth model

plot(expgrowth_ZH)
# Plot of exponential growth model

expgrowth_ZH@GrowthModel_OLS$exp_gr
# Doubling rate (OLS fit)
expgrowth_ZH@GrowthModel_NLS$exp_gr
# Doubling rate (NLS fit)


data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

Corona_nbmat <- 
  nbmatrix (
    RKI_Corona_counties, 
    ID_col="AGS",
    verbose = TRUE
  )
# Creating neighborhood matrix

Corona_nbstat <- nbstat(
  Corona_nbmat,
  link_data = RKI_Corona_counties, 
  ID_col = "AGS", 
  data_col = "EWZ", 
  func = "sum",
  verbose = TRUE
  )
# Sum of population (EWZ) of neighboring counties

plot(
  Corona_nbstat,
  main = "No. of inhabitants"
)
# Plot simple map of "EWZ"

Corona_nbstat2 <- nbstat(
  Corona_nbmat,
  link_data = RKI_Corona_counties, 
  ID_col = "AGS", 
  data_col = "cases_per_", 
  func = NULL,
  verbose = TRUE
)
# New nbmatrix instance
# Defining "cases_per_" (COVID-19 cases per 100,000 inhabitants) as analysis column

Corona_SpatialStatistics <- moran(
  Corona_nbstat2,
  verbose = TRUE
  )
# Calculating Global Moran's I

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object

Corona_SpatialStatistics <- getisord(
  Corona_SpatialStatistics,
  verbose = TRUE
)
# Calculating Global Getis-Ord statistic

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object

Corona_SpatialStatistics <- gstar(
  Corona_SpatialStatistics,
  verbose = TRUE
)
# Calculating Local Getis-Ord Gi*

summary(Corona_SpatialStatistics)
# Summary of nbmatrix object

plot(
  Corona_SpatialStatistics,
  statistic = "gstar",
  attribute = "Cluster category",
  main = "Hotspot high vs. low",
  pal = c("blue", "red")
  )


data(did_fatalities_splm_coef)
# Results of a difference-in-differences model

plot_coef_ci(
  point_estimates = did_fatalities_splm_coef$Estimate,
  confint_lower = did_fatalities_splm_coef$CI_lower_Bonferroni,
  confint_upper = did_fatalities_splm_coef$CI_upper_Bonferroni,
  coef_names = did_fatalities_splm_coef$Var,
  skipvars = c(
    "Alpha_share", 
    "lambda",
    "rho",
    "log(D_Infections_daily_7dsum_per100000_lag2weeks)",
    "vacc_cum_per100000_lag2weeks"
    ),
  lwd = 13,
  pch = 19,
  auto_color = TRUE
)
# Plot with point estimates and confidence intervals

plot_coef_ci(
  point_estimates = did_fatalities_splm_coef$Estimate,
  confint_lower = did_fatalities_splm_coef$CI_lower_Bonferroni,
  confint_upper = did_fatalities_splm_coef$CI_upper_Bonferroni,
  coef_names = did_fatalities_splm_coef$Var,
  p = did_fatalities_splm_coef$Pr_t_Bonferroni,
  skipvars = c(
    "Alpha_share", 
    "lambda",
    "rho",
    "log(D_Infections_daily_7dsum_per100000_lag2weeks)",
    "vacc_cum_per100000_lag2weeks"
  ),
  lwd = 13,
  pch = 19,
)
# Plot with point estimates and confidence intervals


data(Infections)
# Confirmed SARS-CoV-2 cases in Germany

breakpoints_infections <- breaks_growth(
  y = Infections$infections_daily,
  t = Infections$day,
  ln = TRUE,
  verbose = TRUE
)
# Breakpoints for time series of infections

summary(breakpoints_infections)
# Summary of breakpoints

plot(breakpoints_infections)
# Plot breakpoints


hawkes_BS <- hawkes_growth(
  y = Infections$infections_daily
)
# Hawkes Process model

summary(hawkes_BS)
# Summary of Hawkes model estimates

plot(hawkes_BS)
# Plot of Hawkes Process model
