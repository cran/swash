#-------------------------------------------------------------------------------
# Name:        test-swash
# Purpose:     testthat tests for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.0.0
# Last update: 2026-10-10 08:42
# Copyright (c) 2026 Thomas Wieland
#-------------------------------------------------------------------------------


library(swash)
library(testthat)


# Data for tests 1-3:

data(COVID19Cases_geoRegion)
data(K4kant20260101gf_ch2007Poly)

COVID19Cases_geoRegion <-
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion %in% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

test_that("#1 Importing infections panel data as infpan instance and performing analyses", {
  
  # Import as infections panel data set
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
  
  expect_s4_class(infpan_CH, "infpan")
  
  
  # Plot cases
  expect_no_error(
    plot(
      infpan_CH,
      plot_rollmean = TRUE
    )
  )
  
  
  # Calculate effective reproduction number
  infpan_CH <- calculate_Rt(
    infpan_CH,
    verbose = TRUE
  )
  
  expect_true("R_t" %in% names(infpan_CH@input_data))
  
  
  # Calculate rolling mean
  infpan_CH <- calculate_rollmean(
    infpan_CH,
    col_name = "RollingMean",
    verbose = TRUE
  )
  
  expect_true("RollingMean" %in% names(infpan_CH@input_data))
  
  
  # Calculate cumulative values
  infpan_CH <- calculate_cum(
    infpan_CH,
    col_name = "cumulatives",
    verbose = TRUE
  )
  
  expect_true("cumulatives" %in% names(infpan_CH@input_data))
  
  
  # Calculate incidence
  infpan_CH <- calculate_incidence(
    infpan_CH,
    col_name = "incidence",
    verbose = TRUE
  )
  
  expect_true("incidence" %in% names(infpan_CH@input_data))
  
  
  # Add geodata
  expect_no_error(
    infpan_CH <- add_geodata(
      infpan_CH,
      K4kant20260101gf_ch2007Poly,
      unit_col = "cant_abbrev"
    )
  )
  
  # Summary
  expect_no_error(summary(infpan_CH))
  
  
  # Timestamps
  expect_no_error(timestamps(infpan_CH))
  
  
  # Map of incidence
  expect_no_error(
    plot_map(
      infpan_CH,
      attribute = "Incidence",
      main = "COVID-19 Incidence 2020-05-31",
      breaks = c(0, 0.05, 0.1, 0.2, 1),
      verbose = TRUE
    )
  )
  
  
  # Moran's I
  moran_incidence <- expect_no_error(
    spatial_statistic(
      infpan_CH,
      statistic = "moran",
      randomization = FALSE
    )
  )
  
  expect_s4_class(moran_incidence, "nbmatrix")
  
  
  # Logistic growth
  CH_covidwave1_growth <- expect_no_error(
    growth(infpan_CH)
  )
  
  expect_no_error(summary(CH_covidwave1_growth))
  
  
  # Initial exponential growth
  CH_covidwave1_initialgrowth_3weeks <- expect_no_error(
    growth_initial(
      infpan_CH,
      time_units = 21
    )
  )
  
  expect_no_error(summary(CH_covidwave1_initialgrowth_3weeks))
  
  
  # Hawkes process
  CH_covidwave1_Hawkes <- expect_no_error(
    growth_hawkes(infpan_CH)
  )
  
  expect_no_error(summary(CH_covidwave1_Hawkes))
  
  
  # Breakpoints
  CH_covidwave1_breaks <- expect_no_error(
    growth_breaks(infpan_CH)
  )
  
  expect_no_error(summary(CH_covidwave1_breaks))
  
  
  # Swash-Backwash model
  CH_covidwave1 <- expect_no_error(
    swash(
      infpan_CH,
      verbose = TRUE
    )
  )
  
  expect_no_error(summary(CH_covidwave1))
  
})
# Warnings are irrelevant



test_that("#2 Swash-Backwash Model including confidence intervals", {
  
  # Test whether COVID19Cases_geoRegion is balanced panel data
  COVID19Cases_geoRegion_balanced <- expect_no_error(
    is_balanced(
      data = COVID19Cases_geoRegion,
      col_cases = "entries",
      col_date = "datum",
      col_region = "geoRegion"
    )
  )
  
  expect_true(COVID19Cases_geoRegion_balanced$data_balanced)
  
  
  # Swash-Backwash Model for Swiss COVID19 cases
  CH_covidwave1 <- expect_no_error(
    swash_backwash(
      data = COVID19Cases_geoRegion,
      col_cases = "entries",
      col_date = "datum",
      col_region = "geoRegion",
      verbose = TRUE
    )
  )
  
  # Summary
  expect_no_error(
    summary(CH_covidwave1)
  )
  
  
  # Plot
  expect_no_error(
    plot(CH_covidwave1)
  )
  
  
  # Bootstrap confidence intervals
  CH_covidwave1_confint <- expect_no_error(
    confint(
      CH_covidwave1,
      iterations = 100
    )
  )
  
  # Summary of confidence intervals
  expect_no_error(
    summary(CH_covidwave1_confint)
  )
  
  # Plot of confidence intervals
  expect_no_error(
    plot(CH_covidwave1_confint)
  )
  
})


test_that("#3 Testing logistic and exponential growth models", {
  
  # COVID cases for Zurich
  COVID19Cases_ZH <-
    COVID19Cases_geoRegion[
      (COVID19Cases_geoRegion$geoRegion == "ZH") &
        (COVID19Cases_geoRegion$sumTotal > 0),
    ]
  
  
  # Logistic growth model
  loggrowth_ZH <- expect_no_error(
    logistic_growth(
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
  )
  
  # Summary
  expect_no_error(
    summary(loggrowth_ZH)
  )
  
  # Plot
  expect_no_error(
    plot(loggrowth_ZH)
  )
  
  
  # Exponential growth model for the first 4 weeks
  expgrowth_ZH <- expect_no_error(
    exponential_growth(
      y = COVID19Cases_ZH$sumTotal[1:28],
      t = COVID19Cases_ZH$datum[1:28]
    )
  )
  
  # Summary
  expect_no_error(
    summary(expgrowth_ZH)
  )
  
  # Plot
  expect_no_error(
    plot(expgrowth_ZH)
  )
  
  # Doubling rates
  expect_true(
    is.numeric(expgrowth_ZH@GrowthModel_OLS$exp_gr)
  )
  
  expect_true(
    is.numeric(expgrowth_ZH@GrowthModel_NLS$exp_gr)
  )
  
})



# Data for test 4:

data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

test_that("#4 Testing nbmatrix spatial statistics", {
  
  # Creating neighborhood matrix
  Corona_nbmat <- expect_no_error(
    nbmatrix(
      RKI_Corona_counties,
      ID_col = "AGS",
      verbose = TRUE
    )
  )

  
  # Sum of population of neighboring counties
  Corona_nbstat <- expect_no_error(
    nbstat(
      Corona_nbmat,
      link_data = RKI_Corona_counties,
      ID_col = "AGS",
      data_col = "EWZ",
      func = "sum",
      verbose = TRUE
    )
  )

  
  # Plot population
  expect_no_error(
    plot(
      Corona_nbstat,
      main = "No. of inhabitants"
    )
  )

  
  # Define cases_per_ as analysis column
  Corona_nbstat2 <- expect_no_error(
    nbstat(
      Corona_nbmat,
      link_data = RKI_Corona_counties,
      ID_col = "AGS",
      data_col = "cases_per_",
      func = NULL,
      verbose = TRUE
    )
  )

  
  # Global Moran's I
  Corona_SpatialStatistics <- expect_no_error(
    moran(
      Corona_nbstat2,
      verbose = TRUE
    )
  )

  
  # Summary
  expect_no_error(
    summary(Corona_SpatialStatistics)
  )
  
  
  # Global Getis-Ord statistic
  Corona_SpatialStatistics <- expect_no_error(
    getisord(
      Corona_SpatialStatistics,
      verbose = TRUE
    )
  )
  
  
  # Summary
  expect_no_error(
    summary(Corona_SpatialStatistics)
  )

  
  # Local Getis-Ord Gi*
  Corona_SpatialStatistics <- expect_no_error(
    gstar(
      Corona_SpatialStatistics,
      verbose = TRUE
    )
  )

  
  # Summary
  expect_no_error(
    summary(Corona_SpatialStatistics)
  )

  
  # Plot hotspot categories
  expect_no_error(
    plot(
      Corona_SpatialStatistics,
      statistic = "gstar",
      attribute = "Cluster category",
      main = "Hotspot high vs. low",
      pal = c("blue", "red")
    )
  )
  
})
# Warnings are irrelevant


# Data for test 5:
data(Infections)
# Confirmed SARS-CoV-2 cases in Germany

test_that("#5 Testing breakpoints and Hawkes growth models", {
  
  # Breakpoints for time series of infections
  breakpoints_infections <- expect_no_error(
    breaks_growth(
      y = Infections$infections_daily,
      t = Infections$day,
      ln = TRUE,
      verbose = TRUE
    )
  )
  
  # Summary of breakpoints
  expect_no_error(
    summary(breakpoints_infections)
  )
  
  # Plot breakpoints
  expect_no_error(
    plot(breakpoints_infections)
  )
  
  
  # Hawkes Process model
  hawkes_BS <- expect_no_error(
    hawkes_growth(
      y = Infections$infections_daily
    )
  )
  
  # Summary of Hawkes model estimates
  expect_no_error(
    summary(hawkes_BS)
  )
  
  # Plot of Hawkes Process model
  expect_no_error(
    plot(hawkes_BS)
  )
  
})