#---------------------------------------------------------------
# Name:        config (swash package)
# Purpose:     Configuration parameters for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.0.0
# Last update: 2026-09-12 11:27
# Copyright (c) 2025-2026 Thomas Wieland
#---------------------------------------------------------------


# Required packages:
library(sf)
library(spdep)
library(sfdep)
library(zoo)
library(strucchange)
library(lubridate)


# Package name and version:
package_name <- "swash"
package_version <- "3.0.0"
package_title <- "swash: Health Geography Toolbox for Model-Based Analysis of Infections Panel Data"

# Class description texts:
infpan_description <- "Infections Panel Data"
nbmatrix_description <- "Neighborhood matrix with spatial statistics"

# Columns description texts in class infpan:
cases_col_description = "Cases"

# Additional data columns:
permitted_other_cols <- c(
  "R_t", 
  # Effective reproduction number
  "Cum. cases", 
  # Cumulative cases
  "Incidence", 
  # Incidence (per xxx pop)
  "Population",
  # Population size of the region
  "Roll. mean",
  # Rolling mean of cases
  "Roll. sum"
  # Rolling sum of cases
)

permitted_cols <-
  c(
    cases_col_description,
    permitted_other_cols
  )

# Descriptions of model-based analyses in the infpan class
model_descriptions <- 
  list(
    "sbm" = "Swash-Backwash Model for the Single Epidemic Wave",
    "loggrowth" = "Logistic Growth Model",
    "expgrowth" = "Exponential Growth Model",
    "hawkes" = "Hawkes Process",
    "breaksgrowth" = "Time Series Model with Breakpoints"
  )

# Descriptions mapping of spatial statistics in the infpan and nbmatrix classes
spatial_statistics_descriptions <-
  list(
    "nbstat" = "Statistic of all neighboring regions",
    "nbcount" = "Number of neighboring regions",
    "moran" = "Global Moran's I",
    "getisord" = "Global Getis-Ord",
    "gstar" = "Local Getis-Ord Gi*"
  )

# Descriptions mapping of columns in result data frame of function sfdep::local_gstar_perm
gstar_df_colnames <-
  c(
    "gi_star" = "Local Gi*",
    "cluster" = "Cluster category",
    "e_gi" = "Expected",
    "var_gi" = "Variance",
    "std_dev" = "Std. Dev.",
    "p_value" = "p",
    "p_sim" = "p sim.",
    "p_folded_sim" = "p folded sim.",
    "skewness" = "Skewness",
    "kurtosis" = "Kurtosis"
  )