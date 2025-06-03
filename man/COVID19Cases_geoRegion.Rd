\name{COVID19Cases_geoRegion}
\docType{data}
\alias{COVID19Cases_geoRegion}
\title{Switzerland Daily COVID-19 cases by region}
\description{
  A dataset containing COVID-19 cases by region (NUTS 3 = cantons) and time periods (days) for Switzerland (Source: Federal Office of Public Health FOPH).
}
\usage{
data(COVID19Cases_geoRegion)
}
\format{
  A data.frame with multiple columns:
  \describe{
    \item{geoRegion}{(character) Region for which the data was collected.}
    \item{datum}{(Date) Date of record.}
    \item{entries}{(integer) Number of reported cases on this date.}
    \item{sumTotal}{(integer) Cumulative case numbers.}
    \item{timeframe_14d}{(logical) Indicates whether the time period covers the last 14 days.}
    \item{timeframe_all}{(logical) Indicates whether the time period covers all previous data.}
    \item{offset_last7d}{(integer) Offset of the last 7 days.}
    \item{sumTotal_last7d}{(integer) Cumulative case numbers of the last 7 days.}
    \item{offset_last14d}{(integer) Offset of the last 14 days.}
    \item{sumTotal_last14d}{(integer) Cumulative case numbers of the last 14 days.}
    \item{offset_last28d}{(integer) Offset of the last 28 days.}
    \item{sumTotal_last28d}{(integer) Cumulative case numbers of the last 28 days.}
    \item{sum7d}{(numeric) Sum of the last 7 days.}
    \item{sum14d}{(numeric) Sum of the last 14 days.}
    \item{mean7d}{(numeric) Average of the last 7 days.}
    \item{mean14d}{(numeric) Average of the last 14 days.}
    \item{entries_diff_last_age}{(integer) Difference from the last age group.}
    \item{pop}{(integer) Population of the region.}
    \item{inz_entries}{(numeric) Incidence of the entries.}
    \item{inzsumTotal}{(numeric) Incidence of cumulative cases.}
    \item{inzmean7d}{(numeric) Incidence of the 7-day average.}
    \item{inzmean14d}{(numeric) Incidence of the 14-day average.}
    \item{inzsumTotal_last7d}{(numeric) Incidence of cumulative cases in the last 7 days.}
    \item{inzsumTotal_last14d}{(numeric) Incidence of cumulative cases in the last 14 days.}
    \item{inzsumTotal_last28d}{(numeric) Incidence of cumulative cases in the last 28 days.}
    \item{inzsum7d}{(numeric) Incidence of the last 7 days.}
    \item{inzsum14d}{(numeric) Incidence of the last 14 days.}
    \item{sumdelta7d}{(numeric) Difference in sums of the last 7 days.}
    \item{inzdelta7d}{(numeric) Difference in incidence of the last 7 days.}
    \item{type}{(character) Type of recorded data (e.g., COVID-19 cases).}
    \item{type_variant}{(character) Variant of the data type.}
    \item{version}{(character) Version of the data collection.}
    \item{datum_unit}{(character) Unit of date specification (e.g., day).}
    \item{entries_letzter_stand}{(integer) Last known count of entries.}
    \item{entries_neu_gemeldet}{(integer) Newly reported entries.}
    \item{entries_diff_last}{(integer) Difference in last entries.}
  }
}
\details{
The data is included as it was published in by the Swiss Federal Office of Public Health (Bundesamt fuer Gesundheit, BAG). Note that the reporting date equals the date of SARS-CoV-2 testing.
}
\source{
  Federal Office of Public Health FOPH (2023) COVID-19 Dashboard Source Data. \url{https://www.covid19.admin.ch/api/data/documentation} (retrieved 2023-06-28)
}
\examples{
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion \%in\% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

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
}