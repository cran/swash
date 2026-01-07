\name{C19dNUTSdata}
\alias{C19dNUTSdata}
\docType{data}
\title{
Regional cumulative COVID-19 deaths
}

\description{
Cumulative COVID-19 deaths absolute and per 100,000 pop at NUTS3 level for 31 EU/EFTA countries
}


\usage{data("C19dNUTSdata")}

\format{
  A data frame with 1,143 observations (each one represents a spatial NUTS unit).
  \describe{
    \item{\code{NUTS_ID}}{NUTS ID of the spatial unit}
    \item{\code{CNTR_CODE}}{Country code (= NUTS 0 ID) of the given spatial unit}
    \item{\code{NUTS_Level}}{NUTS level of the given spatial unit (0 = national, 1, 2, 3)}
    \item{\code{NUTS2_ID}}{NUTS 2 ID of the spatial unit}
    \item{\code{NUTS1_ID}}{NUTS 1 ID of the spatial unit}
    \item{\code{NUTS_Name}}{Latin name of the spatial unit}
    \item{\code{C19deaths}}{Cumulative COVID-19 deaths [persons]}
    \item{\code{pop2020}}{Population in 2020 [persons]}
    \item{\code{C19deaths_per100000}}{Cumulative COVID-19 deaths [per 100,000]}
    \item{\code{annotation}}{Annotation}
)}
}

\details{

\emph{Note: This data was originally released in the author's package C19dNUTS in 2022 (https://cran.r-project.org/package=C19dNUTS). Some of the URLs referred to here were moved or deleted.}

The dataset contains cumulative COVID-19 deaths at the regional level (mostly NUTS 3, N=1,143) for 31 EU/EFTA countries (AT, BE, BG, CH, CY, CZ, DE, DK, EE, EL, ES, FI, FR, HR, HU, IE, IS, IT, LT, LU, LV, MT, NL, NO, PL, PT, RO, SE, SI, SK, UK). The \code{C19deaths} variable contains the absolute number of COVID-19 related deaths, and the variable \code{C19deaths_per100000} equals the death numbers relative to the population (per 100,000). 

Unless otherwise noted, data includes all reported COVID-19 related deaths since the beginning of the COVID-19 pandemic through June 2022. Please refer to the source section below for the exact date on which each raw dataset was retrieved. The spatial level is the current NUTS 2021 classification of the European Union (see `https://ec.europa.eu/eurostat/web/nuts/background`), with one slight modification (see "Technical details" below). The variable \code{NUTS_Level} documents the spatial level for which the numbers apply (mostly \code{NUTS_Level} = 3 for NUTS3). 

\emph{Technical details:}

This dataset contains cumulative numbers and no time series, as many countries only publish cumulative data on COVID-19 deaths. In cases where countries only publish COVID-19 deaths in the form of daily data, the numbers were summed up over the entire period under consideration at the respective spatial level.

The definition of a COVID-19 death may vary between countries. The respective definition can usually be found on the website of the national health authority. In some countries, data is reported based on different definitions. For example, Lithuania uses three different definitions, namely a) based on the main cause of death in the death certificate, b) based on a mention in the death certificate and c) died within 28 days of a positive SARS-CoV-2 test (\url{https://open-data-sets-ls-osp-sdg.hub.arcgis.com/datasets/ba35de03e111430f88a86f7d1f351de6_0/about}). In England, for example, a distinction is made between the deceased who tested positive and those who died from COVID-19 based on the death certificate (https://coronavirus.data.gov.uk/details/deaths). In these cases, the definition used has always been the equivalent of the total number of COVID-19 deaths as reported by the national figures from Johns Hopkins University (\url{https://coronavirus.jhu.edu/data/cumulative-cases}).

In some cases, countries publish regional COVID-19 data directly at NUTS3 level (e.g., Germany) or NUTS2 level (e.g., Italy). In most cases, the regional level had to be linked manually using the name of the region (e.g., Bulgaria, Norway, Switzerland). Some countries even publish the relevant data on a smaller scale, i.e. below NUTS3 (e.g., Austria, Netherlands, Poland, England). In these cases, where a reference table (subnational spatial unit <-> NUTS3) was available, the lower level was linked to the NUTS3 level (e.g., England). If no reference table but geodata (shapefiles) for the lower spatial level was available (e.g., Austria, Netherlands, Poland), the lower level was linked to the NUTS3 level via a spatial join (Polygon centroids; in cases where the centroid was outside the polygon, it was placed inside the polygon manually). In these cases, the numbers were then summed up at NUTS3 level. 

The spatial reference used here is the current EU NUTS Shapefile (https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts; accessed 2022-06-23). The dataset can be linked directly to this shapefile, where the unique id field to which the link can be made is the column \code{NUTS_ID}. However, there is one exception: To ensure data compatibility, the UK NUTS3 regions UKM61 and UKM63 were aggregated into one region (UKM61).

The data reflects 1,309,326 COVID-19 related deaths in the 31 countries in the investigated time period. The variable \code{C19deaths_per100000} is non-normally distributed (Shapiro-Wilk test: W = 0.92284, p < 0.01). The natural log of \code{C19deaths_per100000} is spatially autocorrelated (Moran's I with queen contiguous spatial weighting: I = 0.65228, p < 0.01).

\emph{Data limitations:}

It can be assumed that there are differences between countries and possibly also over time in the definition of a COVID-19 death (see "Technical details" above). Please check the definition on the website of the respective national health authority.

Data on COVID-19 deaths are incomplete for the following EU/EFTA countries: Bulgaria, France, Poland. In Bulgaria, regional COVID-19 deaths were only published for the years 2020 and 2021 (36,142 COVID-19 related deaths in total), i.e. the cases for 2022 are missing. France only publishes the COVID-19 patients which died in a hospital at the regional level, which equals 120,630 COVID-19 related deaths over the period under consideration (as of 2022-06-30). The total number of COVID-19 related deaths in France for the same time is equal to 149,533, which means that there is a lack of 28,903 COVID-19 fatalities (19.3 \%, e.g., people which died in nursing homes). Polish deaths are missing COVID-19 deaths from the first pandemic wave. Therefore only the COVID-19 fatalities from the date 2020-11-24 are included, which equals 102,449 deaths. In the previous period, 13,780 COVID-19 deaths were reported, which are not included in the data set, i.e. 11.9\% of the deaths are missing.

Of the 31 EU/EFTA countries included, regional data are only available for 24 countries. The following countries have not published sub-national data for COVID-19 deaths: Cyprus, Finland, Island, Hungary, Estonia, Latvia, Malta. The values for Finland, Hungary, Estonia and Latvia refer to the national level (NUTS 0), which is indicated by the variable \code{NUTS_Level} = 0. In the cases of Cyprus, Malta and Iceland (which are rather small countries), the NUTS 0 level also corresponds to the NUTS 2 level, which is why they are marked here in the dataset with \code{NUTS_Level} = 2. It is comparatively difficult to compare the data with Belgium because COVID-19 death figures are only published there at NUTS 1 level (3 regions; \code{NUTS_Level} = 1).

Some countries report separately persons who died of/with COVID-19 who live outside the country or cannot be assigned to a region (e.g., Greece, Norway). These cases are shown separately in the dataset, but cannot be related to population numbers and cannot be linked to the NUTS shapefile.

Norway does not provide COVID-19 data for the NUTS3 regions NO0B1 and NO0B2. 

In the UK, each country (England, Wales, Scotland, and Northern Ireland) is independently responsible for publishing COVID-19 data. Therefore the data are not all available at the same spatial aggregation level (e.g. England: NUTS 3, Wales: NUTS 2). 
}

\source{
Raw data of COVID-19 deaths:

\emph{Note: Some of the URLs have been moved or deleted.}

AT: https://covid19-dashboard.ages.at/data/CovidFaelle_Timeline_GKZ.csv (accessed 2022-06-23)

BE: \url{https://epistat.sciensano.be/Data/COVID19BE_MORT.csv} (accessed 2022-06-21)

BG: \url{https://www.nsi.bg/sites/default/files/files/data/table/COVID_2020_2021_EN.xls} (accessed 2022-06-29)

CH: \url{https://www.covid19.admin.ch/api/data/20220621-t6j901v4/downloads/sources-csv.zip} (accessed 2022-06-21)

CY: \url{https://covid19.who.int/region/euro/country/cy} (accessed 2022-06-30)

CZ: \url{https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19/umrti.csv} (accessed 2022-06-24)

DE: \url{https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/917fc37a709542548cc3be077a786c17_0/about} (accessed 2022-06-23)

DK: \url{https://files.ssi.dk/covid19/overvagning/dashboard/overvaagningsdata-dashboard-covid19-28062022-byis} (accessed 2022-06-29), folder: "Regionalt_DB", file: "07_antal_doede_pr_dag_pr_region"

EE: https://www.terviseamet.ee/en/coronavirus/coronavirus-dataset (accessed 2022-07-11)

EL: \url{https://github.com/Sandbird/covid19-Greece} (accessed 2022-07-02), file "regions"

ES: \url{https://cnecovid.isciii.es/covid19/resources/casos_hosp_uci_def_sexo_edad_provres.csv} (accesssed 2022-06-28)

FI: \url{https://covid19.who.int/region/euro/country/fi} (accessed 2022-07-01)

FR: \url{https://www.data.gouv.fr/fr/datasets/synthese-des-indicateurs-de-suivi-de-lepidemie-covid-19/} (accessed 2022-07-01), file "table-indicateurs-open-data-dep-2022-06-30-19h00"

HR: https://www.koronavirus.hr/zupanije/139 (accessed 2022-06-28)

HU: \url{https://covid19.who.int/region/euro/country/hu} (accessed 2022-07-02)

IE: \url{https://epi-covid-19-hpscireland.hub.arcgis.com/} (accessed 2022-06-29)

IS: https://www.covid.is/data (accessed 2022-06-27)

IT: \url{https://github.com/pcm-dpc/COVID-19/tree/master/dati-regioni} (accessed 2022-06-24), file "dpc-covid19-ita-regioni-latest_raw"

LV: \url{https://covid19.gov.lv/en/node/16387} (accessed 2022-07-27)

LT: \url{https://open-data-ls-osp-sdg.hub.arcgis.com/datasets/ba35de03e111430f88a86f7d1f351de6_0/explore} (accessed 2022-06-27)

LU: \url{https://covid19.public.lu/fr/graph.html} (accessed 2022-06-27)

MT: \url{https://covid19.who.int/table} (accessed 2022-07-01)

NL: \url{https://data.rivm.nl/covid-19/COVID-19_aantallen_gemeente_per_dag.csv} (accessed 2022-06-27)

NO: https://www.fhi.no/contentassets/8a971e7b0a3c4a06bdbf381ab52e6157/vedlegg/2022/ukerapport-uke-20-16.05---22.05.22.pdf (accessed 2022-07-07)

PL: \url{https://www.gov.pl/web/koronawirus/wykaz-zarazen-koronawirusem-sars-cov-2} (accessed 2022-06-23)

PT: https://github.com/dssg-pt/covid19pt-data/blob/master/data.csv (accessed 2022-06-29)

RO: https://covid19.geo-spatial.org/?map=decese (accessed 2022-07-01)

SE: \url{https://experience.arcgis.com/experience/19fc7e3f61ec4e86af178fe2275029c5} (accessed 2022-06-23)

SI: \url{https://www.nijz.si/sites/www.nijz.si/files/uploaded/tedenski_prikaz_umrli20220627.xlsx} (accessed 2022-06-28)

SK: \url{https://github.com/Institut-Zdravotnych-Analyz/covid19-data} (accessed 2022-06-28), folder "Deaths", file "OpenData_Slovakia_Covid_Deaths_AgeGroup_District"

UK - England: https://coronavirus.data.gov.uk/details/deaths (accessed 2022-06-24), file "ltla_2022_06_23_cumDeaths60DaysByDeathDate_ref"

UK - Northern Ireland: \url{https://www.nisra.gov.uk/system/files/statistics/Weekly_Deaths\%20-\%20w\%20e\%2017th\%20June\%202022.XLSX} (accessed 2022-07-01)

UK - Scotland: https://www.nrscotland.gov.uk/files//statistics/covid19/covid-deaths-22-data-week-25.xlsx (data for 2021-2022) and https://www.nrscotland.gov.uk/files//statistics/covid19/covid-deaths-20-data-final.xlsx (data for 2022) (accessed 2022-07-01)

UK - Wales: \url{https://public.tableau.com/app/profile/public.health.wales.health.protection/viz/COVID-19Rapidmortalitydata/Summary} (accessed 2022-07-04)

Population data:

\url{https://ec.europa.eu/eurostat/databrowser/view/DEMO_R_PJANGRP3/default/table?lang=en&category=reg.reg_dem.reg_dempoar} (accessed 2022-06-22)
}


\examples{
data(C19dNUTSdata)

# Summary:
summary(C19dNUTSdata)

# Check for normal distribution:
hist(C19dNUTSdata$C19deaths_per100000)
shapiro.test(C19dNUTSdata$C19deaths_per100000)

# no. of regions for each country:
table(C19dNUTSdata$CNTR_CODE)
# only for countries with data on at least NUTS 2 level:
table(C19dNUTSdata[C19dNUTSdata$NUTS_Level > 1,]$CNTR_CODE)

}

\keyword{datasets}