\encoding{UTF-8}
\name{RKI_Corona_counties}
\docType{data}
\alias{RKI_Corona_counties}
\title{German Counties with COVID-19 Cases}
\description{
  A dataset containing German counties (NUTS 3) with COVID-19 cases (Source: Robert Koch Institute).
}
\usage{
data(RKI_Corona_counties)
}
\format{
  A data.frame with multiple columns:
  \describe{
    \item{OBJECTID}{unknown/not necessary}
    \item{ADE}{(unknown/not necessary}
    \item{GF}{unknown/not necessary}
    \item{BSG}{unknown/not necessary}
    \item{RS}{(character) County code 1}
    \item{AGS}{(character) County code 2}
    \item{SDV_RS}{(character) County code 3}
    \item{GEN}{(character) County name}
    \item{BEZ}{(character) County type}
    \item{IBZ}{unknown/not necessary}
    \item{BEM}{unknown/not necessary}
    \item{NBD}{unknown/not necessary}
    \item{SN_L}{unknown/not necessary}
    \item{SN_R}{unknown/not necessary}
    \item{SN_K}{unknown/not necessary}
    \item{SN_V1}{unknown/not necessary}
    \item{SN_V2}{unknown/not necessary}
    \item{SN_G}{unknown/not necessary}
    \item{FK_S3}{unknown/not necessary}
    \item{NUTS}{(character) NUTS 3 code}
    \item{RS_0}{unknown/not necessary}
    \item{AGS_0}{unknown/not necessary}
    \item{WSK}{unknown/not necessary}
    \item{EWZ}{(numeric) Population}
    \item{KFL}{(numeric) Area in sq. km}
    \item{DEBKG_ID}{unknown/not necessary}
    \item{Shape__Are}{unknown/not necessary}
    \item{Shape__Len}{unknown/not necessary}
    \item{death_rate}{}
    \item{cases}{(numeric) COVID-19 cases}
    \item{deaths}{(numeric) COVID-19 associated deaths}
    \item{cases_per_}{(numeric) COVID-19 cases per 100,000 inhabitants}
    \item{cases_pe_1}{unknown/not necessary}
    \item{BL}{(character) Federal state}
    \item{BL_ID}{(integer) Federal state ID}
    \item{county}{(character) County name}
    \item{last_updat}{Date of last update}
    \item{geometry}{Geometry}
  }
}
\details{
The data is included as it was published in by the Robert Koch Institute (Robert Koch-Institut, RKI) but extended by the geometry column (Original data: shapefile).
}
\source{
  RKI (2020) RKI Corona Landkreise. Robert Koch-Institut (RKI), dl-de/by-2-0. Attribution: Robert Koch-Institut, Bundesamt für Kartographie und Geodäsie. \url{https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/917fc37a7095} (retrieved 2020-03-30)
}
\examples{
data(RKI_Corona_counties)
# German counties (Source: Robert Koch Institute)

Corona_nbstat <- 
  nbstat (
    RKI_Corona_counties, 
    ID_col="AGS",
    link_data = RKI_Corona_counties, 
    data_ID_col = "AGS", 
    data_col = "EWZ", 
    func = "sum"
  )
Corona_nbstat$nbmat_data_aggregate
# Sum of population (EWZ) of neighbouring counties
}