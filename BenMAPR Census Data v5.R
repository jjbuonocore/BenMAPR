library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(broom)
library(maps)
library(rmapshaper)
library(sf)
library(tidyverse)
library(RCurl)
library(tidycensus)
library(tigris)

#### Get Variable Listing files from US Census ####

# censuskey <- put your Census API key here

# test <- get_acs("county subdivision", product = "population", year = 2017, state = "CA", key = censuskey, table = "S0101")

acs5s <- tidycensus::load_variables(2018, "acs5/subject", cache = T)
acs5p <- load_variables(2018, "acs5/profile", cache = T)
acs5 <- load_variables(2018, "acs5", cache = T)

#### Get relevant states ####

data(fips_codes)

fips_codes <- fips_codes %>% unite(col = "stco_code", state_code, county_code, sep = "", remove = F)

fips_codes %>% distinct(state_name)

fips.benmapr <- fips_codes %>% 
  filter(!(state_name %in% c("Hawaii", "Alaska", "American Samoa", "Guam", "Northern Mariana Islands", "Puerto Rico", "U.S. Minor Outlying Islands", "U.S. Virgin Islands")))

benmapr.states <- fips.benmapr %>% distinct(state_name) %>% pull(state_name)

tcistates <- c("Maine", "Vermont", "New Hampshire", "New York", "Massachusetts", "Connecticut", "Rhode Island", "Pennsylvania", "Maryland", "Delaware", "District of Columbia", "New Jersey", "Virginia")
tcistatesabb <- c("ME", "VT", "NH", "NY", "MA", "CT", "RI", "PA", "MD", "DE", "DC", "NJ", "VA")

#### Filter variable lists from US Census ####

# just age

# acs5s %>% filter(concept == "AGE AND SEX") %>% 
#   separate(label, c("x1", "x2", "x3", "x4", "agegroup"), sep = "!!") %>% View() # just looking at the full list

popagevarlist <- acs5s %>% filter(concept == "AGE AND SEX") %>% 
  separate(label, sep = "!!", into = c("x1", "x2", "popvar", "age", "agegroup")) %>% 
  filter(x2 == "Total" & age == "AGE") # grab variables for total population by age

totalpopvar <- acs5s %>% filter(concept == "AGE AND SEX") %>% 
  separate(label, sep = "!!", into = c("x1", "x2", "popvar", "age", "agegroup")) %>% 
  filter(x2 == "Total" & is.na(age)) # grab total population

#### Download Census data from API ####

us_tract_sf_tot <- map_dfr(
  benmapr.states, 
  ~ get_acs("tract", variables = totalpopvar$name, year = 2018, geometry = T, key = censuskey, state = ., keep_geo_vars = T)
)
# filter to needed variables
us_tract_sf_tot <- us_tract_sf_tot %>% select(-AFFGEOID, -NAME.x, -LSAD, -ALAND, -AWATER, -moe) %>% rename(NAME = NAME.y)

us_tract_df <- map_dfr(
  benmapr.states, 
  ~ get_acs("tract", variables = popagevarlist$name, year = 2018, geometry = F, key = censuskey, state = .)
)

gc()

# broad plan is to 1) start with census raw data, 2) join in variable definitions, 3) do pivot_wider() to convert to wide format, 4) do a mutate function to get sums in the bigger age bins (below), 5) join to the tract spdf, and 6) convert to sf

# map of total population for a state

us_tract_sf_tot <- us_tract_sf_tot %>% select(-variable) %>% rename(pop_allages = estimate)
ggplot(us_tract_sf_tot %>% filter(str_detect(NAME, "Pennsylvania"))) + geom_sf(aes(fill = pop_allages), size = 0.01)

#### work with data just broken down by age ####

# this joins in the variable names, widens the dataset, and uses the names of the age bins as the new variable names. Original census variable code gets dropped.
us.tract.age <- us_tract_df %>% left_join(., popagevarlist, by = c("variable" = "name")) %>% select(-x1, -x2)
us.tract.age.wide <- us.tract.age %>% select(-concept, -variable, -age, -popvar, -moe) %>% pivot_wider(id_cols = c(GEOID, NAME), names_from = agegroup, values_from = estimate) 

# aggregate from census age bins to CRF-relevant age bins. pop_total_calculated serves as a check
AgeAggregation <- function(uscensusdata.widebyage) {
  res <- uscensusdata.widebyage %>% 
    mutate(pop_infants = `Under 5 years` * (1/5),
           pop_0to4 = `Under 5 years`,
           pop_5to17 = `5 to 9 years` + `10 to 14 years` + (`15 to 19 years` * (3/5)),
           pop18andunder = `Under 5 years` + `5 to 9 years` + `10 to 14 years` + (`15 to 19 years` * (4/5)),
           pop18andover = (`15 to 19 years` * (2/5)) + `20 to 24 years` + `25 to 29 years` + `30 to 34 years` + `35 to 39 years` + `40 to 44 years` + `45 to 49 years` + `50 to 54 years` + `55 to 59 years` + `60 to 64 years` + `65 to 69 years` + `70 to 74 years` + `75 to 79 years` + `80 to 84 years` + `85 years and over`,
           pop25plus = `25 to 29 years` + `30 to 34 years` + `35 to 39 years` + `40 to 44 years` + `45 to 49 years` + `50 to 54 years` + `55 to 59 years` + `60 to 64 years` + `65 to 69 years` + `70 to 74 years` + `75 to 79 years` + `80 to 84 years` + `85 years and over`,
           pop65plus = `65 to 69 years` + `70 to 74 years` + `75 to 79 years` + `80 to 84 years` + `85 years and over`,
           poptotal_calculated = `Under 5 years` + `5 to 9 years` + `10 to 14 years` + `15 to 19 years` + `20 to 24 years` + `25 to 29 years` + `30 to 34 years` + `35 to 39 years` + `40 to 44 years` + `45 to 49 years` + `50 to 54 years` + `55 to 59 years` + `60 to 64 years` + `65 to 69 years` + `70 to 74 years` + `75 to 79 years` + `80 to 84 years` + `85 years and over`) %>% 
    select(-(`Under 5 years`:`85 years and over`))
  return(res)
}

us.tract.age.wide <- AgeAggregation(us.tract.age.wide)

# join the age aggregated dataset to the shapefile, and reproject
us.tract <- geo_join(spatial_data = us_tract_sf_tot, data_frame = us.tract.age.wide, by_sp = c("GEOID", "NAME"), by_df = c("GEOID", "NAME")) %>% st_transform(., "+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97")

# us.tract.sf <- st_as_sf(us.tract) 

ggplot(us.tract %>% filter(str_detect(NAME, "Massachusetts"))) + geom_sf(aes(fill = pop65plus), size = 0.01)

ggplot(us.tract %>% left_join(fips_codes, by = c("STATEFP" = "state_code", "COUNTYFP" = "county_code")) %>% filter(state_name %in% tcistates)) + geom_sf(aes(fill = pop18andover), size = 0.01) + scale_fill_gradient(low = "white", high = "dark blue")

save(us.tract, fips_codes, fips.benmapr, benmapr.states, tcistates, tcistatesabb, censuskey, popagevarlist, totalpopvar, file = "2018 pop by age processed.RData")


