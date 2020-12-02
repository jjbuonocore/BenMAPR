load(file = "healthdata_full.RData")
load("2018 pop by age processed.RData")

#===============#
#### BenMAPR ####
#===============#

#### build a map suitable for cartography ####

# this file (us.map.sf) is a low resolution file for mapping
us.map.spdf <- readOGR(dsn = "cb_2018_us_county_20m/", layer = "cb_2018_us_county_20m")
us.map.sf <- st_as_sf(us.map.spdf)

us.map.sf <- us.map.sf %>% 
  inner_join(fips_codes, by = c("STATEFP" = "state_code", "COUNTYFP" = "county_code")) %>% 
  filter(state_name != "Alaska", state_name != "Hawaii", state_name != "Puerto Rico")

#==================#
#### data time! ####
#==================#

us_tract_sf2 <- us.tract %>% unite("stco_code", sep = "", STATEFP, COUNTYFP) %>% left_join(morbidity_full2, by = c("stco_code" = "fips")) %>% left_join(mr, by = c("stco_code" = "fips")) #us_tract_sf comes from Census Data v4, mr and morbidity_full come from health data

us_tract_sf2 <- st_transform(us_tract_sf2, "+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97") # reproject the tract sf

# ggplot() + 
#   geom_raster(aes(x = x_lcc, y = y_lcc, fill = pm), data = aqdat.annual.proj) + scale_fill_gradient(low = hsv(0.6, 0.3, 0.9), high = hsv(0, 1, 0.6), trans = "log10") +
#   geom_sf(data = st_transform(us.map.sf, "+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97"), alpha = 0, color = "black", size = 0.1)
# 
# Sys.time() - tstamp
#todo - set this up somewhere in the looping

us_tract_sf2 <- us_tract_sf2 %>% 
  mutate(tractarea = st_area(us_tract_sf2)) %>% 
  mutate(tract.popinfants.dens = pop_infants / as.numeric(tractarea),
         tract.pop0to4.dens = pop_0to4 / as.numeric(tractarea),
         tract.pop5to17.dens = pop_5to17 / as.numeric(tractarea),
         tract.pop18andunder.dens = pop18andunder / as.numeric(tractarea),
         tract.pop18andover.dens = pop18andover / as.numeric(tractarea),
         tract.pop25plus.dens = pop25plus / as.numeric(tractarea),
         tract.pop65plus.dens = pop65plus / as.numeric(tractarea))

# actually do the overlay
tstamp <- Sys.time()
cmaq_tract_int <- st_intersection(cmaq.sf, us_tract_sf2)

# write_csv(cmaq_tract_int, "tract to grid intersection.csv")

# cmaq_tract_int %>% select(x_lcc, y_lcc, LAT, LON, COL, ROW, GEOID) %>% st_write(., "tract to grid intersection.shp", update = T)

Sys.time() - tstamp

cmaq_tract_int <- cmaq_tract_int %>% 
  mutate(polyarea = st_area(cmaq_tract_int)) %>% 
  mutate(poly.popinfants.dens = tract.popinfants.dens,
         poly.pop0to4.dens = tract.pop0to4.dens,
         poly.pop5to17.dens = tract.pop5to17.dens,
         poly.pop18under.dens = tract.pop18andunder.dens,
         poly.pop18over.dens = tract.pop18andover.dens,
         poly.pop25over.dens = tract.pop25plus.dens,
         poly.pop65over.dens = tract.pop65plus.dens)

cmaq_tract_int_tbl <- cmaq_tract_int %>% as_tibble() %>% select(-geometry)

#### actual health calculation ####

health.calculation <- function(cmaqdata, cmaq_tract_int_tbl) {
  
  tract.res <- cmaq_tract_int_tbl %>% left_join(cmaqdata, by = c("COL", "ROW")) %>% 
    mutate(deathspm_mid = as.numeric(polyarea) * poly.pop25over.dens * 0.0129 * pm * mr.25over, # Vodonos
           deathspm_low = as.numeric(polyarea) * poly.pop25over.dens * 0.0109 * pm * mr.25over,
           deathspm_high = as.numeric(polyarea) * poly.pop25over.dens * 0.015 * pm * mr.25over,
           
           deathso3d_mid = as.numeric(polyarea) * poly.pop25over.dens * 0.0004 * o3 * mr.25over, # Levy & Penn
           deathso3d_low = as.numeric(polyarea) * poly.pop25over.dens * 0.0002 * o3 * mr.25over,
           deathso3d_high = as.numeric(polyarea) * poly.pop25over.dens * 0.0006 * o3 * mr.25over,
           
           deathso3a_mid = as.numeric(polyarea) * poly.pop25over.dens * 0.002 * o3 * mr.25over, # Turner
           deathso3a_low = as.numeric(polyarea) * poly.pop25over.dens * 0.001 * o3 * mr.25over,
           deathso3a_high = as.numeric(polyarea) * poly.pop25over.dens * 0.004 * o3 * mr.25over,
           
           deathsno2faustini_mid = as.numeric(polyarea) * poly.pop25over.dens * 0.004 * no2 * (1.88) * mr.25over, # Faustini
           deathsno2faustini_low = as.numeric(polyarea) * poly.pop25over.dens * 0.002 * no2 * (1.88) * mr.25over,
           deathsno2faustini_high = as.numeric(polyarea) * poly.pop25over.dens * 0.006 * no2 * (1.88) * mr.25over,
           # # 
           # deathsno2eum_mid = as.numeric(polyarea) * poly.pop65over.dens * 0.0052 * no2 * mr.65over, # Eum
           # deathsno2eum_low = as.numeric(polyarea) * poly.pop65over.dens * 0.0051 * no2 * mr.65over,
           # deathsno2eum_high = as.numeric(polyarea) * poly.pop65over.dens * 0.0054 * no2 * mr.65over,
           # 
           # deathsno2atkinson_mid = as.numeric(polyarea) * poly.pop25over.dens * 0.0011 * no2 * mr.25over, # Atkinson
           # deathsno2atkinson_low = as.numeric(polyarea) * poly.pop25over.dens * 0.00053 * no2 * mr.25over,
           # deathsno2atkinson_high = as.numeric(polyarea) * poly.pop25over.dens * 0.0015 * no2 * mr.25over,
           
           cvhosppm_mid = as.numeric(polyarea) * poly.pop65over.dens * 0.00094 * pm * cvhosp_rate65over, # old CPP
           cvhosppm_low = as.numeric(polyarea) * poly.pop65over.dens * (0.00094 + (0.00015 * 1.96)) * pm * cvhosp_rate65over,
           cvhosppm_high = as.numeric(polyarea) * poly.pop65over.dens * (0.00094 - (0.00015 * 1.96)) * pm * cvhosp_rate65over,
           
           rsphosppm_mid = as.numeric(polyarea) * poly.pop65over.dens * 0.0011 * pm * resphosp_rate65over, # old CPP
           rsphosppm_low = as.numeric(polyarea) * poly.pop65over.dens * (0.0011 + (0.00027 * 1.96)) * pm * resphosp_rate65over,
           rsphosppm_high = as.numeric(polyarea) * poly.pop65over.dens * (0.0011 - (0.00027 * 1.96)) * pm * resphosp_rate65over,
           
           amipm_mid = as.numeric(polyarea) * poly.pop18over.dens * 0.0025 * pm * ami_rate18over, # Mustafic
           amipm_low = as.numeric(polyarea) * poly.pop18over.dens * 0.0015 * pm * ami_rate18over,
           amipm_high = as.numeric(polyarea) * poly.pop18over.dens * 0.0036 * pm * ami_rate18over,
           
           amino2_mid = as.numeric(polyarea) * poly.pop18over.dens * 0.0011 * 1.88 * no2 * ami_rate18over, # Mustafic
           amino2_low = as.numeric(polyarea) * poly.pop18over.dens * 0.0006 * 1.88 * no2 * ami_rate18over,
           amino2_high = as.numeric(polyarea) * poly.pop18over.dens * 0.0016 * 1.88 * no2 * ami_rate18over,
           
           resphospo3_mid = as.numeric(polyarea) * poly.pop65over.dens * 0.0016 * o3 * resphosp_rate65over, # Ji - paper gives % change / 10 ppb, this is converted to change per 1 ppb
           resphospo3_low = as.numeric(polyarea) * poly.pop65over.dens * 0.00058 * o3 * resphosp_rate65over,
           resphospo3_high = as.numeric(polyarea) * poly.pop65over.dens * 0.00263 * o3 * resphosp_rate65over,
           
           asthmapminc0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * 0.0296 * pm * asthmaincrate_0to4, # Khreis
           asthmapminc0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * 0.00995 * pm * asthmaincrate_0to4,
           asthmapminc0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * 0.0488 * pm * asthmaincrate_0to4,
           
           asthmano2inc0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * 0.0122 * 1.88 * no2 * asthmaincrate_0to4, # Khreis
           asthmano2inc0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * 0.00495 * 1.88 * no2 * asthmaincrate_0to4,
           asthmano2inc0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * 0.0169 * 1.88 * no2 * asthmaincrate_0to4,
           
           asthmapminc5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * 0.0296 * pm * asthmaincrate_5to17, # Khreis
           asthmapminc5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * 0.00995 * pm * asthmaincrate_5to17,
           asthmapminc5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * 0.0488 * pm * asthmaincrate_5to17,
           
           asthmano2inc5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * 0.0122 * 1.88 * no2 * asthmaincrate_5to17, # Khreis
           asthmano2inc5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * 0.00495 * 1.88 * no2 * asthmaincrate_5to17,
           asthmano2inc5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * 0.0169 * 1.88 * no2 * asthmaincrate_5to17,
           
           asthmapmha0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.022) / 10) * pm * asthmaharate_0to4, # Orellano
           asthmapmha0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * (log(1) / 10) * pm * asthmaharate_0to4,
           asthmapmha0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.045) / 10) * pm * asthmaharate_0to4,
           
           asthmano2ha0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.040) / 10) * no2 * asthmaharate_0to4, # Orellano
           asthmano2ha0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.001) / 10) * no2 * asthmaharate_0to4,
           asthmano2ha0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.081) / 10) * no2 * asthmaharate_0to4,
           
           asthmapmha5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.022) / 10) * pm * asthmaharate_5to17, # Orellano
           asthmapmha5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * (log(1) / 10) * pm * asthmaharate_5to17,
           asthmapmha5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.045) / 10) * pm * asthmaharate_5to17,
           
           asthmano2ha5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.040) / 10) * no2 * asthmaharate_5to17, # Orellano
           asthmano2ha5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.001) / 10) * no2 * asthmaharate_5to17,
           asthmano2ha5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.081) / 10) * no2 * asthmaharate_5to17,
           
           asthmapmed0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.022) / 10) * pm * asthmaedrate_0to4, # Orellano
           asthmapmed0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * (log(1) / 10) * pm * asthmaedrate_0to4,
           asthmapmed0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.045) / 10) * pm * asthmaedrate_0to4,
           
           asthmano2ed0to4_mid = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.040) / 10) * no2 * asthmaedrate_0to4, # Orellano
           asthmano2ed0to4_low = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.001) / 10) * no2 * asthmaedrate_0to4,
           asthmano2ed0to4_high = as.numeric(polyarea) * poly.pop0to4.dens * (log(1.081) / 10) * no2 * asthmaedrate_0to4,
           
           asthmapmed5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.022) / 10) * pm * asthmaedrate_5to17, # Orellano
           asthmapmed5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * (log(1) / 10) * pm * asthmaedrate_5to17,
           asthmapmed5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.045) / 10) * pm * asthmaedrate_5to17,
           
           asthmano2ed5to17_mid = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.040) / 10) * no2 * asthmaedrate_5to17, # Orellano
           asthmano2ed5to17_low = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.001) / 10) * no2 * asthmaedrate_5to17,
           asthmano2ed5to17_high = as.numeric(polyarea) * poly.pop5to17.dens * (log(1.081) / 10) * no2 * asthmaedrate_5to17,
           
           lbwpm_mid = as.numeric(polyarea) * poly.popinfants.dens * 0.0086 * pm * lbw_rate, # Sun
           lbwpm_low = as.numeric(polyarea) * poly.popinfants.dens * 0.0031 * pm * lbw_rate,
           lbwpm_high = as.numeric(polyarea) * poly.popinfants.dens * 0.014 * pm * lbw_rate,
           
           ptbpm_mid = as.numeric(polyarea) * poly.popinfants.dens * 0.012 * pm * ptb_rate, # Sun
           ptbpm_low = as.numeric(polyarea) * poly.popinfants.dens * 0.0029 * pm * ptb_rate,
           ptbpm_high = as.numeric(polyarea) * poly.popinfants.dens * 0.021 * pm * ptb_rate,
           
           asdpmlam_mid = as.numeric(polyarea) * poly.popinfants.dens * 0.084 * pm * asd_rate, # Lam
           asdpmlam_low = as.numeric(polyarea) * poly.popinfants.dens * 0.076 * pm * asd_rate,
           asdpmlam_high = as.numeric(polyarea) * poly.popinfants.dens * 0.092 * pm * asd_rate,
           
           asdpmbecerra_mid = as.numeric(polyarea) * poly.popinfants.dens * 0.01445 * pm * asd_rate, # Becerra
           asdpmbecerra_low = as.numeric(polyarea) * poly.popinfants.dens * 0 * pm * asd_rate,
           asdpmbecerra_high = as.numeric(polyarea) * poly.popinfants.dens * 0.02986 * pm * asd_rate
           
    )
  return(tract.res)
}

healthres.raw <- map(cmaq, health.calculation, cmaq_tract_int_tbl)
