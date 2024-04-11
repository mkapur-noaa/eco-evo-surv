## build_data.R
## script to take outputs from evOSMOSE model runs and build a data frame for analysis in WHAM
## UNIT assumption is single scenario, single species, single replicate
## also save some intermediate summaries and figures
## M S Kapur 2024-04-11


## scenario: 1 = noCC_noEvo, 2 = noCC_Evo, 3 = CC_noEvo, 4 = CC_Evo
## sppIdx: 1-15
## replicate: 1-28; doesn't necessarily correspond to the file name
## survey_sampler: a fixed survey design array with years 2010-2099. Used to lookup cells.


build_Data<-function(scenario = 1, sppIdx = 1, replicate = 1,
 yrs_use = 2010:2099, units = 'numbers'){
dirtmp <- here::here('data','Ev-OSMOSE outputs',paste0('Ev-osmose_',scenLabs2[scenario,2]),'output')

## strip and format catches (yr x Age)
yield_files <- list.files(dirtmp, pattern = 'ns_yield*', recursive = T, full.names = TRUE)

## sample and build index inputs  (year, index as biom/numbers, cv, vector of ages in numbers/biomass, inputN for comps)
units_use <- ifelse(units == 'numbers','abundance','biomass') ## how the files are labeled
age_spatial_path <- list.files(dirtmp, pattern = paste0('spatial_',units_use,'byAge-',sppLabs2[sppIdx,2]), recursive = T, full.names = TRUE)[replicate]
repID <-  as.numeric(stringr::str_extract(age_spatial_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

abundance0 <- ncvar_get(nc_open(age_spatial_path),"abundance") ## this might need to switch with units_use

## this is age (26) lat (25) lon (52) timestep (24/year 2010-2099)
## redefine the timesteps into real year, week
## and then only take out the July populations
abundance <- reshape2::melt(abundance0) %>%
  mutate(year = 2010 + (Var4 - 1) %/% 24,
    month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7 & year %in% yrs_use) %>%
    group_by(year, age, lat, long, month) %>%
    summarise(value = mean(value)) %>% ## average over the month
    ungroup() %>%
    select(-month)


## biomass
## plot map of biomass
## plot raw survey data
## age comps


  ## load the data
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'.RData')))
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'_catch.RData')))
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'_biomass.RData')))
}