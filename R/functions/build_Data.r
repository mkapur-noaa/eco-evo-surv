## build_data.R
## script to take outputs from evOSMOSE model runs and build a data frame for analysis in WHAM
## UNIT assumption is single scenario, single species, single replicate
## also save some intermediate summaries and figures
## M S Kapur 2024-04-11


## scenario: 1 = noCC_noEvo, 2 = noCC_Evo, 3 = CC_noEvo, 4 = CC_Evo
## sppIdx: 1-15
## replicate: 1-28; doesn't necessarily correspond to the file name
## yrs_use: years to extract data for
## obs_error: whether or not survey biomass has observation error (NULL or value)
## units: 'numbers' or 'biomass'
## units_scalar: optional factor to decrease the scale of the index data (applies to comps as well)


build_Data<-function(scenario = 1, 
sppIdx = 8,
 replicate = 1,
 yrs_use = 2010:2099, 
 obs_error = 0.2, ## whether or not survey biomass has observation error
 units = 'numbers',
  units_scalar = 1e6){

dirtmp <- here::here('data','Ev-OSMOSE outputs',paste0('Ev-osmose_',scenLabs2[scenario,2]),'output')

## strip and format catches (yr x Age)
yield_files <- list.files(dirtmp, pattern = 'ns_yield*', recursive = T, full.names = TRUE)

## sample and build index inputs  (year, index as biom/numbers, cv, vector of ages in numbers/biomass, inputN for comps)
units_use <- ifelse(units == 'numbers','abundance','biomass') ## how the files are labeled

age_spatial_path <- list.files(dirtmp, pattern = paste0('spatial_',units_use,'byAge-',sppLabs2[sppIdx,2]), recursive = T,
 full.names = TRUE)[replicate]
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


total_area <- length(unique(abundance$lat))*length(unique(abundance$long))


# Initialize an empty list to store the results
results_age <- results_index <- list()

# For each timestep
for (timestep in unique(abundance$year)) {
  cat(timestep,"\n")
  # Filter the data for the current timestep
  timestep_data <- abundance[abundance$year == timestep, ]

# pull out the survey stations for this year
  selected_cells <- survey_array[survey_array$year == timestep, ]

# expand the mean and sd of the abundance for the selected cells
# this approach has NO survey selectivity
survey_abund <- semi_join(timestep_data, selected_cells, by = c("lat", "long")) %>%
filter(!is.na(value)) %>% ## remove NAs
group_by(lat, long) %>%
summarise(station_abund = ifelse(is.null(obs_error), 
sum(value),
rnorm(1,sum(value),obs_error*mean(value)))) %>% ## sum over all ages
ungroup() %>%
summarise(
abund_mean = mean(station_abund)*total_area, 
abund_sd = sqrt(var(station_abund, na.rm = T)*total_area^2 / nrow(selected_cells)),
abund_cv = abund_sd/abund_mean ) %>%
mutate(year = timestep, replicate = repID, scenario = scenario, species = sppLabs2[sppIdx,2])

results_index[[paste(timestep)]] <- survey_abund


## perform multinomial sampling of ages at each station
  for (i in 1:nrow(selected_cells)) {
    # Filter the data for the current cell (all ages present in population)
    cell_data <- timestep_data[timestep_data$long == selected_cells$long[i] & 
                               timestep_data$lat == selected_cells$lat[i], ]
    # Perform the age sampling and store the results
  results_age[[paste(timestep, selected_cells$long[i], selected_cells$lat[i], sep = "_")]] <- 
      sample_ages(cell_data, timestep, long = selected_cells$long[i], lat = selected_cells$lat[i])

  } ## end stations loop for age comps
} ## end timesteps loop for survey data

# Combine the results into a single data frame
results_df_age_spatial <- do.call(rbind, results_age)
results_df_index <- do.call(rbind, results_index)

# reduce the scale of the index data 
results_df_index[1:2,] <- round(results_df_index[1:2,],2)

# Aggregate the agecomp data by timestep and age, summing the count (collapse space)
results_df_age <- results_df_age_spatial %>%
  group_by(timestep, age) %>%
  summarise(count = sum(count)) %>%
  select(year = timestep, age, count)


# rescale, reshape to WHAM format and save
survey_results <- results_df_index %>%
mutate(abund_mean= abund_mean/units_scalar)%>%
select(year, abund_mean, abund_cv) %>%
merge(., results_df_age %>% mutate(count = count/units_scalar) %>% tidyr::pivot_wider(., names_from = age, values_from = count), by = 'year') %>%
mutate(inputN = 50) 

# save the results
write.csv(survey_results,  
here::here('data','wham_inputs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_survey.csv')), 
row.names = FALSE)


## summary figures

## maps of true abundance thru time
map <- abundance %>%
group_by(year, lat, long) %>%
summarise(tot_val = sum(value)) %>%
mutate(abundance_rescale =  rescale(tot_val, to=c(0,1), na.rm = T)) %>% 
ungroup() %>%
filter(year %in% seq(2015,2075,length.out = 4))  %>% 
ggplot(data = ., aes(x = lat, y = long, 
fill = abundance_rescale))+
theme_void()+
geom_raster()+
scale_fill_viridis_c(na.value = NA)+
theme(strip.text = element_text(size = 25), 
strip.text.y = element_blank(),
legend.position = 'none')+
facet_wrap(~year, ncol = 2)


## survey index data
true_abund <- abundance %>% 
group_by(year) %>% 
summarise(abund_mean=sum(value,na.rm = T)/units_scalar) %>% 
ungroup()


index <- ggplot(survey_results, aes(x = year, y = abund_mean)) +
geom_point()+
#geom_line(data = true_abund, color = 'red')+
geom_errorbar(aes(ymin = abund_mean - abund_cv*abund_mean, 
ymax = abund_mean + abund_cv*abund_mean), width = 0) +
scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 5))+
labs(x = 'Year', y = paste0('Abundance ',
ifelse(units == 'numbers','(millions)','(tons)'))) 

## survey age comps
comps <- results_df_age %>%
group_by(year) %>%
mutate(frequency = count / sum(count)) %>%
ungroup() %>%
ggplot(., aes(x = age, y = frequency, group = year, color = factor(year))) +
geom_line()+ 
scale_x_continuous(breaks = seq(0, max(results_df_age$age), by = 5))+
labs(x = 'Age', y = 'Frequency', color = 'Year') + 
scale_color_manual(values =  monochromeR::generate_palette("purple", modification = "go_lighter", 
                 n_colours = length(yrs_use), view_palette = TRUE))+
theme(legend.position = 'none')



png(file = here::here('figs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',
scenLabs2[scenario,2],'-abundance-',units,'-',Sys.Date(),'.png')),
height = 8, width = 12, unit = 'in',res = 520)

Rmisc::multiplot(map, index, comps, cols = 3)

dev.off()


ggsave(last_plot(),
file=here::here('figs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',
scenLabs2[scenario,2],'-abundance-',units,'-',Sys.Date(),'.png')),
height = 4, width = 8, unit = 'in',dpi = 520)



# Calculate the frequency for each age group within each unique timestep




## biomass
## plot map of biomass
## plot raw survey data
## age comps


  ## load the data
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'.RData')))
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'_catch.RData')))
  #load(here::here('data','evosmose_outputs',paste0('evosmose_',scenario,'_biomass.RData')))
}