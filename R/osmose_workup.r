## workup of evosmose outputs
#0=AtlanticHerring 1=AtlanticMackerel 2=Sandeel 3=EuropeanSprat 4=NorwayPout
# 5=EuropeanPlaice 6=CommonSole 7=Saithe 8=AtlanticCod 9=Haddock 10=HorseMackerel
##11=Whiting 12=CommonDab 13=GreyGurnard 14=Hake 15=Shrimp
#https://pjbartlein.github.io/REarthSysSci/netCDF.html#reading-restructuring-and-writing-netcdf-files-in-r
#install.packages('ncdf4')
library(ncdf4)
library(dplyr)
library(here)
library(ggplot2)
library(scales)

scenPal <- c('#d8a6a6','#a00000','#9fc8c8','#298c8c') ## scenario colors: reds have climate change, blues do not
scenLabs <- c('No Climate Change, No Evolution',
'No Climate Change, Evolution',
'Climate Change, No Evolution',
'Climate Change, Evolution')

#spp_of_interest <- c(0,1,3,5,8,11,13)
spp_of_interest <- c(8,11,0,3)
spp_of_interest_idx <- spp_of_interest+1 ## for indexing into arrays



spp.labs <-c("Atlantic Herring" , "Atlantic Mackerel" ,"Sandeel" ,"European Sprat","Norway Pout" ,
"European Plaice","Common Sole","Saithe","Atlantic Cod","Haddock","Horse Mackerel" ,
"Whiting","Common Dab","Grey Gurnard" ,"Hake","Shrimp")
names(spp.labs) <- as.character(0:15)

spp.labs[spp_of_interest_idx]

#In this file, you have some vectors and matrices (x, y, species, trophic level, size,
# abundance, age, genotype, ...).
nc_path <- here('data','ns_snapshot_step2399.nc')
ncin <- nc_open(nc_path)
print(ncin)

## extract lat, long, abundance, and age for key species and reshape
## this works for the snapshot
get_dat <- function(ncdf_file){
spp <- ncvar_get(ncin,"species") ## index spp of interest

lon <- ncvar_get(ncin,"x")[spp %in% spp_of_interest]
lat <- ncvar_get(ncin,"y")[spp %in% spp_of_interest]
age <- ncvar_get(ncin,"age")[spp %in% spp_of_interest]
abundance <- ncvar_get(ncin,"abundance")[spp %in% spp_of_interest]

tmp <- data.frame(species = spp[spp %in% spp_of_interest],
lon, lat, age, abundance)
return(tmp)
}
get_dat(ncdf_file = ncin) %>%
group_by(species) %>% summarise(median(abundance), max(abundance))



step2399 <- get_dat(ncdf_file = ncin) %>%
filter(abundance< 250000) %>% ## there are some weirdly high values offmap?
group_by(species) %>%
mutate(abund_rescale =  rescale(abundance, to=c(0,1))) %>%
ungroup()

ggplot(step2399, aes(x = lat, y = lon, fill = abund_rescale))+
geom_tile()+coord_equal() +
facet_wrap(~species, labeller = labeller(species = spp.labs))


nc_path <- here('data','ns_spatialized_Simu0.nc')
ncin <- nc_open(nc_path)
print(ncin)
names(ncin$var)
dim(ncin)

## lat, lon, spp, yr
biomass <- ncvar_get(ncin,"biomass")[,,,]  %>%
reshape2::melt(.) %>%
filter(Var3 %in% spp_of_interest_idx)
names(biomass) <- c('lat','lon','species','timestep','true_biomass')
biomass$species <- biomass$species-1 ## revert to original ids, so it works with labels
data.frame(spp.labs)

biomass$species_name <- merge(biomass, data.frame(spp.labs), by.x = 'species', by.y = 'row.names')$spp.labs
save(biomass, file = here('data','ns_spatialized_simu0_biomass.rda'))


bio_filter_5054 <- biomass%>%
group_by( species) %>%
mutate(biomass_rescale =  rescale(true_biomass, to=c(0,1))) %>%
ungroup() %>%
filter(timestep %in% 50:54)

maps <- ggplot(data = bio_filter_5054, aes(x = lat, y = lon,
fill = biomass_rescale))+
theme_void()+
geom_raster()+
scale_fill_viridis_c(na.value = NA)+
theme(strip.text = element_text(size = 25),
strip.text.y = element_blank(),
legend.position = 'none')+
facet_grid(species~timestep, labeller = labeller(species = spp.labs))

bio_filter_avg <- biomass%>%
group_by(species, lat, lon) %>%
summarise(mtb = mean(true_biomass,na.rm = TRUE)) %>%
ungroup() %>%
mutate(biomass_rescale =  rescale(mtb, to=c(0,1)))
head(bio_filter_avg)

maps <- ggplot(bio_filter_avg, aes(x = lat, y = lon,
fill = biomass_rescale))+
theme_void()+
geom_raster()+
scale_fill_viridis_c(na.value = NA)+
theme(strip.text = element_text(size = 25),
strip.text.y = element_blank(),
legend.position = 'none')+
facet_grid(species, labeller = labeller(species = spp.labs))

#facet_grid(~species , labeller = labeller(species = spp.labs))
ggsave(last_plot(),file=here('figs','Simu0_5y6spp_biomass.png'),
height = 12, width = 12, unit = 'in',dpi = 520)


## experimental sampler----
## there are 1300 cells
## key question is a) how many stations and b) what to do when we know it doesn't cover survey extent?
## use full extent for now
domain <- unique(biomass[c("lat","lon")])
n_cells <- nrow(domain)
n_stations <- 0.5*n_cells
station_id <- sample(1:n_cells,n_stations) ##  random sampling sites within domain (don't change)
stations <- domain[station_id,] ## lat and lon of sampling sites
total_area <- n_cells  ## assumes each station has an area of 1


## get biomass @ stations
sample_raw <- right_join(biomass, stations, by = c("lat","lon")) %>% filter(!is.na(true_biomass))
## confirm that we have a record for each spp  (did not sample an NA cell)
sample_raw %>% group_by(species) %>% summarise(sum(!is.na(true_biomass)))
nrow(sample_raw)/nrow(biomass) <= n_stations/n_cells ## should return correct ratio of data
write.csv(sample_raw, here('data','sampler_raw.csv'))
sample_raw %>% group_by(species,timestep) %>% summarise(n=n())

sampler_results <- merge(sample_raw %>%
group_by(species, timestep) %>%
summarise(sim_mean = mean(true_biomass)*total_area,
sim_sd = sqrt(var(true_biomass, na.rm = T)*total_area^2 / n_stations),sim_cv = sim_sd/sim_mean ) %>%
ungroup() ,
 biomass %>%
 group_by(species, timestep) %>%
summarise(true_mean = sum(true_biomass, na.rm = T),
true_sd = sqrt(var(true_biomass, na.rm = T) / n_stations),true_cv = true_sd/true_mean ) %>%
ungroup(),
by = c('species','timestep'))  %>%
mutate(rel_err_mean = 100*(sim_mean-true_mean)/true_mean,
rel_err_cv =100*(sim_cv-true_cv)/true_cv,ratio = sim_mean/true_mean) %>%
group_by(species) %>%
mutate(true_mean_rescale =  rescale(true_mean, to=c(0,1)),
sim_mean_rescale =  rescale(sim_mean, to=c(0,1)),
true_cv_rescale =  rescale(true_cv, to=c(0,1)),
sim_cv_rescale =  rescale(sim_cv, to=c(0,1))) %>%
ungroup()



## Population vs Observations, with error
rawobs <- ggplot(sampler_results, aes(x = timestep)) +
theme_minimal() +
#geom_line(aes(y = true_mean, color = 'Population')) +
geom_point(aes(y = sim_mean, color = 'Observed')) +
scale_color_manual(values = c('grey22','indianred')) +
geom_errorbar(aes(ymin = true_mean-true_mean*true_cv,
 ymax = true_mean+true_mean*true_cv,color = 'Population'), width = 0) +
geom_errorbar(aes(ymin = sim_mean-sim_mean*sim_cv,
ymax = sim_mean+sim_mean*sim_cv,color = 'Observed'), width = 0) +
annotate('rect', ymin =-Inf, ymax = Inf, xmin = 50, xmax = 55, alpha = 0.5, fill = '#B8DE29FF') +
facet_wrap(~species, scales = 'free_y',labeller =  labeller(species = spp.labs), ncol = 1)

write.csv(sampler_results, here('data','sampler_results.csv'))
## log scale
ggplot(sampler_results, aes(x = timestep)) +
theme_bw() +
geom_line(aes(y = log(true_mean), color = 'Population')) +
geom_point(aes(y = log(sim_mean), color = 'Observed')) +
#geom_errorbar(aes(ymin = log(true_mean-true_mean*true_cv),
# ymax = true_mean+true_mean*true_cv,color = 'Population'), width = 0) +
#geom_errorbar(aes(ymin = sim_mean-sim_mean*sim_cv,
#ymax = sim_mean+sim_mean*sim_cv,color = 'Observed'), width = 0) +
facet_wrap(~species, scales = 'free_y',labeller =  labeller(species = spp.labs))

## scaled to mean
mre <- sampler_results %>%
ggplot(., aes(x = timestep)) +
theme_minimal() +
scale_color_manual(values = c('grey22','indianred')) +
geom_line(aes(y = true_mean_rescale, color = 'Population')) +
geom_point(aes(y =sim_mean_rescale, color = 'Observed')) +
geom_errorbar(aes(ymin = sim_mean_rescale-sim_mean_rescale*sim_cv_rescale,
ymax = sim_mean_rescale+sim_mean_rescale*sim_cv_rescale,color = 'Observed'), width = 0) +
annotate('rect', ymin = 0, ymax = Inf, xmin = 50, xmax = 55, alpha = 0.5, fill = '#B8DE29FF') +
facet_wrap(~species,labeller =  labeller(species = spp.labs), ncol = 1) +
labs(x = 'Year', y = 'Abundance, Scaled to Mean', color = '') +

theme(legend.position = 'bottom')


ggsave(last_plot(),file=here('figs','scaled_obs_0.5coverage.png'),
height = 12, width = 12, unit = 'in',dpi = 520)


## scaled to mean rel error
 sampler_results %>%
mutate(sim_lower_rescale =sim_mean_rescale-sim_mean_rescale*sim_cv_rescale,
 sim_upper_rescale =sim_mean_rescale+sim_mean_rescale*sim_cv_rescale) %>%
ggplot(., aes(x = timestep)) +
theme_bw() +
scale_color_manual(values = c('grey22','indianred')) +
geom_line(aes(y = (sim_mean_rescale-true_mean_rescale)/true_mean_rescale, color = 'Population')) +
geom_ribbon(aes(ymin = (sim_lower_rescale-true_mean_rescale)/true_mean_rescale,
ymax = (sim_upper_rescale-true_mean_rescale)/true_mean_rescale))+
facet_wrap(~species, labeller =  labeller(species = spp.labs), ncol = 1) +
labs(x = 'Year', y = 'Abundance, Scaled to Mean', color = '') +

theme(legend.position = 'bottom')


ggsave(last_plot(),file=here('figs','scaled_obs_0.5coverage.png'),
height = 12, width = 12, unit = 'in',dpi = 520)

png(here('figs','maps_scaled_obs_raw_0.5coverage.png'),
height = 16, width = 14, unit = 'in',res = 520)

Rmisc::multiplot(maps, mre, rawobs, cols = 3)

dev.off()


## age comps ----
nc_path <- here::here('data','ns_spatial_abundancebyAge_Simu0_withF',
'ns_spatial_abundancebyAge-Whiting_Simu0.nc')

nc_path <- here::here('data','ev-osmose outputs','ev-osmose_cc_evo','output','spatial',
'ns_spatial_abundancebyAge-AtlanticCod_Simu0.nc')

nc_path <- here::here('data','ev-osmose outputs','ev-osmose_noCC_noEvo','output','spatial',
'ns_spatial_abundancebyAge-AtlanticCod_Simu0.nc')

ncin <- nc_open(nc_path)

abundance0 <- ncvar_get(ncin,"abundance")
## this is age (26) lat (25) lon (52) timestep (24/year 2010-2099)
## redefine the timesteps into real year, week
## and then only take out the July populations
abundance <- reshape2::melt(abundance0) %>%
  mutate(year = 2010 + (Var4 - 1) %/% 24,
    month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7) %>%
    group_by(year, age, lat, long, month) %>%
    summarise(value = mean(value)) %>% ## average over the month
    ungroup() %>%
    select(-month)

head(abundance)

# Define a function to perform the sampling
sample_individuals <- function(df, timestep, long, lat) {
  # Exclude NAs
  df <- df[!is.na(df$value), ]

# If there are no positive values, return a data frame with zero counts for all age groups
  if (all(df$value <= 0)) {
    return(data.frame(timestep = timestep, long = long, lat = lat,
    age = 1:max(abundance$age), count = integer(max(abundance$age))))
  }
  # Calculate the proportions of each age group
  props <- df$value / sum(df$value)

  # Sample 500 individuals using a multinomial distribution
  sampled <- rmultinom(1, 500, props)

   # Create a data frame with zero counts for all age groups
  result <- data.frame(timestep = timestep, long = long, lat = lat,
  age = 1:max(abundance$age), count = integer(max(abundance$age)))

  # Fill in the sampled counts for the age groups that were present in the input data
  result[result$age %in% df$age, "count"] <- sampled

  # Return the number of individuals in each age group
  return(result)
    #return(data.frame(timestep = timestep, long = long, lat = lat, age = df$age, count = result))
}

# Identify unique combinations of long and lat
spatial_cells <- unique(abundance[, c("long", "lat")])
spatial_cells <- spatial_cells[!apply(spatial_cells, 1, function(cell) {
  all(is.na(abundance[abundance$long == cell["long"] & abundance$lat == cell["lat"], "value"]))
}), ]
# Initialize an empty list to store the results
results <- list()

# For each timestep
for (timestep in unique(abundance$year)) {
  cat(timestep,"\n")
  # Filter the data for the current timestep
  timestep_data <- abundance[abundance$year == timestep, ]

  # Randomly select 50% of the spatial cells
  selected_cells <- spatial_cells[sample(nrow(spatial_cells), nrow(spatial_cells) / 2), ]

  # For each selected cell
  for (i in 1:nrow(selected_cells)) {
    # Filter the data for the current cell
    cell_data <- timestep_data[timestep_data$long == selected_cells$long[i] &
                               timestep_data$lat == selected_cells$lat[i], ]

    # Perform the sampling and store the results
  results[[paste(timestep, selected_cells$long[i], selected_cells$lat[i], sep = "_")]] <-
      sample_individuals(cell_data, timestep, long = selected_cells$long[i], lat = selected_cells$lat[i])
  }
}

# Combine the results into a single data frame
results_df <- do.call(rbind, results)

# Aggregate the data by timestep and age, summing the count
results_df <- results_df %>%
  group_by(timestep, age) %>%
  summarise(count = sum(count))

# Calculate the frequency for each age group within each unique timestep
results_df <- results_df_age %>%
  group_by(year) %>%
  mutate(frequency = count / sum(count))


results_df %>% filter(age == 3)
# Create a new column for the cell identifier
#results_df$cell <- rownames(results_df)
#results_df$cell <- interaction(results_df$lat, results_df$long, sep = "_")

# Convert the timestep and cell columns to factors
##results_df$timestep <- as.factor(results_df$timestep)
#results_df$cell <- as.factor(results_df$cell)

#results_df %>% group_by(timestep, cell) %>% summarise(n = n())

# Create the plot
ggplot(results_df,
aes(x = age, y = frequency, color = timestep, group = timestep)) +
  geom_line() +
  #facet_wrap(~ timestep) +
  labs(x = "Age", y = "Frequency") +
  theme_minimal()+  theme(legend.position='none')

ggplot(results_df_1,
aes(x = age, y = frequency, color = timestep, group = timestep)) +
  geom_line() +
  #facet_wrap(~ timestep) +
  labs(x = "Age", y = "Frequency") +
  theme_minimal()+  theme(legend.position='none')

comps <- abundance %>%
filter(!is.na(value)) %>%
    tidytable::mutate(tot = sum(value),
    .by = timestep) %>%
    tidytable::summarise(numbers = mean(tot),
    age_tot = sum(value),
    .by = c(timestep, age)) %>%
    tidytable::mutate(prop = age_tot / numbers)




comps %>% filter( age == 4)

ggplot(subset(comps, age < 10), aes(x = age, y = prop)) +
geom_point()+geom_line() +
facet_wrap(~timestep)


biom <- list.files(here::here('data','ev-osmose outputs'),
pattern = "biomass_Simu0.csv", recursive = TRUE, full.names = TRUE)
basename(dirname(dirname(biom)))

dat2 <- lapply(biom, function(x) {
  dat <- read.csv(x, skip = 1) %>%
  mutate(
         year = 2010 + floor(Time-70.04166),   # Adjust the Time so that 0 corresponds to January 1
         month =cut(Time %% 1, breaks = seq(0, 1, 1/12),
         labels = 1:12, include.lowest = TRUE),
         scenario = basename(dirname(dirname(x))))


  #dat$species <- gsub('biomass_Simu0.csv','',x)
  return(dat)
}) %>% bind_rows()

## plot average annual biomass
dat2 %>%
select(-Time) %>%
reshape2::melt(id = c('year','month','scenario')) %>%
filter(!is.na(value)) %>%
group_by(year, variable, scenario) %>%
summarise(sd_v = sd(value),  mvalue = mean(value),
lci = mvalue -1.96*sd_v, uci = mvalue+1.96*sd_v) %>%
#filter(variable %in% c('NorwayPout','AtlanticHerring','AtlanticCod','AtlanticMackerel','EuropeanSprat')) %>%
ggplot(., aes(x = year, y = mvalue, color = scenario)) +
geom_line() +
#geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
scale_color_manual(values = scenPal)+
facet_wrap(~variable, scales = 'free_y') +
ggsidekick::theme_sleek() +
labs(x = 'Year', y = 'Mean Biomass (t)')

ggsave(last_plot(),file=here('figs','scenario_compare_biomass.png'),
height = 12, width = 12, unit = 'in',dpi = 520)
