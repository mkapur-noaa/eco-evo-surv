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
#spp_of_interest <- c(0,1,3,5,8,11,13) 
spp_of_interest <- c(8,11,0,3)
spp_of_interest_idx <- spp_of_interest+1 ## for indexing into arrays

spp.labs <-c("Atlantic Herring" , "Atlantic Mackerel" ,"Sandeel" ,"European Sprat","Norway Pout" ,
"European Plaice","Common Sole","Saithe","Atlantic Cod","Haddock","Horse Mackerel" ,
"Whiting","Common Dab","Grey Gurnard" ,"Hake","Shrimp")
names(spp.labs) <- as.character(0:15)

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

nc_path <- here::here('data','ns_spatial_abundancebyAge_Simu0_withF',
'ns_spatial_abundancebyAge-AtlanticCod_Simu0.nc')
ncin <- nc_open(nc_path)
get_dat(ncdf_file = ncin) %>% head()

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


## age comps
nc_path <- here('data','ns_spatialized_Simu0.nc')
ncin <- nc_open(nc_path)
print(ncin)
names(ncin$var)
dim(ncin)
