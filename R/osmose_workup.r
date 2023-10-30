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
spp_of_interest <- c(6,8,0,12) ## sole, cod, herring, dab
spp_of_interest_idx <- spp_of_interest+1 ## for indexing into arrays
spp.labs <- c("Common Sole", "Atlantic Cod", "Atlantic Herring", "Common Dab")
names(spp.labs) <- c("6","8","0","12") 

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
biomass <- ncvar_get(ncin,"biomass")[,,,]  %>% reshape2::melt(.) %>%
filter(Var3 %in% spp_of_interest_idx)
names(biomass) <- c('lat','lon','species','timestep','biomass')
biomass$species <- biomass$species-1 ## revert to original ids, so it works with labels


biomass%>%
group_by(species) %>%
filter(timestep > 50) %>%
mutate(biomass_rescale =  rescale(biomass, to=c(0,1))) %>%
summarise(sd(biomass_rescale, na.rm =T))

biomass%>%
group_by(timestep, species) %>%
mutate(biomass_rescale =  rescale(biomass, to=c(0,1))) %>%
ungroup() %>%
filter(timestep %in% 50:54) %>%

ggplot(., aes(x = lat, y = lon, fill = biomass_rescale))+
theme_void()+
geom_raster()+
scale_fill_viridis_c(na.value = NA)+
theme(strip.text = element_text(size = 25), legend.position = 'none')+
facet_grid(species~timestep , labeller = labeller(species = spp.labs))

#facet_grid(~species , labeller = labeller(species = spp.labs))
ggsave(last_plot(),file=here('figs','Simu0_5y_biomass.png'),
height = 12, width = 12, unit = 'in',dpi = 520)


## experimental sampler----
## there are 1300 cells
## key question is a) how many stations and b) what to do when we know it doesn't cover survey extent?
## use full extent for now
n_cells <- length(unique(biomass$lat)) *  length(unique(biomass$lon))
domain <- expand.grid(unique(biomass$lat),unique(biomass$lon))
n_stations <- 0.1*n_cells
station_id <- sample(1:n_cells,n_stations) ##  random sampling sites within domain (don't change)
stations <- domain[station_id,] ## lat and lon of sampling sites
total_area <- length(domain) ## the number of cells; assume each cell = 1



biomass %>%
merge(., station_id) %>% head()
filter(~lat == stations$Var1 & ~lon == stations$Var2)

## get biomass @ stations
sample_raw <- biomass[biomass[,"lat"]==stations[,'Var1'] & biomass[,"lon"]==stations[,'Var2'],] 
## confirm that we have a record for each spp
sample_raw %>% group_by(species) %>% summarise(sum(!is.na(biomass)))

sample_raw %>%
group_by(species, timestep) %>%
mutate(sample_obs_error= rnorm(n_stations, mean = biomass, sd = 0.2) )

## apply obs error
sample_obs_error <- rnorm(n_stations, mean = sample_raw$biomass, sd = 0.2)
#sim_sample <- rescale(sim_sample$x.data, to=c(1,2)) ## do away with negs
sim_mean <- mean(sample_obs_error) * total_area ## observed index
sim_sd <- sqrt(var(sample_obs_error)*total_area^2 / n_stations)
sim_cv <- sim_sd/sim_mean

rel_err_mean <- 100*(sim_mean-true_mean)/true_mean
rel_err_cv <- 100*(sim_cv-true_cv)/true_cv