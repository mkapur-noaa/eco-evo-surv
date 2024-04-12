## evOsmose-Master
## execute the full study in parallel

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

#1,7,9,13

build_Data(scenario = 2,
sppIdx = 7, 
replicate = 1,
yrs_use = 2010:2044)

## for one species-replicate combo, four scenarios takes about 20 seconds
foreach(scenario=1:4) %dopar% {
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

for(species in c(1,7,9,13)) {
build_Data(scenario,
sppIdx = species,
 repuse = 1,
 yrs_use = 2010:2044, 
 obs_error = 0.2, ## whether or not survey biomass has observation error
 units = 'numbers',
  units_scalar = 1e6)
  #
  #sim <- replicate(1,  {retry_function(function()run_MSE(nprojyrs, neqyrs,spaceID), 1000)})
} ## end species loop
} ## end scenario loop (parallel)

# When you're done, clean up the cluster
stopImplicitCluster()
# stopCluster()
