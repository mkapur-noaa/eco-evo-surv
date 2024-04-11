## evOsmose-Master
## execute the full study in parallel

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

#1,7,9,13

build_Data(scenario = 2,
sppIdx = 9, 
replicate = 1,
yrs_use = 2010:2044)

foreach(scenario=1:4) %dopar% {
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

build_Data(scenario,
sppIdx = 1,
 repuse = 1,
 yrs_use = 2010:2044, 
 obs_error = 0.2, ## whether or not survey biomass has observation error
 units = 'numbers',
  units_scalar = 1e6)
  #
  #sim <- replicate(1,  {retry_function(function()run_MSE(nprojyrs, neqyrs,spaceID), 1000)})
}
# When you're done, clean up the cluster
stopImplicitCluster()
# stopCluster()
