## evOsmose-Master
## execute the full study in parallel

rm(list = ls()) ## clear workspace
lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source) ## load all functions and presets

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

parallelized_sim <- foreach(replic = reps_run, .combine = cbind) %dopar% {
  #source(here::here("R",'load_files_EM.R'))
  #expID = z # HCRS 1 = status quo, 2 = eHCR, 3 = mean catch x fleet 2015-2019
  #spaceID = s  # spatial structure
  ## 1 = like OM, 
  ## 2 = 6 subareas, stocks = mgmt, with movement 
  ## 3 = 3 subareas = stocks = mgmt, with movement
  ## 4 = 3 subareas = stocks = mgmt, no movement
  ## 5 = 1 subarea, stock, mgmt (panmixia)
  #sim <- replicate(1,  {retry_function(function()run_MSE(nprojyrs, neqyrs,spaceID), 1000)})
}
# When you're done, clean up the cluster
stopImplicitCluster()
# stopCluster()
