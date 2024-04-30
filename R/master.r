## evOsmose-Master
## execute the full study in parallel
## M S Kapur Feb 2024 ++
## maia.kapur@noaa.gov

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

## Create input data ----
## Only need to do this if the simulations themselves or sampling protocol have changed
## (e.g., obs error, survey array, years or frequency thereof, etc)

## for one species-replicate combo, four scenarios takes about 20 seconds
# foreach(ff = 1:nrow(F0_sim)) %:%
#   #   foreach(s = 1) %dopar%
# foreach(scenario=1:4) %dopar% {
#   for(species in c(sppLabs2$Var3[sppLabs2$Var4]+1)) {


foreach(scenario=1:4)%:%
  foreach(species = c(sppLabs2$Var3[sppLabs2$Var4]+1)) %dopar%  {
    invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

    build_Data(scenario,
               sppIdx = species,
               repID = 1, ## replicates 1:29
               yrs_use = 2010:2080, ## years to extract data for
               srv_selex = 11, ## age at 50% selex
               obs_error = 0.2, ## observation error for surveys
               units = 'biomass',
               units_scalar = 1,
               do_GAM = FALSE)
  } ## end nested scenario loop

# When you're done, clean up the cluster
stopImplicitCluster()
# stopCluster()

## Run WHAM models ----
files_to_run <-   list.dirs(path = here::here('outputs','wham_runs'),recursive = FALSE)

foreach(files_to_run) %dopar% {
  invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets
  run_WHAM(yrs_use = 2010:2099, ## years to run the assessment
           file_suffix = files_to_run)
} ## end files loop


## Summarize results across all species, scenarios, simulations ----
