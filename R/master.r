## evOsmose-Master
## execute the full study in parallel
## M S Kapur Feb 2024 ++
## maia.kapur@noaa.gov

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets


## Create input data ----
## Only need to do this if the simulations themselves or sampling protocol have changed
## (e.g., obs error, survey array, years or frequency thereof, etc)

#* brute ----
for(scenario in 1){
  for(species in  c(sppLabs2$Var3[sppLabs2$Var4]+1)){
    for(replicate in 1:3){
      build_Data(scenario,
                 sppIdx = species,
                 repID = replicate,
                 yrs_use = 2010:2080, ## years to extract data for
                 srv_selex = 11, ## age at 50% selex
                 obs_error = 0.2, ## observation error for surveys
                 units = 'biomass',
                 units_scalar = 1,
                 do_GAM = FALSE)
    } ## end replicate
  } ## end species
} ## end scenario
#* parallel ----
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)
## for one species-replicate combo, four scenarios takes about 20 seconds
foreach(scenario=1:4) %:%
  foreach(species = c(sppLabs2$Var3[sppLabs2$Var4]+1))%:%
  foreach(replicate=1:5)  %dopar%  {
    invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

    scen_use = scenario; sp_use = species; replicate_use = replicate

    build_Data(scenario=scen_use,
               sppIdx = sp_use,
               repID = replicate_use,
               yrs_use = 2010:2080, ## years to extract data for
               srv_selex = 11, ## age at 50% selex
               obs_error = 0.2, ## observation error for surveys
               units = 'biomass',
               units_scalar = 1,
               do_GAM = FALSE)
  } ## end nested loop
# When you're done, clean up the cluster
stopImplicitCluster();stopCluster()

## Run WHAM model(s) ----
## list all the folders with outputs; can grep() or select from here
files_to_run <-   list.dirs(path = here::here('outputs','wham_runs'),recursive = FALSE)

foreach(files_to_run) %dopar% {
  invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets
  run_WHAM(yrs_use = 2010:2099, ## years to run the assessment
           file_suffix = basename(files_to_run[2]))
} ## end files loop


## Summarize results across all species, scenarios, simulations ----

list.files(here::here('outputs','wham_runs'), pattern = 'true_biomass.csv', recursive = TRUE, full.names = TRUE)[1:10] %>%
  lapply(., FUN = read.csv) %>% bind_rows() %>% head()
