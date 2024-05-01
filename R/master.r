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
# for(scenario in 1){
#   for(species in  c(sppLabs2$Var3[sppLabs2$Var4]+1)){
#     for(replicate in 1:3){
#       build_Data(scenario,
#                  sppIdx = species,
#                  repID = replicate,
#                  yrs_use = 2010:2080, ## years to extract data for
#                  srv_selex = 11, ## age at 50% selex
#                  obs_error = 0.2, ## observation error for surveys
#                  units = 'biomass',
#                  units_scalar = 1,
#                  do_GAM = FALSE)
#     } ## end replicate
#   } ## end species
# } ## end scenario
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

foreach(file_use = files_to_run[21]) %dopar% {
  invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

  run_WHAM(yrs_use = 2010:2080, ## years to run the assessment
           file_suffix = basename(file_use))
} ## end files loop


## Summarize results across all species, scenarios, simulations ----

tabund <- list.files(here::here('outputs','wham_runs'),
                     pattern = 'true_biomass.csv',
                     recursive = TRUE,
                     full.names = TRUE) %>%
  lapply(., FUN = read.csv) %>% bind_rows() %>%
  group_by(year,scenario, species) %>%
  summarize(med=median(abund_mean)/1000,
            lwr50=quantile(abund_mean, .25)/1000,
            upr50=quantile(abund_mean, .75)/1000,
            lwr95=quantile(abund_mean, .0275)/1000,
            upr95=quantile(abund_mean, .975)/1000)

ggplot(tabund, aes(x = year, y = med, fill = scenario, color = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr50, ymax = upr50),
              alpha = 0.1, color = NA) +
  scale_fill_manual(values = scenPal, labels = scenLabs)+
  scale_color_manual(values = scenPal, labels = scenLabs)+
  facet_wrap(~species, scales = 'free_y', labeller = as_labeller(sppLabs)) +
  labs(x = 'Year', y = 'True Biomass (kmt)', color = '', fill = '')+
  theme(legend.position = 'bottom')

ggsave(last_plot(),
       file =here('figs','biomass_by_scenario_50ci.png'),
       width = 9, height = 6, unit = 'in', dpi = 400)
