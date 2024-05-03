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
foreach(scenario=1) %:%
  foreach(species = c(sppLabs2$Var3[sppLabs2$Var4]+1))%:%
  foreach(replicate=1)  %dopar%  {
    invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

    scen_use = scenario; sp_use = species; replicate_use = replicate

    build_Data(scenario=scen_use,
               sppIdx = sp_use,
               repID = replicate_use,
               yrs_use = 2010:2080, ## years to extract data for
               srv_selex = 7, ## age at 50% selex
               obs_error = 0.2, ## observation error for surveys
               units = 'biomass',
               units_scalar = 1,
               do_GAM = FALSE)
  } ## end nested loop
# When you're done, clean up the cluster
stopImplicitCluster();stopCluster()

## Run WHAM model(s) ----
## list all the folders with outputs; can grep() or select from here
files_to_run <- list.dirs.depth.n( here::here('outputs','wham_runs'), n = 3)%>%
  .[grepl('2024-05-03/rep',.)]

foreach(file_use = files_to_run[9]) %dopar% {
  invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

  run_WHAM(yrs_use = 2010:2080, ## years to run the assessment
           file_suffix = files_to_run[2])
} ## end files loop


## Summarize results across all species, scenarios, simulations ----

#* 2060 map of biomass by spp ----
#*
tabund_spatial <- list.files(here::here('outputs','wham_runs'),
                             pattern = 'true_biomass_spatial.csv',
                             recursive = TRUE,
                             full.names = TRUE) %>%
  lapply(., FUN = read.csv) %>% bind_rows() #%>%
summarise(abundance_rescale_year =  rescale(tot_val, to=c(0,1), na.rm = T),
          .by = c(year, species, scenario))

maplist = list()
for(s in 1:4){
  maplist[[s]] <-  tabund_spatial %>%
    filter(year %in% c(2045) &
             scenario == scenLabs2$Var2[s],
           species == "AtlanticCod",
           replicate == 0)  %>%
    ggplot(data = ., aes(x = lat, y = long,
                         fill = abundance_rescale))+
    theme_void()+
    geom_raster()+
    # scale_fill_viridis_c(na.value = NA)+
    scale_fill_gradient2(low = "#efeee7", mid = alpha(scenPal[s],0.2), high = scenPal[s])+
    theme(
      legend.background = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position =  'none')
  # facet_grid( scenario~.,
  #             labeller = as_labeller(scenLabs))

}

png(here('figs','biomass_maps_by_scenario_2045.png'),
    width = 3, height = 10, unit = 'in', res = 400)
Rmisc::multiplot(plotlist = maplist, cols = 1)
dev.off()

ggsave(Rmisc::multiplot(plotlist = maplist, cols = 1),
       file =here('figs','biomass_maps_by_scenario_2045.png'),
       width = 3, height = 10, unit = 'in', dpi = 400)

#* biomass timeseries by scenarios ----
tabund <- list.files(here::here('outputs','wham_runs'),
                     pattern ="*true_biomass\\b.csv",
                     recursive = TRUE,
                     full.names = TRUE) %>%
  .[grepl('2024-04-30',.)] %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  group_by(year,scenario, species) %>%
  summarize(med=median(abund_mean)/1000,
            lwr50=quantile(abund_mean, .25)/1000,
            upr50=quantile(abund_mean, .75)/1000,
            lwr95=quantile(abund_mean, .0275)/1000,
            upr95=quantile(abund_mean, .975)/1000)



ggplot(tabund, aes(x = year, y = med, fill = scenario, color = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr50, ymax = upr50),
              alpha = 0.2, color = NA) +
  scale_fill_manual(values = scenPal, labels = scenLabs)+
  scale_color_manual(values = scenPal, labels = scenLabs)+
  facet_wrap(~species, scales = 'free_y', labeller = as_labeller(sppLabs)) +
  labs(x = 'Year', y = 'True Biomass (kmt)', color = '', fill = '') +
  theme(legend.position = 'top')
  # theme_bw()

ggsave(last_plot(),
       file =here('figs','biomass_by_scenario_50ci.png'),
       width = 8, height = 3, unit = 'in', dpi = 400)

## survey time series by scenario ---

sobs <- list.files(here::here('outputs','wham_runs'),
           pattern = 'survey_obs_biomass.csv',
           recursive = TRUE,
           full.names = TRUE) %>%
  lapply(., FUN = read.table, sep = ' ', header = TRUE) %>%
  bind_rows() %>%
  mutate(abund_mean = abund_mean/1000)


ggplot(sobs, aes(x = year, y = abund_mean , fill = scenario, color = scenario)) +
  geom_point() + geom_line()+
  geom_errorbar(aes(ymin = abund_mean-abund_mean*abund_cv,
                    ymax =  abund_mean+abund_mean *abund_cv,
                    color = scenario ),
              alpha = 0.2, width = 0) +
  scale_fill_manual(values = scenPal, labels = scenLabs)+
  scale_color_manual(values = scenPal, labels = scenLabs)+
  facet_wrap( ~ species, scales = 'free_y') +
  theme(  strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position ='none')+
  labs(x = 'Year', y = 'Survey Biomass (kmt)', color = '', fill = '')


ggsave(last_plot(),
       file =here('figs','survobs_by_scenario_95ci.png'),
       width = 8, height = 3, unit = 'in', dpi = 400)
