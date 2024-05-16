## evOsmose-Master
## execute the full study in parallel
## M S Kapur Feb 2024 ++
## maia.kapur@noaa.gov

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets


## Create input data ----
## Only need to do this if the simulations themselves or sampling protocol have changed
## (e.g., obs error, survey array, years or frequency thereof, etc)

#* parallel ----
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)
## one species, one replicate, one scenario, two fc scenarios = 2 mins
## three spp, one replicate, four scenarios  = 4 mins
foreach(scenario=1:4) %:%
  foreach(species = c(1)) %:%
  # foreach(species = c(sppLabs2$Var3[sppLabs2$Var4]+1)) %:%
  foreach(repl=1:10)  %dopar%  {

    invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets

    scen_use = scenario;
    sp_use = species;
    replicate_use = repl;

    for(fc_use in c(1,0.15,0.15001)){
      build_Data(scenario=scen_use,
                 sppIdx = sp_use,
                 repID = replicate_use,
                 yrs_use = 2010:2080, ## years to extract data for
                 # srv_selex = ifelse(fc_use == 0.15001, 8, NA), ## age at 50% selex
                 # obs_error = ifelse(fc_use == 0.15001, 0.1, NA), ## observation error for surveys
                 fractional_coverage_use = fc_use, ## fractional coverage of survey
                 srv_selex = NA, ## age at 50% selex
                 obs_error = NA, ## observation error for surveys
                 units = 'biomass',
                 units_scalar = 1,
                 do_GAM = FALSE)
    } ## end fractional coverage loop

  } ## end nested loop
# When you're done, clean up the cluster
stopImplicitCluster();stopCluster()

## Run WHAM model(s) ----
## list all the folders with outputs; can grep() or select from here
# greppy <- paste0('rep',c(0,"1\\b","2\\b",10,11,12), collapse = "|")

files_to_run <- list.dirs.depth.n( here::here('outputs','wham_runs'), n = 3) %>%
  .[grepl('2024-05-16/rep',.)] %>%
  .[grepl('Herring',.)] #%>%

# .[!grepl(greppy, .)]

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach(file_use = files_to_run[6:12]) %dopar% {
  invisible(lapply(list.files(here::here('R','functions'), full.names = TRUE), FUN=source)) ## load all functions and presets
  for(fc_use in c(1,0.15)){
    run_WHAM(yrs_use = 2010:2080, ## years to run the assessment
             fractional_coverage_use = fc_use, ## which survey setup to read from
             file_suffix = file_use)
  }
} ## end files loop
stopImplicitCluster();stopCluster()

## Summarize results across all species, scenarios, simulations ----

#* MRE by rep, scenario, species ----
files_to_run <- list.dirs.depth.n( here::here('outputs','wham_runs'), n = 4) %>%
  .[grepl('2024-05-16/rep',.)] %>%
  .[grepl('Herring',.)] %>%
  .[grepl('perfect_information',.)]


mre_all <- list.files(files_to_run,
                      pattern = 'ssb_mre.csv',
                      recursive = TRUE,
                      full.names = TRUE) %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  group_by(year,scenario, species) %>%
  summarize(med=median(MRE_scaled),
            lwr50=quantile(MRE_scaled, .25),
            upr50=quantile(MRE_scaled, .75),
            lwr95=quantile(MRE_scaled, .0275),
            upr95=quantile(MRE_scaled, .975))


ggplot(mre_all, aes(x = year, y = med,
                    color = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr95, ymax = upr95),
              alpha = 0.15, color = NA) +
  scale_y_continuous(limits = c(-100,100)) +
  scale_fill_manual(values = scenPal, labels = scenLabs)+
  scale_color_manual(values = scenPal, labels = scenLabs)+
  labs(x = 'Year', y = 'MRE SSB, %', color = '', fill = '')+
  geom_hline(yintercept = 0, color = 'pink') +
  theme(legend.position = 'none')+
  facet_wrap(~species,labeller = as_labeller(sppLabs))


ggsave(last_plot(),
       file =here('figs',paste0(Sys.Date(),'-MRE_by_scenario-95ci-perfectinfo.png')),
       width = 5, height = 5, unit = 'in', dpi = 400)


#* 2060 map of biomass by spp ----
tabund_spatial <- list.files(files_to_run,
                             pattern = 'true_biomass_ys\\b.csv',
                             recursive = TRUE,
                             full.names = TRUE) %>%
  .[grepl('rep0',.)] %>%
  lapply(., FUN = read.csv) %>%
  bind_rows()
# summarise(abundance_rescale_year =  rescale(tot_val, to=c(0,1), na.rm = T),
#           .by = c(year, species, scenario))

maplist = list()
for(spp in c(sppLabs2$Var2[sppLabs2$Var4])){
for(s in 1:4){

  p<- tabund_spatial %>%
    filter(year %in% c(2060) &
             scenario == scenLabs2$Var2[s] &
             species == spp)  %>%
    ggplot(data = ., aes(x = lat, y = long,
                         fill = abundance_rescale)) +
    theme_void()+
    geom_raster()+
    # scale_fill_viridis_c(na.value = NA)+
    scale_fill_gradient2(mid = "#FFF5D1",
                         low = "#efeee7",
                         high = scenPal[s])+
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
  ggsave(p,
         file =here('figs',paste0('biomass_maps_by_scenario_2060-',scenLabs2$Var2[s],"-",spp,'.png')),
         width = 2, height = 2, unit = 'in', dpi = 400)
  # invisible(mapply(ggsave,
  #                  file = here('figs',paste0('biomass_maps_by_scenario_2060-', scenLabs2$Var2[s],"-",spp,'.png')),
  #                  plot=l))

} ## end scenario

# png(here('figs','biomass_maps_by_scenario_2045.png'),
#     width = 3, height = 10, unit = 'in', res = 400)
# Rmisc::multiplot(plotlist = maplist, cols = 1)
# dev.off()
#
# ggsave(Rmisc::multiplot(plotlist = maplist, cols = 1),
#        file =here('figs',paste0('biomass_maps_by_scenario_2060-',species,'.png')),
#        width = 3, height = 10, unit = 'in', dpi = 400)

} ## end species



#* biomass timeseries by scenarios ----
for(species in c(sppLabs2$Var2[sppLabs2$Var4])){

  tabund <- list.files(files_to_run,
                       pattern ="*true_biomass_y\\b.csv",
                       recursive = TRUE,
                       full.names = TRUE) %>%
    .[grepl(species,.)] %>%
    lapply(., FUN = read.csv) %>%
    bind_rows() %>%
    group_by(year,scenario, species) %>%
    summarize(med=median(total_biomass  )/1000,
              lwr50=quantile(total_biomass  , .25)/1000,
              upr50=quantile(total_biomass  , .75)/1000,
              lwr95=quantile(total_biomass  , .0275)/1000,
              upr95=quantile(total_biomass  , .975)/1000)
    # summarize(med=median(ssb_true )/1000,
    #           lwr50=quantile(ssb_true , .25)/1000,
    #           upr50=quantile(ssb_true , .75)/1000,
    #           lwr95=quantile(ssb_true , .0275)/1000,
    #           upr95=quantile(ssb_true , .975)/1000)

  ggplot(tabund, aes(x = year, y = med, fill = scenario, color = scenario)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr50, ymax = upr50),
                alpha = 0.2, color = NA) +
    scale_fill_manual(values = scenPal, labels = scenLabs)+
    scale_color_manual(values = scenPal, labels = scenLabs)+
    {if(species == 'AtlanticHerring') scale_y_continuous(limits = c(1500,4000))} +
    {if(species == 'AtlanticCod') scale_y_continuous(limits = c(250,1000))} +
    {if(species == 'EuropeanSprat') scale_y_continuous(limits = c(2000,9000))} +
    facet_wrap(~species, scales = 'free_y', labeller = as_labeller(sppLabs)) +
    labs(x = 'Year', y = 'True SSB (kmt)', color = '', fill = '') +
    theme(legend.position = 'none')
  # theme_bw()

  ggsave(last_plot(),
         file =here('figs',paste0('biomass_by_scenario-95ci-',species,'.png')),
         width = 3, height = 4, unit = 'in', dpi = 400)
} ## end species

#* survey time series by scenario ----
for(species in c(sppLabs2$Var2[sppLabs2$Var4])){
  for(fc in c(1,0.15)){

    ## strip one species and coverage-specific replicate
    sobs <- list.files(files_to_run,
                       pattern = paste0(fc,'-survey_obs_biomass.csv'),
                       recursive = TRUE,
                       full.names = TRUE) %>%
      .[grepl(species,.)] %>%
      lapply(., FUN = read.csv) %>%
      bind_rows() %>%
      mutate(abund_mean = abund_mean/1000) %>%
      # filter(replicate == 0) %>%
      group_by(year,scenario, species) %>%
      summarize(med=median(abund_mean  ),
                lwr50=quantile(abund_mean  , .25),
                upr50=quantile(abund_mean  , .75),
                lwr95=quantile(abund_mean  , .0275),
                upr95=quantile(abund_mean  , .975))

    # sobs$scenario <- factor(sobs$scenario, levels = scenLabs2$Var2[c(2,1,3,4)])

    ggplot(sobs, aes(x = year,
                     y = med ,
                     group= interaction(replicate, scenario),
                     fill = scenario,
                     color = scenario)) +
      geom_point() +
      {if(species == 'AtlanticHerring') scale_y_continuous(limits = c(1500,4000))} +
      {if(species == 'AtlanticCod') scale_y_continuous(limits = c(250,1000))} +
      {if(species == 'EuropeanSprat') scale_y_continuous(limits = c(2000,9000))} +
      # geom_errorbar(aes(ymin = abund_mean-abund_mean*abund_cv,
      #                   ymax =  abund_mean+abund_mean *abund_cv,
      #                   color = scenario ),
      #               alpha = 0.2, width = 0) +
      {if(fc != 1) geom_errorbar(aes(ymin = lwr50,   ymax =  upr50, color = scenario ),  alpha = 0.2, width = 0)} +
      {if(fc==1) geom_errorbar(aes(ymin = med-0.05*med,   ymax =   med+0.05*med, color = scenario ),  alpha = 0.2, width = 0)} +
      scale_fill_manual(values = scenPal, labels = scenLabs)+
      scale_color_manual(values = scenPal, labels = scenLabs)+

      theme(  strip.background = element_blank()  ,
              legend.position = 'none'
      )+
      labs(x = 'Year', y = 'Survey Biomass (kmt)', color = '', fill = '') +
    facet_wrap( ~ species, scales = 'free_y', labeller = as_labeller(sppLabs)) #+
    # facet_wrap(~scenario, ncol = 1, scales = 'fixed')


    ggsave(last_plot(),
           file =here('figs',paste0('survobs_by_scenario_95ci-',species,'-',fc,'.png')),
           width = 3, height = 4, unit = 'in', dpi = 400)
  } ## end fractional coverage
} ## end species

# deprecated ----

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
