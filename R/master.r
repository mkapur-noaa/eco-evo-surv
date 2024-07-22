## evOsmose-Master
## execute the full study in parallel
## M S Kapur Feb 2024 ++
## maia.kapur@noaa.gov

rm(list = ls()) ## clear workspace
invisible(lapply(list.files(
  here::here('R', 'functions'), full.names = TRUE
), FUN = source)) ## load all functions and presets


## Create input data ----
## Only need to do this if the simulations themselves or sampling protocol have changed
## (e.g., obs error, survey array, years or frequency thereof, etc)

#* parallel ----
## one species, one replicate, one scenario, two fc scenarios = 2 mins
## three spp, one replicate, four scenarios  = 4 mins
## four spp, one replicate, four scenarios = 7 mins
## four spp, five replicates, one scenario = CRASH
## one spp, eight replicates, one scenario = 6 mins
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach(scenario = 1:4) %:%
  foreach(species = c(sppLabs2$Var3[sppLabs2$Var4] + 1)[2:4]) %:%
  # foreach(species = c(sppLabs2$Var3[sppLabs2$Var4] + 1)) %:%
  foreach(repl = c(1:27))  %dopar%  {
    invisible(lapply(list.files(
      here::here('R', 'functions'), full.names = TRUE
    ), FUN = source)) ## load all functions and presets

    scen_use = scenario
    sp_use = species
    replicate_use = repl

    for (fc_use in c(1)) {
      build_Data(
        scenario = scen_use,
        sppIdx = sp_use,
        repID = replicate_use,
        yrs_use = 2010:2080,
        ## years to extract data for
        fractional_coverage_use = fc_use,
        ## fractional coverage of survey
        srv_selex = ifelse(fc_use == 0.15001, 'mat', NA),
        ## how survey selex is specified
        obs_error = ifelse(fc_use == 0.15001, 0.1, NA),
        ## observation error for surveys
        # srv_selex = NA, ## age at 50% selex
        # obs_error = NA, ## observation error for surveys
        units = 'biomass',
        units_scalar = 1,
        date.use = "2024-07-09",
        do_GAM = FALSE
      )
    } ## end fractional coverage loop

  } ## end nested loop
# When you're done, clean up the cluster
stopImplicitCluster()
stopCluster()

#* check how many were made ----
dat_files <-
  list.files(here::here('outputs', 'wham_runs'),
             recursive = TRUE,
             pattern = '*wham_survey.csv', full.names = TRUE) %>%
  .[grepl('2024-07-09/rep', .)] %>%
  basename(.) %>%
  gsub("-wham_survey.csv","",.) %>%
  stringr::str_split_fixed(., "-", 4) %>%
  data.frame()

## species, scenario should be 56 each
dat_files %>%
  summarise(n=n(), .by = c(X1,X2)) %>%
  arrange(X1,n) %>%
  filter(n!=56)

## Run WHAM model(s) ----
## list all the folders with outputs; can grep() or select from here
# greppy <- paste0('rep',c(0,"1\\b","2\\b",3,4,5), collapse = "|")
# greppy <- paste0('rep',15:20, collapse = "|")

files_to_run <-
  list.dirs.depth.n(here::here('outputs', 'wham_runs'), n = 3) %>%
  .[grepl('2024-07-09/rep', .)]  %>%
  .[!grepl('Saithe',.)] #%>%
  # .[grepl('Herring',.)] #%>%
  # .[grepl(greppy, .)] %>%
  # .[!grepl('noCC_noEvo', .)] %>%
# .[grepl('\\bCC_Evo',.)] %>%
# .[!grepl(paste(files_to_run_done,collapse ="|"), .)]


# files_to_run_done <- c(files_to_run_done, files_to_run[1:20])
# file_suffix =file_use= files_to_run[2];yrs_use = 2010:2080;fractional_coverage_use=fc_use = 0.15

cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach(file_use = files_to_run[315:329]) %:%
  foreach(ewaa_use = c('perfect','averaged')[1])  %:%
  foreach(q_treatment_use = c('fixed','estimated')[1])  %:%
  foreach(fc_use = c(1, 0.15)) %dopar%  {
    invisible(lapply(list.files(
      here::here('R', 'functions'), full.names = TRUE
    ), FUN = source)) ## load all functions and presets

    run_WHAM(
      yrs_use = 2010:2080,## years to run the assessment
      fractional_coverage_use = fc_use, ## which survey setup to read from
      ewaa  = ewaa_use, ## which ewaa input to read from
      q_treatment = q_treatment_use, ## testing only; whether or not Q is estimated
      file_suffix = file_use
    )
  } ## end dopar loop
stopImplicitCluster()
stopCluster()

## Summarize results across all species, scenarios, simulations ----

#* MRE by rep, scenario, species ----
# files_to_run <- list.dirs.depth.n( here::here('outputs','wham_runs'), n = 4) %>%
#   .[grepl('2024-05-16/rep',.)] %>%
#   .[grepl('Herring',.)] %>%
#   .[grepl('perfect_information',.)]

mre_all0 <- list.files(
  files_to_run,
  pattern = 'mre.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  filter(q_treatment == 'fixed' & ewaa == 'perfect' &fc != 0.15001) %>%
  distinct() ## drop old runs and duplicates

## check how many ran (should be 28 each)
summarise(mre_all0,
          nrep = n() / 71,
          .by = c(scenario, species, fc))

## IDs of fails
summarise(mre_all0,
          nrep = n() / 71,
          .by = c(scenario, replicate, species)) %>%
  tidyr::complete(replicate, scenario, species) %>%
  filter(is.na(nrep)) %>%
  arrange(species, scenario, replicate) %>%
  View()

mre_all <- mre_all0 %>%
  group_by(year, scenario, species, fc, ewaa) %>%
  summarize(
    med2 = median(MRE_ssb*100),
    med  = median(MRE_scaled),
    lwr50 = quantile(MRE_scaled, .25),
    upr50 = quantile(MRE_scaled, .75),
    lwr95 = quantile(MRE_scaled, .0275),
    upr95 = quantile(MRE_scaled, .975)
  )

mre_all$ewaa <- factor(mre_all$ewaa, levels = c('perfect','averaged'))
mre_all$fc <- factor(mre_all$fc, levels = c(1, 0.15))

for (spp in unique(mre_all$species)) {
  # for(fcc in unique(mre_all$fc)){
  # ggplot(mre_all,
  ggplot(subset(mre_all,  species == spp),
         aes(
           x = year,
           y = med,
           color = scenario,
           fill = scenario
         )) +
    geom_hline(yintercept = 0,
               color = 'grey50',
               linetype = 'dotted') +
    geom_line() +
    geom_ribbon(aes(ymin = lwr95, ymax = upr95),
                alpha = 0.15,
                color = NA) +
    scale_y_continuous(limits = c(-60, 60), breaks = seq(-50, 50, 10)) +
    scale_fill_manual(values = scenPal, labels = scenLabs) +
    scale_color_manual(values = scenPal, labels = scenLabs) +
    labs(
      x = 'Year',
      y = 'MRE SSB, %',
      color = '',
      fill = ''
    ) +
    # facet_wrap(   ~ fc, labeller = as_labeller(fcLabs), ncol = 1) +

    facet_grid( fc ~ ewaa) +
    theme(
      # legend.position = 'none',
      # strip.text.x = element_blank(),
      axis.title = element_blank()
    )


  # ggsave(  Rmisc::multiplot(plotlist = list(plist[[1]],plist[[2]]), cols = 2),
  #        file =here('figs',paste0(Sys.Date(),'-',fcc,'-MRE_by_scenario-95ci.png')),
  #        width = 3, height = 4, unit = 'in', dpi = 400)

  # png(here('figs',paste0(Sys.Date(),'-',fcc,'-MRE_by_scenario-95ci.png')),
  #     width = 5, height = 3, unit = 'in', res = 400)
  # Rmisc::multiplot(plotlist = list(plist[[1]],plist[[2]]), cols = 2)
  # dev.off()

  # ggsave(
  #   last_plot(),
  #   file = here(
  #     'figs',
  #     paste0(Sys.Date(), '-', spp, '-MRE_by_scenario-95ci.png')
  #   ),
  #   width = 4,
  #   height = 6,
  #   unit = 'in',
  #   dpi = 400
  # )

}

# }
write.csv(
  mre_all0,
  file = here('outputs', 'summary_data', paste0(Sys.Date(), '-mre_all.csv')),
  row.names = FALSE
)

## OM Outputs ----
#* 2060 map of biomass by spp ----
tabund_spatial <- list.files(
  files_to_run,
  pattern = 'true_biomass_ys\\b.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  .[grepl('rep0', .)] %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  mutate(
    abundance_rescale_year =  rescale(tot_val, to = c(0, 1), na.rm = T),
    .by = c(scenario, species)
  )

# tabund_spatial$scenario <- factor(tabund_spatial$scenario,
#                                   levels = scenLabs2$Var2[])
for (spp in c(sppLabs2$Var2[sppLabs2$Var4])) {
  maplist = list()
  for (s in 1:4) {
    maplist[[s]] <- tabund_spatial %>%
      filter(year %in% c(2020, 2040, 2060) &
               scenario == scenLabs2$Var2[s] &
               species == spp)  %>%
      ggplot(data = ., aes(x = lat, y = long,
                           fill = abundance_rescale_year)) +
      theme_void() +
      geom_raster() +
      # scale_fill_viridis_c(na.value = NA)+
      scale_fill_gradient2(
        low = "#FFF5D1",
        mid = "#efeee7",
        high = scenLabs2$Pal2[s],
        # high = 'black'
      ) +
      # high = scenLabs2$Pal[s])+
      theme(
        strip.text.y = element_text(angle = 90),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position =  'none'
      ) +
      facet_grid(scenario ~ year)
    png(
      here(
        'figs',
        paste0(
          Sys.Date(),
          '-biomass_maps_by_scenario_timestep-',
          spp,
          '.png'
        )
      ),
      width = 4,
      height = 6,
      unit = 'in',
      res = 400
    )
    Rmisc::multiplot(plotlist = maplist, cols = 1)
    dev.off()
    #
    # ggsave(p,
    #        file =here('figs',paste0('biomass_maps_by_scenario_timestep-',
    #                                 scenLabs2$Var2[s],"-",
    #                                 spp,'.png')),
    #        width =3, height = 6/4, unit = 'in', dpi = 400)
    # invisible(mapply(ggsave,
    #                  file = here('figs',paste0('biomass_maps_by_scenario_2060-', scenLabs2$Var2[s],"-",spp,'.png')),
    #                  plot=l))
  }
} ## end scenario



# Rmisc::multiplot(plotlist = maplist, cols = 1)
# ggsave(,
#        file =here('figs',paste0('biomass_maps_by_scenario_2060-',species,'.png')),
#        width = 3, height = 10, unit = 'in', dpi = 400)

# } ## end species



#* biomass timeseries by scenarios ----
for (species in c(sppLabs2$Var2[sppLabs2$Var4])) {
  tabund <- list.files(
    files_to_run,
    pattern = "*true_biomass_y\\b.csv",
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    .[grepl(species, .)] %>%
    lapply(., FUN = read.csv) %>%
    bind_rows() %>%
    group_by(year, scenario, species) %>%
    summarize(
      med = median(total_biomass) / 1000,
      lwr50 = quantile(total_biomass  , .25) / 1000,
      upr50 = quantile(total_biomass  , .75) / 1000,
      lwr95 = quantile(total_biomass  , .0275) / 1000,
      upr95 = quantile(total_biomass  , .975) / 1000
    )
  # summarize(med=median(ssb_true )/1000,
  #           lwr50=quantile(ssb_true , .25)/1000,
  #           upr50=quantile(ssb_true , .75)/1000,
  #           lwr95=quantile(ssb_true , .0275)/1000,
  #           upr95=quantile(ssb_true , .975)/1000)

  ggplot(tabund, aes(
    x = year,
    y = med,
    fill = scenario,
    color = scenario
  )) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr50, ymax = upr50),
                alpha = 0.2,
                color = NA) +
    scale_fill_manual(values = scenPal, labels = scenLabs) +
    scale_color_manual(values = scenPal, labels = scenLabs) +
    # {if(species == 'AtlanticHerring') scale_y_continuous(limits = c(1500,4000))} +
    # {if(species == 'AtlanticCod') scale_y_continuous(limits = c(250,1000))} +
    # {if(species == 'EuropeanSprat') scale_y_continuous(limits = c(2000,9000))} +
    facet_wrap( ~ species, scales = 'free_y', labeller = as_labeller(sppLabs)) +
    labs(
      x = 'Year',
      y = 'True SSB (kmt)',
      color = '',
      fill = ''
    ) +
    theme(legend.position = 'none')
  # theme_bw()

  ggsave(
    last_plot(),
    file = here(
      'figs',
      paste0('biomass_by_scenario-95ci-', species, '.png')
    ),
    width = 3,
    height = 3,
    unit = 'in',
    dpi = 400
  )
} ## end species

#* survey time series by scenario ----

## summarize cv trends
for (species in c(sppLabs2$Var2[sppLabs2$Var4])) {
  for (fc in c(1, 0.15, 0.15001)) {
    ## strip one species and coverage-specific replicate
    sobs <- list.files(
      files_to_run,
      pattern = paste0(fc, '-survey_obs_biomass.csv'),
      recursive = TRUE,
      full.names = TRUE
    ) %>%
      .[grepl(species, .)] %>%
      lapply(., FUN = read.csv) %>%
      bind_rows() %>%
      mutate(abund_mean = abund_mean / 1000) %>%
      # filter(replicate == 0) %>%
      group_by(year, scenario, species) %>%
      summarize(
        med = median(abund_mean),
        cv_median = median(abund_cv),
        lwr50 = quantile(abund_mean  , .25),
        upr50 = quantile(abund_mean  , .75),
        lwr95 = quantile(abund_mean  , .0275),
        upr95 = quantile(abund_mean  , .975)
      )

    # sobs$scenario <- factor(sobs$scenario, levels = scenLabs2$Var2[c(2,1,3,4)])

    ggplot(sobs,
           aes(
             x = year,
             y = med ,
             group = interaction(scenario),
             fill = scenario,
             color = scenario
           )) +
      geom_point() +
      # {if(species == 'AtlanticHerring') scale_y_continuous(limits = c(1000,4000))} +
      # {if(species == 'AtlanticCod') scale_y_continuous(limits = c(0,1000))} +
      # {if(species == 'EuropeanSprat') scale_y_continuous(limits = c(1500,8000))} +
      # geom_errorbar(aes(ymin = abund_mean-abund_mean*abund_cv,
      #                   ymax =  abund_mean+abund_mean *abund_cv,
      #                   color = scenario ),
      #               alpha = 0.2, width = 0) +
      {
        if (fc != 1)
          geom_errorbar(
            aes(
              ymin = lwr50,
              ymax =  upr50,
              color = scenario
            ),
            alpha = 0.2,
            width = 0
          )
      } +
      {
        if (fc == 1)
          geom_errorbar(
            aes(
              ymin = med - 0.05 * med,
              ymax =   med + 0.05 * med,
              color = scenario
            ),
            alpha = 0.2,
            width = 0
          )
      } +
      scale_fill_manual(values = scenPal, labels = scenLabs) +
      scale_color_manual(values = scenPal, labels = scenLabs) +

      theme(strip.background = element_blank()  ,
            legend.position = 'none') +
      labs(
        x = 'Year',
        y = 'Survey Biomass (kmt)',
        color = '',
        fill = ''
      ) +
      facet_wrap(~ species,
                 scales = 'free_y',
                 labeller = as_labeller(sppLabs)) #+
    # facet_wrap(~scenario, ncol = 1, scales = 'fixed')


    ggsave(
      last_plot(),
      file = here(
        'figs',
        paste0('survobs_by_scenario_95ci-', species, '-', fc, '.png')
      ),
      width = 3,
      height = 3,
      unit = 'in',
      dpi = 400
    )
  } ## end fractional coverage
} ## end species

#* input WAA at maturity across scenario----
strip_waa <- function(x) {
  tmpwaa <- read.table(x, sep = ' ', header = TRUE)
  names(tmpwaa) <- 1:ncol(tmpwaa)
  tmpwaa$year <- 2010:(2009 + nrow(tmpwaa))

  filen2 = basename(x)
  # filen2 <- basename(dirname(dirname(file_suffix)))
  tmpwaa$scenario <- strsplit(filen2, '-')[[1]][2]
  tmpwaa$rep <- strsplit(filen2, '-')[[1]][3]
  tmpwaa$species <-   strsplit(filen2, '-')[[1]][1]
  waa <-
    reshape2::melt(tmpwaa, id = c('year', 'rep', 'scenario', 'species'))
  names(waa)[5:6] <- c('age', 'weight_kg')
  waa$age <- as.numeric(waa$age)
  return(waa)
}

waa_index <- list.files(
  files_to_run,
  pattern = 'wham_waa_ssb_perfect.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  # .[grepl('rep0',.)] %>%
  lapply(., FUN = strip_waa) %>%
  bind_rows() %>%
  filter(!is.na(weight_kg)) %>%
  group_by(year, scenario, species, age) %>%
  summarize(
    med = median(weight_kg),
    lwr50 = quantile(weight_kg, .25),
    upr50 = quantile(weight_kg, .75),
    lwr95 = quantile(weight_kg, .0275),
    upr95 = quantile(weight_kg, .975)
  )

for (spp in c(sppLabs2$Var2[sppLabs2$Var4])) {
  ggplot(
    subset(waa_index, age == 4 &
             year < 2081 & species == spp),
    aes(
      x = year,
      y = med,
      group = interaction(scenario),
      fill = scenario,
      color = scenario
    )
  ) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr50, ymax = upr50),
                alpha = 0.2,
                color = NA) +
    scale_fill_manual(values = scenPal, labels = scenLabs) +
    scale_color_manual(values = scenPal, labels = scenLabs) +
    theme(strip.background = element_blank()  ,
          legend.position = 'none') +
    labs(
      x = 'Year',
      y = 'EWAA @ ~50% maturity, kg',
      color = '',
      fill = ''
    )
  # ggsave(
  #   last_plot(),
  #   file = here('figs', paste0(Sys.Date(),'-ewaa4_by_scenario_95ci-', spp, '.png')),
  #   width = 3,
  #   height = 4,
  #   unit = 'in',
  #   dpi = 400
  # )
}

#* input agecomps by scenario ----
acomps_index <- list.files(
  files_to_run,
  pattern = 'survey_obs_agecomp',
  recursive = TRUE,
  full.names = TRUE
) %>%
  # .[grepl('rep0',.)] %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  filter(fc %in% c(1, 0.15)) %>%
  mutate(scount = sum(count),
         .by = c(year, replicate, fc, scenario)) %>%
  mutate(frequency = count / scount)  %>%
  group_by(year, fc, scenario, species, age) %>%
  summarize(
    med = median(frequency),
    lwr50 = quantile(frequency, .25),
    upr50 = quantile(frequency, .75),
    lwr95 = quantile(frequency, .0275),
    upr95 = quantile(frequency, .975)
  )

ggplot(
  subset(acomps_index, year == 2040 & fc == 0.15),
  aes(
    x = age,
    group = interaction(fc, year, scenario),
    color = scenario,
    fill = scenario
  )
) +
  geom_line(aes(y = med)) +
  geom_ribbon(aes(ymin = lwr50, ymax = upr50),
              alpha = 0.2,
              color = NA) +
  scale_x_continuous(breaks = seq(0, 17, by = 2), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = scenPal) +
  scale_color_manual(values = scenPal) +
  labs(x = 'Age',
       y = 'Frequency',
       color = '',
       fill = '') +
  theme(legend.position = 'none', axis.text.y = element_blank()) +
  facet_wrap(. ~ fc, labeller = as_labeller(fcLabs))

ggsave(
  last_plot(),
  file = here(
    'figs',
    paste0('acompssimplesurvey-', 'atlanticherring', '.png')
  ),
  width = 3,
  height = 4,
  unit = 'in',
  dpi = 400
)


# deprecated ----
# mre_all %>%
#   group_by(scenario,fc) %>%
#   summarise(mam = median(abs(med))) %>%
#   arrange(mam)
# ggplot(mre_all0, aes(x = year, y = MRE_scaled,
#                      group = interaction(replicate, scenario),
#                      color = scenario, fill = scenario)) +
#   geom_hline(yintercept = 0, color = 'grey50', linetype = 'dotted') +
#   geom_line() +
#   # geom_ribbon(aes(ymin = lwr50, ymax = upr50),
#   #             alpha = 0.15, color = NA) +
#   scale_y_continuous(limits = c(-35,35)) +
#   scale_fill_manual(values = scenPal, labels = scenLabs)+
#   scale_color_manual(values = scenPal, labels = scenLabs)+
#   labs(x = 'Year', y = 'MRE SSB, %', color = '', fill = '')+
#   theme(legend.position = 'none')+
#   facet_grid(~fc)
#
#
# mre_all1 <- mre_all0 %>%
#   dplyr::select(year, replicate, scenario, fc, ssb_true, ssb_est) %>%
#   reshape2::melt(id = c('year','replicate','scenario','fc'))
#
# ggplot(mre_all1, aes(x = year, y = value,
#                      group = interaction(replicate,scenario, variable),
#                      linetype = variable,
#                      color = variable)) +
#   geom_line() +
#   # scale_fill_manual(values = scenPal, labels = scenLabs)+
#   # scale_color_manual(values = scenPal, labels = scenLabs)+
#   labs(x = 'Year', y = 'MRE SSB, %', color = '', fill = '')+
#   theme(legend.position = 'none')+
#   facet_grid(scenario~fc)
#* brute ----
#* stuff for megsie----
 # mre_all_mod$scenario <- factor(mre_all_mod$scenario)

# mre_all_megsie0 <- mre_all0 %>%
#   dplyr::select(species, scenario, fc, replicate,
#                 total_biomass,
#                 total_biomass_cv, year, MRE, MRE_scaled) %>%
#   mutate(fc = as.factor(fc)) %>%
#   dplyr::filter(fc %in% c("1","0.15"))
#
#
# mre_all_megsie1 <- merge(mre_all_megsie0,sppLabs2[,c('Var2','Var5')],
#                      by.x = 'species',
#                      by.y = 'Var2')
# names(mre_all_megsie1)[names(mre_all_megsie1)=='Var5'] <- 'trophic'

# sobs_cv <- list.files(files_to_run,
#                       pattern = paste0('-survey_obs_biomass.csv'),
#                       recursive = TRUE,
#                       full.names = TRUE) %>%
#   lapply(., FUN = read.csv) %>%
#   bind_rows() %>%
#   mutate(fc = as.factor(fc)) %>%
#   dplyr::filter(fc %in% c("1","0.15")) %>%
#   mutate(abund_mean = abund_mean/1000) %>%
#   dplyr::select(fc, year, replicate, species, scenario, abund_cv)
#
# ## merge dfs, collapses replicate
# mre_all_megsie <- merge(mre_all_megsie1,
#                         sobs_cv,
#                         by = c('year','replicate','fc','scenario','species'),
#                         all.x=TRUE)
#
# View(subset(mre_all_megsie, is.na(abund_cv)))
#
# write.csv(mre_all_megsie, file = here::here('outputs','summary_data','mre_for_megsie-2024-05-24.csv'),row.names = FALSE)
# mre_all_megsie$trophic <- mre_all_megsie$Var5
# mre_all_mod$species <- factor(mre_all_mod$species, levels = sppLabs2$Var5[sppLabs2$Var4])
# mre_all_mod$replicate <- factor(mre_all_mod$replicate)
# mre_all_mod$climate_change <- factor(grepl('\\bCC',mre_all_mod$scenario))
# mre_all_mod$evolution <- factor(grepl('_Evo',mre_all_mod$scenario))
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
