## Figures for manuscript & supplementary material
## assumes you have mre_all with final results (including survey data & diagnostics)

## build plotting summaries ----
bobs <- full_df %>%
  filter(fc == 1) %>% ## OM biomass is identical regardles of coverage
  mutate(ssb_true = ssb_true/1000) %>%
  group_by(year, scenario, species, fc) %>%
  dplyr::summarize(
    med = median(ssb_true),
    lwr80 = quantile(ssb_true, .10),
    upr80 = quantile(ssb_true, .90),
    lwr95 = quantile(ssb_true, .0275),
    upr95 = quantile(ssb_true, .975)
  )

sobs <- full_df %>%
  filter(!is.na(Obs)) %>%
  group_by(year, scenario, species, fc) %>%
  mutate(Obs = Obs/1000) %>%
  dplyr::summarize(
    med = median(Obs),
    lwr80 = med-med*quantile(abund_cv, .10),
    upr80 = med+med*quantile(abund_cv, .90),
    lwr95 = quantile(Obs, .0275),
    upr95 = quantile(Obs, .975)
  )

mre_summary <- full_df %>%
  group_by(year, scenario, species, fc,ewaa, q_treatment) %>%
  mutate(MRE_scaled = MRE_totbio*100) %>%
  summarize(
    med = median(MRE_scaled),
    lwr80 = quantile(MRE_scaled, 0.1),
    upr80 = quantile(MRE_scaled, 0.9),
    lwr95 = quantile(MRE_scaled, .0275),
    upr95 = quantile(MRE_scaled, .975)
  )
mre_summary$fc <- factor(mre_summary$fc,
                         levels = c(1.01,1,0.5,0.15,0.05))



## OM Maps ----
files_to_plot <- files_to_run[grepl('rep0',files_to_run)]
tabund_spatial <- list.files(
  files_to_plot,
  pattern = 'true_biomass_ys\\b.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  mutate(
    abundance_rescale_year =  rescale(tot_val, to = c(0, 1), na.rm = T),
    .by = c(scenario, species)
  )


tibble(tabund_spatial) %>%
  filter(year %in% c( 2040))  %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(x = lat, y = long,
                   fill = abundance_rescale_year)) +
      theme_void() +
      geom_raster() +
      theme(strip.text.y = element_text(angle = 270, size = 10))+
      facet_grid(species ~ scenario ,
                 labeller = labeller(species = sppLabs,
                                     scenario =  scenLabs))

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  theme(legend.position = 'bottom') &
  labs(fill = 'relative abundance in year 2040') &
  scale_fill_gradientn(colours = c("#95cec7","#2a9d8f","#e9c46a","#f4a261","#e76f51"))

ggsave(last_plot(),
       file = here::here('figs','manuscript','OM_Maps.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')

##  OM SB & Survey Obs ----
dobs <- bobs %>%
  mutate(variable = 'True OM Biomass') %>%
  rbind(., sobs %>% mutate(variable = 'Observed Survey Biomass')) %>%
  arrange(species, variable, fc)

dobs$variable <- factor(dobs$variable,levels = c('True OM Biomass',
                                                 'Observed Survey Biomass'))
dobs$fc <- factor(dobs$fc,levels = c(1.01,1,0.15,.5,0.05))
dobs$fc[dobs$variable=='True OM Biomass']<- 1.01

tibble(dobs) %>%
  filter(fc %in% c(1.01,1,0.15) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      theme(strip.text.y = element_text(angle = 270, size = 10))+

      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      labs(x = 'Year',y = '1000 mt', color = '', fill = '') +
      facet_grid(species ~ variable + fc,
                 labeller = labeller(species = sppLabs, fc = fcLabs))
  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  guides(fill=guide_legend(nrow=1,byrow=TRUE),
         color =guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') &
  scale_fill_manual(values = scenPal, labels = scenLabs) &
  scale_color_manual(values = scenPal, labels = scenLabs)

ggsave(last_plot(),
       file = here::here('figs','manuscript','OM_bio_survey.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')



## MRE ----
tibble(mre_summary) %>%
  filter(fc %in% c(1.01,1,0.15) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey80')+

      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      theme(strip.text.y = element_text(angle = 270, size = 10))+

      labs(x = 'Year',y = 'Relative Error', color = '', fill = '') +
      facet_grid(species ~ fc,
                 labeller = labeller(species = sppLabs, fc = fcLabs))
  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  guides(fill=guide_legend(nrow=1,byrow=TRUE),
         color =guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') &
  scale_fill_manual(values = scenPal, labels = scenLabs) &
  scale_color_manual(values = scenPal, labels = scenLabs)


ggsave(last_plot(),
       file = here::here('figs','manuscript','MRE_main.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')

## diagnostic pass rates ----
diags %>%
  filter(!is.na(runs_test_passed) ) %>%
  summarise(prop_passed = sum(runs_test_passed)/n(), .by = c('species','scenario','fc')) %>%
  ggplot(aes(y = prop_passed, x = factor(fc), fill = scenario)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = scenPal, labels = scenLabs) +
  theme(legend.position = 'none', strip.text.x = element_text(size = 5))+
  labs(x = 'Survey Coverage (%)', y = 'Runs test pass rate (%)')+
  facet_grid(species ~ scenario,
             labeller = labeller(species = sppLabs, scenario = scenLabs))

ggsave(last_plot(),
       file = here::here('figs','manuscript','diagnostics.png'),
       width = 8, height = 4, dpi = 500, unit = 'in')


## Supplementary Figures ----
#* diagnostic runs test vlaues
diags <- read.csv(here::here('outputs','summary_data',paste0("2025-01-14-diags.csv")))

diags %>%
  filter(runs_test_val < 0.07 & fc == 0.15 ) %>%
  # summarise(prop_passed = sum(runs_test_passed)/n(), .by = c('species','scenario','fc')) %>%
  ggplot(aes(x = runs_test_val, fill = scenario)) +
  # geom_bar(stat = 'identity', position = 'dodge') +
  geom_histogram(position = 'dodge') +
  scale_fill_manual(values = scenPal, labels = scenLabs) +
  theme(legend.position = 'none', strip.text.x = element_text(size = 5))+
  labs(x = 'Survey Coverage (%)', y = 'Runs test value, pass if < 0.05')+
  facet_grid(species ~ fc,
             labeller = labeller(species = sppLabs, scenario = scenLabs))

ggsave(last_plot(),
       file = here::here('figs','manuscript','diagnostics.png'),
       width = 8, height = 4, dpi = 500, unit = 'in')

#* OM Compositions ----
acomps_index <- list.files(
  files_to_run,
  pattern = 'survey_obs_agecomp',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = read.csv) %>%
  bind_rows() %>%
  filter(fc %in% c(1, 0.15)) %>%
  mutate(scount = sum(count),
         .by = c(year, replicate, fc, scenario)) %>%
  mutate(frequency = count / scount)  %>%
  group_by(fc, scenario, species, age) %>%
  # group_by(year, fc, scenario, species, age) %>%
  summarize(
    med = median(frequency),
    lwr80 = quantile(frequency, .1),
    upr80 = quantile(frequency, 0.9),
    lwr95 = quantile(frequency, .0275),
    upr95 = quantile(frequency, .975)
  ) %>%
  filter(fc == 0.15) %>%
  ungroup() %>%
  dplyr::select(-fc) %>%
  mutate(variable = 'survey')

acomps_catch <- list.files(
  files_to_run,
  pattern = 'wham_catch_at_age.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = function(x){
    tmp <- NULL
    tmp <- read.csv(x, sep = ' ')
    tmp$year = 2010:2099
    tmp$replicate <- as.numeric(gsub('rep','',basename(dirname(x))))
    tmp$scenario  <- strsplit(basename(dirname(dirname(x))),'-')[[1]][2]
    tmp$species <- strsplit(basename(dirname(dirname(x))),'-')[[1]][1]
    tmp <- tmp %>%
      dplyr::select(-total) %>%
      reshape2::melt(., id = c('year','replicate','scenario','species')) %>%
      mutate(age = as.numeric(gsub('X','',variable))) %>%
      dplyr::select(year, replicate, scenario, species, age, count = value)
    return(tmp)
  }) %>%
  bind_rows() %>%
  mutate(scount = sum(count),
         .by = c(year, replicate, scenario)) %>%
  mutate(frequency = count / scount)  %>%
  group_by(scenario, species, age) %>%
  # group_by(year,scenario, species, age) %>%
  summarize(
    med = median(frequency),
    lwr80 = quantile(frequency, .1),
    upr80 = quantile(frequency, 0.9),
    lwr95 = quantile(frequency, .0275),
    upr95 = quantile(frequency, .975)
  )%>%
  mutate(variable ='fishery')


tibble(bind_rows(acomps_index,acomps_catch)) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(
      x = age,
      group = interaction(scenario,species),
      color = scenario,
      fill = scenario
    )
    ) +
      geom_line(aes(y = med)) +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80),
                  alpha = 0.2,
                  color = NA) +
      # scale_x_continuous(breaks = seq(0, 17, by = 2), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = scenPal, labels = scenLabs) +
      scale_color_manual(values = scenPal, labels = scenLabs) +
      labs(x = 'Age',
           y = 'Frequency',
           color = '',
           fill = '') +
      theme(legend.position = 'none') +
      theme(strip.text.x = element_text(size =10))+
      facet_grid(species ~ variable,
                 labeller = labeller(species = sppLabs),
                 scales = "free_y")

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  theme(legend.position = 'bottom')&
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color =guide_legend(nrow=2,byrow=TRUE)) &
  theme(legend.position = 'bottom')

ggsave(last_plot(),
       file = here::here('figs','manuscript','OM_acomps.png'),
       width = 6, height = 8, dpi = 500, unit = 'in')

#* OM Catches ----
catch_index <- list.files(
  files_to_run,
  pattern = 'wham_catch_at_age.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = function(x){
    tmp <- NULL
    tmp <- read.csv(x, sep = ' ')
    tmp$year = 2010:2099
    tmp$replicate <- as.numeric(gsub('rep','',basename(dirname(x))))
    tmp$scenario  <- strsplit(basename(dirname(dirname(x))),'-')[[1]][2]
    tmp$species <- strsplit(basename(dirname(dirname(x))),'-')[[1]][1]

    return(tmp)
    }) %>%
  # mutate(year = 2010:2090) %>%
  bind_rows() %>%
  group_by(year, species, scenario) %>%
  mutate(total = total/1000) %>%
  summarize(
    med = median(total),
    lwr80 = quantile(total, 0.1),
    upr80 = quantile(total, 0.9),
    lwr95 = quantile(total, .0275),
    upr95 = quantile(total, .975)
  )

tibble(catch_index) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(
      x = year,
      # group = interaction(fc, year, scenario),
      color = scenario,
      fill = scenario
    )
    ) +
      geom_line(aes(y = med)) +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80),
                  alpha = 0.2,
                  color = NA) +
      # scale_x_continuous(breaks = seq(0, 17, by = 2), expand = c(0, 0)) +
      # scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = scenPal, labels = scenLabs) +
      scale_color_manual(values = scenPal, labels = scenLabs) +
      labs(x = 'Age',
           y = 'Catches (kmt)',
           color = '',
           fill = '') +
      theme(legend.position = 'none') +
      theme(strip.text.x = element_text(size =10))+
      facet_grid(species ~ .,
                 labeller = labeller(species = sppLabs),
                 scales = "free_y")

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  theme(legend.position = 'bottom')&
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color =guide_legend(nrow=2,byrow=TRUE))

ggsave(last_plot(),
       file = here::here('figs','manuscript','catches.png'),
       width = 6, height = 8, dpi = 500, unit = 'in')


#* OM maturity and mortality ----

maa_index <- list.files(
  files_to_run,
  pattern = 'wham_mortality.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = strip_waa) %>%
  bind_rows() %>%
  filter(!is.na(weight_kg)) %>%
  group_by(year, scenario, species, age) %>%
  summarize(
    med = median(weight_kg),
    lwr80 = quantile(weight_kg, .1),
    upr80 = quantile(weight_kg, .90),
    lwr95 = quantile(weight_kg, .0275),
    upr95 = quantile(weight_kg, .975)
  )%>%
  mutate(variable = 'mortality')

mataa_index <- list.files(
  files_to_run,
  pattern = 'wham_maturity.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = strip_waa) %>%
  bind_rows() %>%
  filter(!is.na(weight_kg)) %>%
  group_by(year, scenario, species, age) %>%
  summarize(
    med = median(weight_kg),
    lwr80 = quantile(weight_kg, .1),
    upr80 = quantile(weight_kg, .90),
    lwr95 = quantile(weight_kg, .0275),
    upr95 = quantile(weight_kg, .975)
  )%>%
  mutate(variable = 'maturity')

tibble(bind_rows(maa_index,mataa_index)) %>%
  filter(age == 3 & year < 2081) %>%
  # mutate(variable = ifelse(variable == 'maturity','age 3 maturity',
  #                          'age 3 mortality')) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(
      x = year,
      y = med,
      group = interaction(scenario),
      fill = scenario,
      color = scenario
    )
    ) +
      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80),
                  alpha = 0.2,
                  color = NA) +
      scale_fill_manual(values = scenPal, labels = scenLabs) +
      scale_color_manual(values = scenPal, labels = scenLabs) +
      theme(strip.background = element_blank(),
            legend.position = 'bottom') +

      guides(fill=guide_legend(nrow=1,byrow=TRUE),
             color =guide_legend(nrow=1,byrow=TRUE))+
      facet_grid( species+variable ~ . ,
                 scales = 'free',
                 labeller = labeller(species = sppLabs,
                                     scenario =  scenLabs))+
      labs(
        x = 'Year',
        y = '',
        color = '',
        fill = ''
      )

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color =guide_legend(nrow=2,byrow=TRUE)) &
  theme(legend.position = 'bottom')


ggsave(last_plot(),
       file = here::here('figs','manuscript','OM_mortality_maturity.png'),
       width = 6, height = 10, dpi = 500, unit = 'in')

#* OM mortality ----


#* OM  WAA  ----
waa_index <- list.files(
  files_to_run,
  pattern = 'wham_waa_ssb_perfect.csv',
  recursive = TRUE,
  full.names = TRUE
) %>%
  lapply(., FUN = strip_waa) %>%
  bind_rows() %>%
  filter(!is.na(weight_kg)) %>%
  group_by(year, scenario, species, age) %>%
  summarize(
    med = median(weight_kg),
    lwr80 = quantile(weight_kg, .1),
    upr80 = quantile(weight_kg, .90),
    lwr95 = quantile(weight_kg, .0275),
    upr95 = quantile(weight_kg, .975)
  ) %>%
  mutate(variable = 'waa')


tibble(waa_index) %>%
  filter(age == 3 & year < 2081) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(
      x = year,
      y = med,
      group = interaction(scenario),
      fill = scenario,
      color = scenario
    )
    ) +
      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80),
                  alpha = 0.2,
                  color = NA) +
      scale_fill_manual(values = scenPal, labels = scenLabs) +
      scale_color_manual(values = scenPal, labels = scenLabs) +
      theme(strip.background = element_blank(),
            legend.position = 'bottom') +

      guides(fill=guide_legend(nrow=1,byrow=TRUE),
             color =guide_legend(nrow=1,byrow=TRUE))+
      facet_grid(species ~ scenario ,
                 labeller = labeller(species = sppLabs,
                                     scenario =  scenLabs))+
      labs(
        x = 'Year',
        y = 'EWAA @ 50% maturity, kg',
        color = '',
        fill = ''
      )

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  theme(legend.position = 'bottom')


ggsave(last_plot(),
       file = here::here('figs','manuscript','OM_waa.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')


#* survey stations ----
ns_domain <- as.matrix(read.csv(here::here('data','northsea_latlong.csv')), ncol = 2)
survey_1 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_1.csv')) %>% filter(year == 2013)
survey_15 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_0.15.csv')) %>% filter(year == 2013)
survey_05 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_0.05.csv')) %>% filter(year == 2013)
survey_50 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_0.5.csv')) %>% filter(year == 2013)

set.seed(731)

pars_sims <- matrix(c(0,0,0.5,0,5,0,5,1),ncol = 2, byrow = TRUE)## sigma, phi by scenario

tmp <- geoR::grf(1e3,
                 grid = ns_domain,
                 cov.pars = pars_sims[4,]
)


sims <-  data.frame(tmp$coords,tmp$data)

baseplot <- ggplot(sims, aes(x = lat, y= long, fill = tmp.data, color =tmp.data)) +
  geom_raster() +
  coord_equal() +
  scale_color_viridis_c()+
  scale_fill_viridis_c()+
  theme_void() +
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme(legend.position = 'none', strip.text = element_blank())


surv1 <- baseplot +
  geom_point(data = survey_1, aes(x = lat, y = long), fill = 'blue', color = 'blue')+
  annotate(
    "text", label = "A)",
    x = 2, y = 50, size = 6
  )

surv15 <- baseplot +
  geom_point(data = survey_15, aes(x = lat, y = long), fill = 'blue', color = 'blue')+
  annotate(
    "text", label = "C)",
    x = 2, y = 50, size = 6
  )

surv05 <- baseplot +
  geom_point(data = survey_05, aes(x = lat, y = long), fill = 'blue', color = 'blue')+
  annotate(
    "text", label = "D)",
    x = 2, y = 50, size = 6
  )

surv50 <- baseplot +
  geom_point(data = survey_50, aes(x = lat, y = long), fill = 'blue', color = 'blue') +
  annotate(
    "text", label = "B)",
    x = 2, y = 50, size = 6
  )

surv1  | surv50 | surv15  | surv05

ggsave(last_plot(),
       file = here::here('figs','manuscript','survey_stations.png'),
       width = 8, height = 4, unit = 'in', dpi = 520)

#* Survey Obs for add'l runs ----
tibble(dobs) %>%
  filter(fc %in% c(1.01,0.5,0.05) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      labs(x = 'Year',y = '1000 mt', color = '', fill = '') +
      facet_grid(species ~ variable + fc,
                 # ncol = 2,
                 # scales = 'free',
                 labeller = labeller(species = sppLabs, fc = fcLabs))
  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  guides(fill=guide_legend(nrow=1,byrow=TRUE),
         color =guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') &
  scale_fill_manual(values = scenPal, labels = scenLabs) &
  scale_color_manual(values = scenPal, labels = scenLabs)


ggsave(last_plot(),
       file = here::here('figs','manuscript','Suppl_OM_bio_survey.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')
#* MRE for add'l runs ----

tibble(mre_summary) %>%
  # filter(fc %in% c(1.01,1,0.15) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey80')+

      geom_line() +
      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      labs(x = 'Year',y = 'Relative Error', color = '', fill = '') +
      facet_grid(species ~ fc,
                 labeller = labeller(species = sppLabs, fc = fcLabs))
  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  guides(fill=guide_legend(nrow=1,byrow=TRUE),
         color =guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') &
  scale_fill_manual(values = scenPal, labels = scenLabs) &
  scale_color_manual(values = scenPal, labels = scenLabs)


ggsave(last_plot(),
       file = here::here('figs','manuscript','Suppl_MRE.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')
