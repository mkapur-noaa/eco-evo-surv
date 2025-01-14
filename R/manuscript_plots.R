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
mre_summary$fc <- factor(mre_summary$fc,levels = c(1.01,1,0.15,.5,0.05))



## Figure 1. OM Maps ----
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


maps <- tibble(tabund_spatial) %>%
  filter(year %in% c( 2040))  %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x,aes(x = lat, y = long,
                   fill = abundance_rescale_year)) +
      theme_void() +
      geom_raster() +
      theme(strip.text.y = element_text(angle = 270))+
      facet_grid(species ~ scenario ,
                 labeller = labeller(species = sppLabs,
                                     scenario =  scenLabs))

  }) %>%
  patchwork::wrap_plots(ncol = 1, guides = 'collect') &
  # guides(fill=guide_legend(nrow=1,byrow=TRUE),
  #        color =guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') &
  labs(fill = 'relative abundance in year 2040') &
  # scale_fill_gradientn(colours = rev(scenPal)) &
  scale_fill_gradientn(colours = c("#95cec7","#2a9d8f","#e9c46a","#f4a261","#e76f51"))

ggsave(last_plot(),
       file = here::here('figs','Figure1.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')

## Figure 2. OM  WAA  ----
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




tibble(bind_rows(waa_index)) %>%
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
       file = here::here('figs','Figure2.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')


## Figure 3. OM SB & Survey Obs ----


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
      # geom_point() +
      geom_line() +

      # geom_errorbar( aes( ymin = lwr80,  ymax =  upr80, color = scenario ), alpha = 0.2, width = 0) +

      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      # scale_fill_manual(values = scenPal, labels = scenLabs) +
      # scale_color_manual(values = scenPal, labels = scenLabs) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      # theme(axis.title.y = element_blank())+
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
       file = here::here('figs','Figure3.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')



## Figure 4. MRE ----


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
      labs(x = 'Year',y = 'Median Relative Error', color = '', fill = '') +
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
       file = here::here('figs','Figure4.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')

ggplot(subset(mre_summary,  fc %in% c(0.15,1)),
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
  facet_grid( species ~ fc) +
  theme(axis.title = element_blank()
  )

## Supplementary Figures ----
#* Survey Obs for add'l runs ----

tibble(dobs) %>%
  filter(fc %in% c(1.01,0.5,0.05) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      # geom_point() +
      geom_line() +

      # geom_errorbar( aes( ymin = lwr80,  ymax =  upr80, color = scenario ), alpha = 0.2, width = 0) +

      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      # scale_fill_manual(values = scenPal, labels = scenLabs) +
      # scale_color_manual(values = scenPal, labels = scenLabs) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      # theme(axis.title.y = element_blank())+
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
       file = here::here('figs','FigureS3.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')
#* MRE for add'l runs ----

tibble(dobs) %>%
  filter(fc %in% c(1.01,1,0.15) ) %>%
  group_split(species) %>%
  purrr::map({
    ~ggplot(.x, aes( x = year,   y = med ,  group = interaction(scenario),
                     fill = scenario,
                     color = scenario)) +
      # geom_point() +
      geom_line() +

      # geom_errorbar( aes( ymin = lwr80,  ymax =  upr80, color = scenario ), alpha = 0.2, width = 0) +

      geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.2, color = NA) +
      # scale_fill_manual(values = scenPal, labels = scenLabs) +
      # scale_color_manual(values = scenPal, labels = scenLabs) +
      theme(strip.background = element_blank()) +
      theme(legend.position = 'none')+
      # theme(axis.title.y = element_blank())+
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
       file = here::here('figs','Figure3.png'),
       width = 10, height = 10, dpi = 500, unit = 'in')

