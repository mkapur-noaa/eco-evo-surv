## Figures for manuscript & supplementary material
## assumes you have mre_all with final results (including survey data & diagnostics)

## build plotting summaries
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
    lwr50 = quantile(MRE_scaled, .25),
    upr50 = quantile(MRE_scaled, .75),
    lwr95 = quantile(MRE_scaled, .0275),
    upr95 = quantile(MRE_scaled, .975)
  )



## Figure X1. OM Maps, WAA and M


## Figure X1. OM SB, 1.0 Survey Obs, and EWAA by SPP ----


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



## MRE, perfect coverage
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



