## Figures for manuscript & supplementary material
## assumes you have mre_all with final results (including survey data & diagnostics)

## build plotting summaries
mre_all_filt <- mre_all %>% filter(runs_test_passed)  ## drop failed diagnostics

bobs <- mre_all_filt %>%
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

sobs <- mre_all_filt %>%
  filter(!is.na(Obs)) %>%
  group_by(year, scenario, species, fc) %>%
  dplyr::summarize(
    med = median(Obs),
    lwr80 = quantile(Obs, .10),
    upr80 = quantile(Obs, .90),
    lwr95 = quantile(Obs, .0275),
    upr95 = quantile(Obs, .975)
  )

mre_summary <- mre_all_filt %>%
  group_by(year, scenario, species, fc,ewaa, q_treatment) %>%
  mutate(MRE_scaled = MRE_totbio*100) %>%
  summarize(
    med = median(MRE_scaled),
    lwr50 = quantile(MRE_scaled, .25),
    upr50 = quantile(MRE_scaled, .75),
    lwr95 = quantile(MRE_scaled, .0275),
    upr95 = quantile(MRE_scaled, .975)
  )



## Figure X1. OM SB, 1.0 Survey Obs, and EWAA by SPP ----

#* OM Biomass ----
ggplot(bobs, aes(x= year, y = med, color = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr80, ymax = upr80), alpha = 0.5, color = NA) +
  facet_grid(~species) +
  scale_fill_manual(values = scenPal, labels = scenLabs) +
  scale_color_manual(values = scenPal, labels = scenLabs) +
  facet_wrap( ~ species, scales = 'free_y', labeller = as_labeller(sppLabs)) +
  labs(
    x = 'Year',
    y = 'True SSB (kmt)',
    color = '',
    fill = ''
  ) +
  theme(legend.position = 'none')

#* survey observations, perfect coverage ----
ggplot(subset(sobs, fc == 1),   aes( x = year,   y = med ,  group = interaction(scenario),
         fill = scenario,
         color = scenario
       )) +
  geom_point() +
  geom_errorbar( aes(  ymin = med - 0.05 * med,    ymax =   med + 0.05 * med, color = scenario), alpha = 0.2,  width = 0  ) +
  scale_fill_manual(values = scenPal, labels = scenLabs) +
  scale_color_manual(values = scenPal, labels = scenLabs) +

  theme(strip.background = element_blank()  ,
        legend.position = 'none') +
  labs(x = 'Year',y = 'Survey Biomass (kmt)', color = '', fill = '') +
  facet_wrap(~ species,
             scales = 'free_y',
             labeller = as_labeller(sppLabs))

#* survey observations, 15% coverage ----
ggplot(subset(sobs, fc == 0.15),   aes( x = year,   y = med ,  group = interaction(scenario),
                                     fill = scenario,
                                     color = scenario)) +
  geom_point() +
  geom_errorbar( aes( ymin = lwr80,  ymax =  upr80, color = scenario ), alpha = 0.2, width = 0) +
  scale_fill_manual(values = scenPal, labels = scenLabs) +
  scale_color_manual(values = scenPal, labels = scenLabs) +

  theme(strip.background = element_blank()  ,
        legend.position = 'none') +
  labs(x = 'Year',y = 'Survey Biomass (kmt)', color = '', fill = '') +
  facet_wrap(~ species,
             scales = 'free_y',
             labeller = as_labeller(sppLabs))

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
  facet_grid( fc ~ species) +
  theme(axis.title = element_blank()
  )



