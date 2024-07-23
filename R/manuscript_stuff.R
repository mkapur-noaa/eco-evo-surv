
rm(list = ls()) ## clear workspace
invisible(lapply(list.files(
  here::here('R', 'functions'), full.names = TRUE
), FUN = source)) ## load all functions and presets


## Supplement X - Survey details ----
ns_domain <- as.matrix(read.csv(here::here('data','northsea_latlong.csv')), ncol = 2)
survey_1 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_1.csv')) %>% filter(year == 2013)
survey_15 <- read.csv(here::here('outputs','wham_runs','2024-05-08-survey_array_0.15.csv')) %>% filter(year == 2013)
set.seed(731)

pars_sims <- matrix(c(0,0,0.5,0,5,0,5,1),ncol = 2, byrow = TRUE)## sigma, phi by scenario

tmp <- grf(1e3,
           grid = ns_domain,
           cov.pars = pars_sims[4,]
           )


sims<-  data.frame(tmp$coords,tmp$data)

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
  geom_point(data = survey_15, aes(x = lat, y = long), fill = 'blue', color = 'blue')

surv15 <- baseplot +
  geom_point(data = survey_1, aes(x = lat, y = long), fill = 'blue', color = 'blue')

require(patchwork)
surv1  | surv15

ggsave(last_plot(),
       file = here('figs','manuscript','supplement_1.png'),
       width = 4, height = 4, unit = 'in', dpi = 520)

# assign(value = tmp, x =paste0('sim',isim))
