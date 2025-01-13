
rm(list = ls()) ## clear workspace
invisible(lapply(list.files(
  here::here('R', 'functions'), full.names = TRUE
), FUN = source)) ## load all functions and presets

## RESULTS ----
#* Assessment ----
## convergence/execution rates, uses mre_all0 and mre_all from master.R
total_runs <- length(unique(mre_all0$replicate)) *
  length(unique(mre_all0$species)) *
  length(unique(mre_all0$scenario)) *
  length(unique(mre_all0$fc))

nrow(failed_to_run)/total_runs
nrow(failed_to_converge)/total_runs

## post filter
length(unique(mre_all_filter$replicate)) *
  length(unique(mre_all_filter$species)) *
  length(unique(mre_all_filter$scenario)) *
  length(unique(mre_all_filter$fc))

## diagnostics
all_mods_use  <- list.files(
  files_to_run,
  pattern = 'model.rdata',
  recursive = TRUE,
  full.names = TRUE ) %>%
  lapply(., FUN = load, .GlobalEnv) %>%
  bind_rows() %>%
  dplyr::select(-inputN) %>%
  filter(q_treatment == 'fixed' & ewaa == 'perfect' & fc != 0.15001 ) %>%
  distinct() ## drop old runs and duplicates

inpt <- mod_use$input
rpt <- mod_use$report()

pred_obs_surv <- inpt$data$obs %>%
  filter(type == 'logindex') %>%
  mutate(year = year+2009,
         yt = val) %>%
  merge(.,data.frame(yhat = rpt$pred_log_indices) %>%
           mutate(year = 2010:2080) %>%
           filter(year %% 2 == 0), by = 'year') %>%
  dplyr::select(year, yt, yhat) %>%
  mutate(resid = (yhat - yt)^2)

sqrt(sum(pred_obs_surv$resid)/nrow(pred_obs_surv))



## SUPPLEMENTARY MATERIAL ----
#* Extra OM info ----


#*Survey details ----
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
    "text", label = "B)",
    x = 2, y = 50, size = 6
  )

surv05 <- baseplot +
  geom_point(data = survey_05, aes(x = lat, y = long), fill = 'blue', color = 'blue')+
  annotate(
    "text", label = "C)",
    x = 2, y = 50, size = 6
  )

surv50 <- baseplot +
  geom_point(data = survey_50, aes(x = lat, y = long), fill = 'blue', color = 'blue') +
  annotate(
    "text", label = "D)",
    x = 2, y = 50, size = 6
  )

require(patchwork)
surv1  | surv50 | surv15  | surv05

ggsave(last_plot(),
       file = here('figs','manuscript','supplement_1.png'),
       width = 8, height = 4, unit = 'in', dpi = 520)

# assign(value = tmp, x =paste0('sim',isim))
