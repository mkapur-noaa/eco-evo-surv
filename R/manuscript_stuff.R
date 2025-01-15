
rm(list = ls()) ## clear workspace
invisible(lapply(list.files(
  here::here('R', 'functions'), full.names = TRUE
), FUN = source)) ## load all functions and presets

## RESULTS ----
#* Data Generation ----
full_df %>%
  filter(!is.na(abund_cv) & fc %in% c(1,0.15)) %>%
  summarise(mean(abund_cv), .by = c('fc','species')) %>%
  arrange(fc, species)

#* Assessment ----
#** diagnostics ----
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

## diagnostic pass rates
diags %>%
  filter(!is.na(runs_test_passed)) %>%
  summarise(sum(runs_test_passed)/n(),.by = c('fc'))


tt0 <- full_df0 %>%
  filter(!(fc %in% c(0.5,0.05))) %>%
  summarise(n=length(unique(replicate)),.by = c('fc','species','scenario'))
sum(tt0$n)
tt <- full_df %>%
  filter(!(fc %in% c(0.5,0.05))) %>%
  summarise(n=length(unique(replicate)),.by = c('fc','species','scenario'))
sum(tt$n)


#** performance metrics ----

## MARE by species x fc
full_df %>%
  filter(!(fc %in% c(0.5,0.05)) & !is.na(MRE_scaled) & year < 2060) %>%
  summarise(median(MRE_scaled),.by = c('fc','species')) %>%
  arrange(fc, species)


# all_mods_use  <- list.files(
#   files_to_run,
#   pattern = 'model.rdata',
#   recursive = TRUE,
#   full.names = TRUE ) %>%
#   lapply(., FUN = load, .GlobalEnv) %>%
#   bind_rows() %>%
#   dplyr::select(-inputN) %>%
#   filter(q_treatment == 'fixed' & ewaa == 'perfect' & fc != 0.15001 ) %>%
#   distinct() ## drop old runs and duplicates
#
# inpt <- mod_use$input
# rpt <- mod_use$report()
#
# pred_obs_surv <- inpt$data$obs %>%
#   filter(type == 'logindex') %>%
#   mutate(year = year+2009,
#          yt = val) %>%
#   merge(.,data.frame(yhat = rpt$pred_log_indices) %>%
#            mutate(year = 2010:2080) %>%
#            filter(year %% 2 == 0), by = 'year') %>%
#   dplyr::select(year, yt, yhat) %>%
#   mutate(resid = (yhat - yt)^2)
#
# sqrt(sum(pred_obs_surv$resid)/nrow(pred_obs_surv))



## SUPPLEMENTARY MATERIAL ----
#* Extra OM info ----


## uncertainty in add'l runs
full_df %>%
  filter(!is.na(abund_cv) & fc %in% c(0.5,0.05)) %>%
  summarise(mean(abund_cv), .by = c('fc','species')) %>%
  arrange(fc, species)
## replicate 8 failed for herring no-no 0.05; drop everywhere for sens

## original dims
tt0 <- full_df0 %>%
  filter(fc %in% c(0.5,0.05)) %>%
  summarise(n=length(unique(replicate)),.by = c('fc','species','scenario'))
sum(tt0$n)

## filter diags and drop rep 8
tt <- full_df %>%
  filter(fc %in% c(0.5,0.05) & replicate != 8 ) %>%
  summarise(n=length(unique(replicate)),.by = c('fc','species','scenario'))
sum(tt$n) ## final total



