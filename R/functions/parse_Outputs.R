parse_Outputs <- function(files_to_run){

  ## combine MLE files
  mre_all0 <- list.files(
    files_to_run,
    pattern = 'mre.csv',
    recursive = TRUE,
    full.names = TRUE) %>%
    lapply(., FUN = read.csv) %>%
    bind_rows() %>%
    dplyr::select(-inputN) %>%
    filter(q_treatment == 'fixed') %>%
    distinct() ## drop old runs and duplicates

  ## calculate convergence
  ## IDs of fails to converge
  failed_to_converge <- mre_all0 %>%
    filter(is.na(lower)) %>%
    dplyr::select(scenario, replicate, species) %>%
    distinct()


  ## IDs of fails to execute
  failed_to_run <- summarise(mre_all0,
                             nrep = n() / 71,
                             .by = c(scenario, replicate, species)) %>%
    tidyr::complete(replicate, scenario, species) %>%
    filter(is.na(nrep)  ) %>%
    arrange(species, scenario, replicate)

  mre_all_filter0 <- mre_all0 %>%
    filter(!(replicate %in% failed_to_run$replicate)) %>% ## drop fails across-the-board
    filter(!(replicate %in% failed_to_converge$replicate)) ## drop fails across-the-board

  ## load model and do diagnostics and pull survey
  # length(unique(mre_all_filter0$replicate)) ## how many unique replicates remained
  # nrow(mre_all_filter0 %>% filter(fc %in% c(1,0.15)))/71 ## final total number of unique models run

  ## now re-build files_to_run after filter
  kept_reps <- paste0('rep',sort(unique(mre_all_filter0$replicate)) )
  files_to_run_filter1 <- files_to_run[grep(paste(kept_reps, collapse= "|"),files_to_run)]

  # Diagnostics ----
  # devtools::install_github("jabbamodel/ss3diags")
  #* run's test ----
  ## do run's test and examine >10x resid > 2

  models_filtered <- list.files(files_to_run_filter1,pattern='*_model.rdata', recursive = TRUE,
                                full.names = TRUE)[406:896]


  ## placeholders for diag pass rates
  diags <- data.frame()
  survey_exp <- data.frame() ## going to fill this at the same time

  for(i in seq_along(models_filtered)){

    ## populate data frame front matter
    reppy <- as.numeric(gsub('rep','',basename(dirname(dirname(models_filtered[i])))))     ## yank the replicate
    diags[i,'replicate'] <- reppy
    diags[i,'fc'] <- as.numeric(stringr::str_sub(basename(dirname(models_filtered[i])),21,-22))
    diags[i,'scenario']  <- strsplit(basename(dirname(dirname(dirname(models_filtered[i])))),'-')[[1]][2]
    diags[i,'species']  <- strsplit(basename(dirname(dirname(dirname(models_filtered[i])))),'-')[[1]][1]
    diags[i,'runs_test_passed'] <-  diags[i,'runs_test_val'] <- NA

    x <- try({
      load(models_filtered[i], .GlobalEnv) ## will load as mod_use
    })
    if (inherits(x, 'try-error')) next()
    rep  <- mod_use$report()
    ## do the run's test
    resids <- data.frame(Obs = mod_use$env$data$agg_indices,
                         Exp = rep$pred_indices,
                         abund_cv = mod_use$env$data$agg_index_sigma) %>%
      mutate(year = 2010:2080) %>%
      filter(Obs > 0 ) %>%
      mutate(residuals = log(Obs)-log(Exp))
    # devtools::install_github("jabbamodel/ss3diags")

    runs_test<-ss3diags::ssruns_sig3(x=as.numeric(resids$residuals),type="resid")
    diags[i,'runs_test_val']  <- runs_test$p.runs
    diags[i,'runs_test_passed']  <- ifelse(runs_test$p.runs <= 0.05, FALSE, TRUE)

    resids[,'replicate'] <- reppy
    resids[,'fc'] <- as.numeric(stringr::str_sub(basename(dirname(models_filtered[i])),21,-22))
    resids[,'scenario']  <- strsplit(basename(dirname(dirname(dirname(models_filtered[i])))),'-')[[1]][2]
    resids[,'species']  <- strsplit(basename(dirname(dirname(dirname(models_filtered[i])))),'-')[[1]][1]

    survey_exp <- bind_rows(survey_exp,resids)

    ## check the agecomp residuals, looking for > 10 outliers (>2)
    # agecomp_resid <- merge(reshape2::melt(mod_use$env$data$catch_paa[1,,] ),
    #                        reshape2::melt(rep$pred_catch_paa[,1,]) , by = c('Var1','Var2')) %>%
    #   dplyr::select(year = Var1, age = Var2, Obs = value.x, Exp = value.y) %>%
    #   mutate(pearson = (Exp-Obs)^2)
    #
    # diags[i,'age_pearson_val']  <-  sum(agecomp_resid$pearson > 2)
    # diags[i,'age_pearson_passed']  <- ifelse( diags[i,'age_pearson_val']>=10,FALSE, TRUE)
    write.csv(diags,
              file = here::here('outputs','summary_data',paste0(Sys.Date(),"-diags.csv")),
              row.names = FALSE)
    write.csv(survey_exp,
              file = here::here('outputs','summary_data',paste0(Sys.Date(),"-survey_exp.csv")),
              row.names = FALSE)
    # setTxtProgressBar(pb,i)

  }
  close(pb)
  ## bind new data frames, this is not filtered by diags
  full_df <- merge(mre_all_filter0, survey_exp,
        by = c('replicate','scenario','species', 'year', 'fc'), all.x = TRUE) %>%
    merge(., diags,  by = c('replicate','scenario','species',  'fc') )


  full_df$q_treatment <- factor(full_df$q_treatment, levels = c('fixed','estimated'))
  full_df$q_treatment[is.na(full_df$q_treatment)] <- 'estimated'

  full_df$ewaa <- factor(full_df$ewaa, levels = c('perfect','averaged'))
  full_df$ewaa[is.na(full_df$ewaa)] <- 'perfect'
  # full_df$fc <- factor(full_df$fc, levels = c(1, 0.15))
  return(full_df)
}

# full_df %>%
#   filter(runs_test_passed & !is.na(Obs)) %>%
#   mutate(abund_mean =Obs) %>%
#   group_by(year, scenario, species) %>%
#   summarize(
#     med = median(abund_mean),
#     cv_median = median(abund_cv),
#     lwr50 = quantile(abund_mean  , .25),
#     upr50 = quantile(abund_mean  , .75),
#     lwr95 = quantile(abund_mean  , .0275),
#     upr95 = quantile(abund_mean  , .975)
#   ) %>%
#
#
#
# ggplot(.,
#        aes(
#          x = year,
#          y = med ,
#          group = interaction(scenario),
#          fill = scenario,
#          color = scenario
#        )) +
#   geom_point() +
#   # {if(species == 'AtlanticHerring') scale_y_continuous(limits = c(1000,4000))} +
#   # {if(species == 'AtlanticCod') scale_y_continuous(limits = c(0,1000))} +
#   # {if(species == 'EuropeanSprat') scale_y_continuous(limits = c(1500,8000))} +
#   # geom_errorbar(aes(ymin = abund_mean-abund_mean*abund_cv,
#   #                   ymax =  abund_mean+abund_mean *abund_cv,
#   #                   color = scenario ),
#   #               alpha = 0.2, width = 0) +
#   {
#     if (fc != 1)
#       geom_errorbar(
#         aes(
#           ymin = lwr50,
#           ymax =  upr50,
#           color = scenario
#         ),
#         alpha = 0.2,
#         width = 0
#       )
#   } +
#   {
#     if (fc == 1)
#       geom_errorbar(
#         aes(
#           ymin = med - 0.05 * med,
#           ymax =   med + 0.05 * med,
#           color = scenario
#         ),
#         alpha = 0.2,
#         width = 0
#       )
#   } +
#   scale_fill_manual(values = scenPal, labels = scenLabs) +
#   scale_color_manual(values = scenPal, labels = scenLabs) +
#
#   theme(strip.background = element_blank()  ,
#         legend.position = 'none') +
#   labs(
#     x = 'Year',
#     y = 'Survey Biomass (kmt)',
#     color = '',
#     fill = ''
#   ) +
#   facet_wrap(~ species,
#              scales = 'free_y',
#              labeller = as_labeller(sppLabs)) #+
# # facet_wrap(~scenario, ncol = 1, scales = 'fixed')


