## build_WHAM.R
## script to take formatted data and build, execute WHAM input files

run_WHAM <-function(yrs_use = 2010:2099, ## years to run the assessment
                    fractional_coverage_use = 1,
                    file_suffix = NULL
){

  filen2 <- basename(dirname(dirname(file_suffix)))
  scen <-   strsplit(filen2,'-')[[1]][2]
  repID2 <- as.numeric(gsub('rep','',basename(file_suffix)))
  spname <-   strsplit(filen2,'-')[[1]][1]

  file_suffix2 <- paste0(filen2,'-',repID2) ## append to filenames
  wham.dir <- file_suffix ## where to read things from

  ## subfolder with fractional coverage
  filen3 <- ifelse(fractional_coverage_use==1,
                   'perfect_information',
                   paste0("fractional_coverage=", fractional_coverage_use))
  wham.dir.save <- paste0(file_suffix,"/",filen3); if(!dir.exists(wham.dir.save)) dir.create(wham.dir.save)

  ## load input data ----
  mortality <- read.table(paste0(wham.dir,"/",file_suffix2,'-wham_mortality.csv'),  skip = 1)[1:length(yrs_use),]
  maturity <-read.table(paste0(wham.dir,"/",file_suffix2,'-wham_maturity.csv'),  skip = 1)[1:length(yrs_use),]
  catch_at_age <- read.table(paste0(wham.dir,"/",file_suffix2,'-wham_catch_at_age.csv'),  skip = 1)[1:length(yrs_use),]
  survey <- read.table(paste0(wham.dir,"/",file_suffix2,'-',fractional_coverage_use,'-wham_survey.csv'),  skip = 1)[1:length(yrs_use),]
  # ncol(survey) == asap3$dat$n_ages+4 #(year, obs, cv, ss)
  waa_catch <- read.table(paste0(wham.dir,"/",file_suffix2,'-wham_waa_catch.csv'),  skip = 1)[1:length(yrs_use),]
  waa_ssb <-read.table(paste0(wham.dir,"/",file_suffix2,'-wham_waa_ssb.csv'),  skip = 1)[1:length(yrs_use),]
  N_init <- read.table(paste0(wham.dir,"/",file_suffix2,'-wham_N_ini.csv'),  skip = 1)[2,]
  true_biomass <- read.csv(paste0(wham.dir,"/",file_suffix2,'-true_biomass.csv'))

  ## load asap3-style data file ----
  # copy templated asap3 data file to working directory
  file.copy(from=file.path(here::here('osmose_wham_template.dat')),
            to=wham.dir, overwrite=TRUE)
  # read asap3 data file and convert to input list for wham
  setwd(wham.dir)
  asap3 <- read_asap3_dat("osmose_wham_template.dat") ## uncorrected template

  ## populate asap3 data file ----
  #* dimensions and spex ----
  asap3$dat$n_years <- nrow(maturity)
  asap3$dat$nfinalyear <- max(yrs_use)
  asap3$dat$year1 <- 2010
  asap3$dat$n_ages <- ncol(maturity)
  asap3$dat$n_fleets <- 1 ## one fishing fleet
  asap3$dat$n_fleet_sel_blocks  <- 1 ## one fishing fleet
  asap3$dat$n_indices <- 1 ## one survey fleet
  asap3$dat$n_WAA_mats <- 2 ## SSB and catch WAA
  asap3$dat$WAA_pointers <- c(1,2,1,2,2,2) #catch, discards, total catch, total discards, ssb, jan1bio
  asap3$dat$index_units <- 1 ## 1 = biomass, 2 = numbers
  asap3$dat$index_acomp_units <- 2
  asap3$dat$index_month <- 7
  asap3$dat$use_index_acomp <- 1
  asap3$dat$use_index <- 1
  asap3$dat$sel_block_assign <- matrix(1,asap3$dat$n_years )

  #* initial pars ----
  asap3$dat$N1_ini <- as.matrix(N_init, nrow = 1)
  ## asap3$dat$F1_ini

  #* data ----
  asap3$dat$maturity <- as.matrix(maturity)
  asap3$dat$M  <-  as.matrix(mortality)
  asap3$dat$WAA_mats[[1]] <- as.matrix(waa_catch)
  asap3$dat$WAA_mats[[2]] <- as.matrix(waa_ssb)
  asap3$dat$CAA_mats[[1]] <- as.matrix(catch_at_age)
  asap3$dat$DAA_mats[[1]] <- matrix(0, nrow = asap3$dat$n_years, ncol = ncol(catch_at_age))
  asap3$dat$prop_rel_mats[[1]] <- matrix(0, nrow = asap3$dat$n_years, ncol = asap3$dat$n_ages)
  asap3$dat$IAA_mats[[1]] <- as.matrix(survey)
  asap3$dat$catch_cv <- matrix(0.05,asap3$dat$n_years )
  asap3$dat$catch_Neff <- matrix(100, nrow = asap3$dat$n_years)

  #* selex ----
  #* no time blocks on selex, no random effects
  # input1 <- prepare_wham_input(asap3,
  #                              recruit_model=2, ## (default) Random about mean, i.e. steepness = 1
  #                              model_name=file_suffix2,
  #                              selectivity=list(model=c('logistic','age-specific'),
  #                                               re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                                               initial_pars=list(c(7,0.9), ## fully selected fishery
  #                                                                 rep(1, asap3$dat$n_ages)), ## fully selected survey
  #                                               fix_pars=list(NULL,
  #                                                             c(1:asap3$dat$n_ages))), ## fix 'em all
  #                              # selectivity=list(model=c('logistic','logistic'),
  #                              #                  re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                              #                  initial_pars=list(c(7,0.9), ## age-specific start pars, fishery
  #                              #                                    c(7,0.9))), ## alpha, b1, survey
  #                              NAA_re = list(sigma="rec", cor="iid"))
  # m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time
  # # exp(m1$par[grep('N1',names(m1$par))])
  # # Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
  # check_convergence(m1)

  input2 <- prepare_wham_input(asap3,
                               recruit_model=2, ## (default) Random about mean, i.e. steepness = 1
                               model_name=file_suffix2,
                               selectivity=list(model=c('double-logistic','age-specific'),
                                                re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
                                                initial_pars=list(c(7,0.9,7,0.9), ## fully selected fishery
                                                                  rep(1, asap3$dat$n_ages)), ## fully selected survey
                                                fix_pars=list(NULL,
                                                              c(1:asap3$dat$n_ages))), ## fix 'em all
                               # selectivity=list(model=c('logistic','logistic'),
                               #                  re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
                               #                  initial_pars=list(c(7,0.9), ## age-specific start pars, fishery
                               #                                    c(7,0.9))), ## alpha, b1, survey
                               NAA_re = list(sigma="rec", cor="iid"))
  m2 <- fit_wham(input2, do.osa = F) # turn off OSA residuals to save time
  check_convergence(m2)

  # input3 <- prepare_wham_input(asap3,
  #                              recruit_model=2, ## (default) Random about mean, i.e. steepness = 1
  #                              model_name=file_suffix2,
  #                              selectivity=list(model=c('age-specific','age-specific'),
  #                                               re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                                               initial_pars=list(rep(0.5, asap3$dat$n_ages), ## fully selected fishery
  #                                                                 rep(1, asap3$dat$n_ages)), ## fully selected survey
  #                                               fix_pars=list(NULL,
  #                                                             c(1:asap3$dat$n_ages))), ## fix 'em all
  #                              # selectivity=list(model=c('logistic','logistic'),
  #                              #                  re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                              #                  initial_pars=list(c(7,0.9), ## age-specific start pars, fishery
  #                              #                                    c(7,0.9))), ## alpha, b1, survey
  #                              NAA_re = list(sigma="rec", cor="iid"))
  # m3 <- fit_wham(input3, do.osa = F) # turn off OSA residuals to save time
  # check_convergence(m3)

  # input4 <- prepare_wham_input(asap3,
  #                              recruit_model=2, ## (default) Random about mean, i.e. steepness = 1
  #                              model_name=file_suffix2,
  #                              selectivity=list(model=c('age-specific','age-specific'),
  #                                               re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                                               initial_pars=list(rep(0.5, asap3$dat$n_ages), ## fully selected fishery
  #                                                                 rep(1, asap3$dat$n_ages)), ## fully selected survey
  #                                               fix_pars=list(NULL,
  #                                                            NULL)),
  #                              # selectivity=list(model=c('logistic','logistic'),
  #                              #                  re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
  #                              #                  initial_pars=list(c(7,0.9), ## age-specific start pars, fishery
  #                              #                                    c(7,0.9))), ## alpha, b1, survey
  #                              NAA_re = list(sigma="rec", cor="iid"))
  # m4 <- fit_wham(input4, do.osa = F) # turn off OSA residuals to save time
  # check_convergence(m4)

  # Save list of all fit models
  # mods <- list(m1=m1, m2=m2)
  # save("mods", file=  paste0(wham.dir,"/",file_suffix2,"-all_models.Rdata")

  # Compare models by AIC and Mohn's rho
  # res <- compare_wham_models(mods, table.opts=list(fname="ex1_table", sort=TRUE))
  # res$best

  mod_use <- m2

  # Project best model, m4,
  # Use default values: 3-year projection, use average selectivity, M, etc. from last 5 years
  # m4_proj <- project_wham(model=mods$m4)


  if(mod_use$is_sdrep){
    std <- summary(mod_use$sdrep)
    ssb.ind <- which(rownames(std) == "log_SSB")[1:length(mod_use$years)]
    mre_table <- true_biomass[1:length(mod_use$years),] %>%
      mutate(ssb_est = exp(std[ssb.ind, 1]) ,
             ssb_est_cv = std[ssb.ind, 2],
             lower = ssb_est - ssb_est*ssb_est_cv,
             upper = ssb_est + ssb_est*ssb_est_cv,
             MRE = (ssb_est-ssb_true)/ssb_true,
             MRE_scaled = 100*MRE)
  } else{
    ssb.ind <- mod_use$report()$SSB[1:length(mod_use$years)]
    mre_table <- true_biomass[1:length(mod_use$years),] %>%
      mutate(ssb_est = ssb.ind ,
             ssb_est_cv =NA,
             lower = ssb_est ,
             upper = ssb_est,
             MRE = (ssb_est-ssb_true)/ssb_true,
             MRE_scaled = 100*MRE)
  }

  cat(mean(mre_table$ssb_est/mre_table$ssb_true),"\n")

  ssb_compare <-  mre_table %>%
    select(year, ssb_true, ssb_est, lower, upper) %>%
    reshape2::melt(id = c('year','lower','upper')) %>%
    ggplot(., aes(x = year, y = value/1e3, color = variable)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, aes(ymin = lower/1e3, ymax = upper/1e3),
                color = NA,
                fill = 'grey22')+
    scale_color_manual(values = c(scenLabs2$Pal[scenLabs2$Var2 == scen] ,
                                  'grey22'),
                       labels = c('evOsmose Operating Model','WHAM Estimation Model')) +
    labs(x = 'Year', y = 'SSB (kmt)') +
    theme(legend.position.inside = c(0.8,0.8)) +
    labs(color = '') +
    scale_y_continuous(limits = c(0,1e-3*1.25*max(mre_table$ssb_true)))

  mre<- ggplot(mre_table, aes(x = year, y = MRE_scaled)) +
    geom_line() +
    # scale_x_continuous(limits = c(2040,2099))+
    scale_y_continuous(limits = c(-100,100)) +
    labs(x = 'Year', y = 'MRE SSB, %')+
    geom_hline(yintercept = 0, color = 'pink')

  write.csv(mre_table,paste0(wham.dir.save,"/",file_suffix2,"-ssb_mre.csv"), row.names = FALSE)

  png(file =  paste0(wham.dir.save,"/",file_suffix2,"-ssb_mre.png"),
      height = 5, width = 12, unit = 'in',res = 520)
  Rmisc::multiplot(ssb_compare, mre, cols = 2)
  dev.off()

  ggsave(ssb_compare, file = paste0(wham.dir,"/",file_suffix,"-ssb_compare.png"),
         width = 6, height = 6, unit = 'in', dpi = 400)

  save(mod_use, file = paste0(wham.dir.save,'/model.rdata'))
  # WHAM output plots for best model with projections ----
  plot_wham_output(mod=mod_use,
                   res = 250,
                   dir.main = wham.dir.save) # default is png
  # plot_wham_output(mod=m4_proj, out.type='html')




}


# write.table(asap3$dat, file= "update_asap.dat")
#write.table(matrix(0.2, nrow = asap3$dat$n_years, ncol = asap3$dat$n_ages),
#            sep = ' ',
#            file = here::here('data','wham_inputs','mortality_26.csv'), row.names = FALSE)

#write.table(matrix(c(asap3$dat$maturity[1,1:5], rep(0.9999, asap3$dat$n_ages-5)),
#   byrow = TRUE,
#           nrow = asap3$dat$n_years,
#              ncol = asap3$dat$n_ages),
#     sep = ' ',
#        file = here::here('data','wham_inputs','maturity_26.csv'), row.names = FALSE)

#write.table(matrix(cbind(asap3$dat$WAA_mats[[1]],
#                      matrix(rep(asap3$dat$WAA_mats[[1]][,6],asap3$dat$n_ages-5),
#                           ncol = asap3$dat$n_ages)),
#                 ncol = asap3$dat$n_ages),
#       sep = ' ',
#         file = here::here('data','wham_inputs','waa_17.csv'), row.names = FALSE)

#write.table(round(matrix(byrow = TRUE, cbind(asap3$dat$CAA_mats[[1]][,1:6],
#   matrix(rep(asap3$dat$CAA_mats[[1]][,6],11),
#                  byrow = TRUE,
#         ncol = 11),
#            asap3$dat$CAA_mats[[1]][,7]),
#          ncol = 1+asap3$dat$n_ages)),

#  sep = ' ',
##      file = here::here('data','wham_inputs','caa_17.csv'), row.names = FALSE)

#write.table(matrix(0, ncol = 1+asap3$dat$n_ages, nrow = asap3$dat$n_years),
#      sep = ' ',
#file = here::here('data','wham_inputs','discards_17.csv'), row.names = FALSE)
#
# asap3$dat$IAA_mats
# asap3$dat$index_acomp_units
# asap3$dat$index_sel_end_age
# asap3$dat$n_ages
