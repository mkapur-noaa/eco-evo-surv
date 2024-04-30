## build_WHAM.R
## script to take formatted data and build, execute WHAM input files

run_WHAM <-function(yrs_use = 2010:2099, ## years to run the assessment
                      file_suffix = NULL
                   ){

  filen <- strsplit(file_suffix,'-')[[1]]
  scen <- filen[5]
  repID2 <- filen[6]
  spname <- filen[4]

  wham.dir <- here::here('outputs','wham_runs',file_suffix)

  ## load input data ----
  mortality <- read.table(paste0(wham.dir,"/",file_suffix,'-wham_mortality.csv'),  skip = 1)[1:length(yrs_use),]
  maturity <-read.table(paste0(wham.dir,"/",file_suffix,'-wham_maturity.csv'),  skip = 1)[1:length(yrs_use),]
  catch_at_age <- read.table(paste0(wham.dir,"/",file_suffix,'-wham_catch_at_age.csv'),  skip = 1)[1:length(yrs_use),]
  survey <- read.table(paste0(wham.dir,"/",file_suffix,'-wham_survey.csv'),  skip = 1)[1:length(yrs_use),]
  # ncol(survey) == asap3$dat$n_ages+4 #(year, obs, cv, ss)
  waa_catch <- read.table(paste0(wham.dir,"/",file_suffix,'-wham_waa_catch.csv'),  skip = 1)[1:length(yrs_use),]
  waa_ssb <-read.table(paste0(wham.dir,"/",file_suffix,'-wham_waa_ssb.csv'),  skip = 1)[1:length(yrs_use),]
  N_init <- read.table(paste0(wham.dir,"/",file_suffix,'-wham_N_ini.csv'),  skip = 1)[2,]
  true_biomass <- read.csv(paste0(wham.dir,"/",file_suffix,'-true_biomass.csv'))

  ## load asap3-style data file ----
  # copy templated asap3 data file to working directory
  file.copy(from=file.path(here::here('data','wham_inputs','osmose_wham_template.dat')),
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
  asap3$dat$catch_cv <- matrix(0.1, nrow = asap3$dat$n_years)
  asap3$dat$catch_Neff <- matrix(100, nrow = asap3$dat$n_years)


  #* selex ----
  #* no time blocks on selex, no random effects
  #* fishery selex is age specific
  #* survey selex is logistic
  input1 <- prepare_wham_input(asap3,
                               recruit_model=2,
                               model_name=file_suffix,
                               selectivity=list(model=c('logistic','logistic'),
                                                re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
                                                initial_pars=list(c(11,0.9), ## age-specific start pars, fishery
                                                                  c(11,0.9))), ## alpha, b1, survey
                                                # fix_pars=list(10,2)), ## fix ages 1:11 for fish and b1 for survey
                               NAA_re = list(sigma="rec", cor="iid"))
  m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time

  # Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
  check_convergence(m1)

  # Save list of all fit models
  # mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)
  # save("mods", file="ex1_models.RData")

  # Compare models by AIC and Mohn's rho
  # res <- compare_wham_models(mods, table.opts=list(fname="ex1_table", sort=TRUE))
  # res$best

  # Project best model, m4,
  # Use default values: 3-year projection, use average selectivity, M, etc. from last 5 years
  # m4_proj <- project_wham(model=mods$m4)

  # WHAM output plots for best model with projections ----
  plot_wham_output(mod=m1,
                   res = 250,
                   dir.main = wham.dir) # default is png
  # plot_wham_output(mod=m4_proj, out.type='html')

  # calculate MRE, save figure and table
  std <- summary(mod$sdrep)
  ssb.ind <- which(rownames(std) == "log_SSB")[1:length(m1$years)]

  mre_table <- true_biomass[1:length(m1$years),] %>%
    mutate(ssb_est = exp(std[ssb.ind, 1]) ,
           ssb_est_cv = std[ssb.ind, 2],
           lower = ssb_est - ssb_est*ssb_est_cv,
           upper = ssb_est + ssb_est*ssb_est_cv,
           MRE = (ssb_est-abund_mean)/abund_mean,
           MRE_scaled = 100*MRE)

  ssb_compare <-  mre_table %>%
    select(year, abund_mean, ssb_est, lower, upper) %>%
    reshape2::melt(id = c('year','lower','upper')) %>%
    ggplot(., aes(x = year, y = value/1e3, color = variable)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, aes(ymin = lower/1e3, ymax = upper/1e3),
                color = NA,
                fill = 'grey22')+
    scale_color_manual(values = c(scenLabs2$Pal[scenLabs2$Var2 == scen] ,
                                  'grey22'),
                       labels = c('evOsmose Operating Model','WHAM Estimation Model')) +
    labs(x = 'Assessment Year', y = 'SSB (kmt)') +
    theme(legend.position = c(0.8,0.8)) +
    labs(color = '') +
    scale_y_continuous(limits = c(0,1e-3*1.25*max(mre_table$abund_mean)))

 mre<- ggplot(mre_table, aes(x = year, y = MRE_scaled)) +
    geom_line() +
    scale_y_continuous(limits = c(-50,50)) +
    labs(x = 'Assessment Year', y = 'MRE SSB, %')+
    geom_hline(yintercept = 0, color = 'pink')

  write.csv(mre_table,paste0(wham.dir,"/",file_suffix,"-ssb_mre.csv"), row.names = FALSE)

  png(file =  paste0(wham.dir,'/plots_png','/results','/ssb_mre.png'),
      height = 5, width = 12, unit = 'in',res = 520)
  Rmisc::multiplot(ssb_compare, mre, cols = 2)
  dev.off()


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
