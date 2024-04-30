## build_WHAM.R
## script to take formatted data and build, execute WHAM input files

build_WHAM <-function(scenario,
                      sppIdx,
                      repID = 1,
                      file_suffix  = NULL,
                      yrs_use = 2010:2099){
  # scen <- scenLabs2[scenario,2]
  # repID2 <- repID-1
  # spname <- sppLabs2[sppIdx,2]
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


  ## load asap3 data file ----
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
  mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)
  save("mods", file="ex1_models.RData")

  # Compare models by AIC and Mohn's rho
  res <- compare_wham_models(mods, table.opts=list(fname="ex1_table", sort=TRUE))
  res$best

  # Project best model, m4,
  # Use default values: 3-year projection, use average selectivity, M, etc. from last 5 years
  m4_proj <- project_wham(model=mods$m4)

  # WHAM output plots for best model with projections
  plot_wham_output(mod=m1,res = 250, dir.main = wham.dir) # default is png
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
