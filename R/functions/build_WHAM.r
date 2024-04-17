## build_WHAM.R
## script to take formatted data and build, execute WHAM input files

build_WHAM <-function(scenario,
                      sppIdx,
                      repuse = 1,
                      yrs_use = 2010:2099){
  # create directory for analysis, E.g.,
  write.dir <- here::here('outputs','wham_runs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],"-",Sys.Date()))
  if(!dir.exists(write.dir)) dir.create(write.dir)
  setwd(write.dir)

  # wham.dir <- find.package("wham")
  # file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=TRUE)
    # file.copy(from=file.path(here::here('data','wham_inputs','osmose_wham_template.dat')), to=write.dir, overwrite=TRUE)


  # copy templated asap3 data file to working directory
  # file.copy(from=file.path(here::here('data','wham_inputs','osmose_wham_template.dat')), to=write.dir, overwrite=TRUE)
  file.copy(from=file.path(here::here('data','wham_inputs','ex1_snemayt.dat')), to=write.dir, overwrite=TRUE)
  # read asap3 data file and convert to input list for wham
  # asap3 <- read_asap3_dat("osmose_wham_template.dat")
  asap3 <- read_asap3_dat("ex1_snemayt.dat")


  asap3$dat$n_ages <- 17
  asap3$dat$maturity <- matrix(0.2, nrow = asap3$dat$n_years, ncol = asap3$dat$n_ages)


  input1 <- prepare_wham_input(asap3,
                               recruit_model=2,
                               model_name="Ex 1: SNEMA Yellowtail Flounder",
                               selectivity=list(model=rep("age-specific",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
                                                re=rep("none",asap3$dat$n_fleet_sel_blocks + asap3$dat$n_indices),
                                                initial_pars=list(c(rep(0.5,17)),
                                                                  c(rep(0.5,17))),
                                                fix_pars=list(4:5,2:4)), ## list of length = number of sel blocks
                               NAA_re = list(sigma="rec", cor="iid"))
  m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time
  # Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradiet should be < 1e-06)
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
  plot_wham_output(mod=m4) # default is png
  # plot_wham_output(mod=m4_proj, out.type='html')


}


# write.table(asap3$dat, file= "update_asap.dat")
# write.table(matrix(0.2, nrow = asap3$dat$n_years, ncol = asap3$dat$n_ages),
#             sep = ' ',
#             file = here::here('data','wham_inputs','mortality_17.csv'), row.names = FALSE)
#
# write.table(matrix(c(asap3$dat$maturity[1,1:5], rep(0.9999, 12)),
#                  byrow = TRUE,
#                  nrow = asap3$dat$n_years,
#                  ncol = asap3$dat$n_ages),
#           sep = ' ',
#           file = here::here('data','wham_inputs','maturity_17.csv'), row.names = FALSE)
#
# write.table(matrix(cbind(asap3$dat$WAA_mats[[1]],
#                        matrix(rep(asap3$dat$WAA_mats[[1]][,6],11), ncol = 11)), ncol = asap3$dat$n_ages),
#           sep = ' ',
#           file = here::here('data','wham_inputs','waa_17.csv'), row.names = FALSE)
#
# write.table(round(matrix(byrow = TRUE, cbind(asap3$dat$CAA_mats[[1]][,1:6],
#                        matrix(rep(asap3$dat$CAA_mats[[1]][,6],11),
#                               byrow = TRUE,
#                               ncol = 11),
#                        asap3$dat$CAA_mats[[1]][,7]),
#                    ncol = 1+asap3$dat$n_ages)),
#
#             sep = ' ',
#           file = here::here('data','wham_inputs','caa_17.csv'), row.names = FALSE)
#
# write.table(matrix(0, ncol = 1+asap3$dat$n_ages, nrow = asap3$dat$n_years),
#             sep = ' ',
#           file = here::here('data','wham_inputs','discards_17.csv'), row.names = FALSE)
#
# asap3$dat$IAA_mats
# asap3$dat$index_acomp_units
# asap3$dat$index_sel_end_age
# asap3$dat$n_ages
