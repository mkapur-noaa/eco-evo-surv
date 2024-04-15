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

  # copy templated asap3 data file to working directory
  file.copy(from=file.path(here::here('data','wham_inputs','osmose_wham_template.dat')), to=write.dir, overwrite=TRUE)
  # read asap3 data file and convert to input list for wham
  asap3 <- read_asap3_dat("osmose_wham_template.dat")
}
