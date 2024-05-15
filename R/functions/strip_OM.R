## A one-time-use function to strip the TRUE values from
## the OSMOSE runs; no need to re-generate all the base files
## upon resampling

strip_OM <- function(scenario,
                     sppIdx,
                     repID = 1,
                     yrs_use = 2010:2099,
                     units_scalar = 1){

  ## string designators
  scen <- scenLabs2[scenario,2]
  repID2 <- sort(as.character(0:28))[repID]
  spname <- sppLabs2[sppIdx,2]

  file_suffix <- paste(spname,scen,repID2,sep = '-')

  ## where the raw evOsmose outputs are stored
  dirtmp <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/",scen,'/output')

  ## where the WHAM outputs are to be stored
  ## Species-Scenario head folder
  head.dir <- here::here('outputs','wham_runs',paste(spname,scen, sep = '-')); if(!dir.exists(head.dir)) dir.create(head.dir)
  if(!dir.exists(here::here(head.dir,Sys.Date()))) dir.create(here::here(head.dir,Sys.Date()))
  wham.dir <- here::here(head.dir,Sys.Date(),paste0('rep',repID2)); if(!dir.exists(wham.dir)) dir.create(wham.dir)


  ## load parameters for this species ----
  lw_pars <- read.csv(here::here('outputs','wham_runs','length2weight.csv')) %>%  filter(species == spname)
  max_age_pop <- read.csv(here::here('outputs','wham_runs','max_age.csv')) %>%  filter(species == spname)
  vonBpars <- read.csv(here::here('outputs','wham_runs','vonBpars.csv')) %>%  filter(species == spname)

  #* populate length-at-age vector
  laa <- data.frame(length_cm = NA, age = NA, len_bin = NA)
  for(a in 1:max_age_pop$value){
    laa[a,'age'] <- a
    laa[a,'length_cm'] <- vonBpars$lInf*(1-exp(-vonBpars$K*(a-vonBpars$t0)))
  }
  write.table(laa,
              sep = ' ',
              paste0(wham.dir,"/",file_suffix,'-laa.csv'),
              row.names = FALSE)

  max_age_pop <- as.numeric(max_age_pop$value)

  #* unfished NAA (for N1_ini) ----
  ## ballpark from earlier runs
  age_spatial_nofish_path <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/one_sim_without_fishing",
                                    "/ns_spatial_abundancebyAge-",spname,"_Simu0.nc")
  abundance_nofish <- ncvar_get(nc_open(age_spatial_nofish_path),"abundance") ## this might need to switch with units_use
  reshape2::melt(abundance_nofish) %>%
    filter(!is.na(value)) %>% ## drop land
    mutate(year = 2010 + (Var4 - 1) %/% 24,
           month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(year == 2010) %>% ## assume equilibrium
    group_by(age) %>%
    summarise(value = mean(value)) %>% ## average
    ungroup() %>%
    filter(age <= max_age_pop) %>%
    t() %>%
    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-wham_N_ini.csv'),
                row.names = FALSE)


  mort_path <- paste0(dirtmp, '/mortality/', 'ns_mortalityRate-',spname,"_Simu",repID2,".csv")
  mort_csv <- read.csv(mort_path, skip = 3, header = F)[,c(1,12,18)] %>% ## timestep, Frecruits, Mrecruits
    mutate(year = floor(as.numeric(stringr::str_split_fixed(V1, "\\.", 1)) -70+2010)) %>%
    filter(year != 2100) %>%
    group_by(year) %>%

    summarise(Frecruits = round(sum(V12),4),
              Mrecruits = round(sum(V18),4)) %>%
    ungroup()


  matrix(rep(mort_csv$Mrecruits, max_age_pop),
         byrow=F, nrow = length(2010:2099)) %>%

    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-wham_mortality.csv'),
                row.names = FALSE)

  ## maturity (year x age) ----

  mat_path <- paste0(dirtmp, '/ageIndicators/', 'ns_maturityDistribByAge',"_Simu",repID2,".csv")
  maturity_data <- read.csv(mat_path, skip = 1, header = T) %>%
    reshape2::melt(id = c('Time','Age')) %>%
    filter(variable == spname & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
    mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
    group_by(year,Age) %>%
    summarise(value = round(mean(value),3)) %>% ## average maturity across year
    ungroup() %>%
    tidyr::complete(year = 2010:2099, Age = 1:max_age_pop,
                    fill = list(value = 1)) %>%
    tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = value) %>%
    select(-year)
  write.table(maturity_data,
              sep = ' ',
              paste0(wham.dir,"/",file_suffix,'-wham_maturity.csv'),
              row.names = FALSE)



  age_spatial_path <- list.files(dirtmp,
                                 pattern = paste0('spatial_abundancebyAge-',spname),
                                 recursive = T,
                                 full.names = TRUE)[repID]
  biom_spatial_path <- list.files(dirtmp,
                                  pattern = paste0('spatial_biomassbyAge-',spname),
                                  recursive = T,
                                  full.names = TRUE)[repID]

  abundance0 <- ncvar_get(nc_open(age_spatial_path),"abundance") ## this might need to switch with units_use
  ## redefine the timesteps into real year, week, and not NA
  ## and then only take out the July populations
  abundance <- reshape2::melt(abundance0) %>%
    filter(!is.na(value) & Var3 <= max_age_pop) %>% ## drop land
    mutate(year = 2010 + (Var4 - 1) %/% 24,
           month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7 & year %in% yrs_use) %>% ## July survey
    group_by(year, age, lat, long, month) %>%
    summarise(value = mean(value)) %>% ## average over the month
    ungroup() %>%
    select(-month)



  # define maximum age above which all entries are NA
  max_age_survey <- abundance %>%
    group_by(age) %>%
    summarise(all_zero = all(value == 0)) %>%
    filter(!all_zero) %>%
    summarise(max_age = max(age))
  max_age_survey <- as.numeric(max_age_survey)

  total_area <- 632 ## total number of marine cells, aka dim (all_cells)
  biomass0 <- ncvar_get(nc_open(biom_spatial_path),"biomass")
  biomass <- reshape2::melt(biomass0) %>%
    filter(!is.na(value) & Var3 <= max_age_pop) %>% ## drop land
    mutate(year = 2010 + (Var4 - 1) %/% 24,
           month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7 & year %in% yrs_use) %>% ## July survey
    group_by(year, age, lat, long, month) %>%
    summarise(value = mean(value)) %>% ## average over the month
    ungroup() %>%
    select(-month)

  #* true biomass, SSB ----
  true_biomass <- biomass %>%
    ## join with maturity data
    merge(.,   maturity_data %>%
            mutate(year = 2010:2100) %>%
            reshape2::melt(id = 'year', value.name = 'maturity'),
          by.x = c('year','age'), by.y = c('year','variable') ) %>%
    ## calculate mature biomass by year, age
    mutate(ssb_age = value*maturity) %>%
    ## collapse across age/lat/long per year
    summarise(ssb_true = sum(ssb_age),
              total_biomass=sum(value,na.rm = T),
              total_biomass_sd = sd(value, na.rm = T)*n()/sqrt(n()),
              total_biomass_cv = total_biomass_sd/total_biomass, .by = year) %>%
    mutate(abund_mean_rescale= rescale(total_biomass, to = c(0,1)),
           replicate = repID2,
           scenario = scen,
           species = spname)

  write.csv(true_biomass,
            paste0(wham.dir,"/",file_suffix,"-true_biomass.csv"),
            row.names = FALSE)


}
