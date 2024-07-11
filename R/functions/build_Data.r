## build_data.R
## script to take outputs from evOSMOSE model runs and build a data frame for analysis in WHAM
## UNIT assumption is single scenario, single species, single replicate
## also save some intermediate summaries and figures
## M S Kapur 2024-04-11


## scenario: 1 = noCC_noEvo, 2 = noCC_Evo, 3 = CC_noEvo, 4 = CC_Evo
## sppIdx: 1-15
## replicate: 1-28; doesn't necessarily correspond to the file name
## yrs_use: years to extract data for
## obs_error: whether or not survey biomass has observation error (NULL or value)
## units: 'numbers' or 'biomass'
## units_scalar: optional factor to decrease the scale of the index data (applies to comps as well)

build_Data<-function(scenario,
                     sppIdx,
                     repID = 1,
                     yrs_use = 2010:2080,
                     srv_selex = 11, ## whether or not survey sampled with age-based selex (logistic with A50=11)
                     obs_error = 0.2, ## whether or not survey sampled with observation error
                     fractional_coverage_use = 1, ## fractional coverage of survey
                     units = 'biomass', ## units for survey observations
                     units_scalar = 1,
                     date.use = NULL, ## if you want to backfill an older folder
                     do_GAM = FALSE) ## whether or not to invoke spatial standardization
{


  ## string designators
  scen <- scenLabs2[scenario,2]
  repID2 <- repID-1 #sort(as.character(0:27))[repID]
  spname <- sppLabs2[sppIdx,2]

  file_suffix <- paste(spname,scen,repID2,sep = '-')

  ## where the raw evOsmose outputs are stored
  dirtmp <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/",scen,'/output')

  ## where the WHAM outputs are to be stored
  ## Species-Scenario head folder
  head.dir <- here::here('outputs','wham_runs',paste(spname,scen, sep = '-')); if(!dir.exists(head.dir)) dir.create(head.dir)
  if(is.null(date.use)){
    if(!dir.exists(here::here(head.dir,Sys.Date()))) dir.create(here::here(head.dir,Sys.Date()))
    wham.dir <- here::here(head.dir,Sys.Date(),paste0('rep',repID2)); if(!dir.exists(wham.dir)) dir.create(wham.dir)
  } else{
    wham.dir <- here::here(head.dir,date.use,paste0('rep',repID2)); if(!dir.exists(wham.dir)) dir.create(wham.dir)
  }

  if(file.exists(paste0(wham.dir,"/",file_suffix,'-',fractional_coverage_use, '-wham_survey.csv'))){
    cat(paste('already found outputs for ',file_suffix,"\n"))
    break()
  }

  total_area <- 632 ## total number of marine cells, aka dim (all_cells)

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

  if(fractional_coverage_use==1){
    write.table(laa,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-laa.csv'),
                row.names = FALSE)
  } ## end if fractional_coverage_use == 1

  max_age_pop <- as.numeric(max_age_pop$value)

  #* unfished NAA (for N1_ini) ----
  ## ballpark from earlier runs
  if(fractional_coverage_use==1){
    age_spatial_nofish_path <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/one_sim_without_fishing",
                                      "/ns_spatial_abundancebyAge-",spname,"_Simu0.nc")
    abundance_nofish <- ncvar_get(nc_open(age_spatial_nofish_path),"abundance") ## this might need to switch with units_use
    reshape2::melt(abundance_nofish) %>%
      filter(!is.na(value)) %>% ## drop land
      mutate(year = 2010 + (Var4 - 1) %/% 24,
             month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
      dplyr::select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
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
  } ## end if fractional_coverage_use == 1




  ## fishing mortality x age ----
  ## OSMOSE outputs the matrix of F by size by timestep, though it is invariant thru time, so just take first row
  ## Don't have Frate_age exactly, rather Frate_size so this also needs to convert on a species-basis
  ## it is only saved ONCE under cc-evo
  # f_rate0 <- read.csv(paste0(dirname(dirname(dirtmp)),'/cc_evo',
  #          '/fishing_rate',"/mortality.fishing.rate.byDt.bySize.file.sp",sppIdx-1,".csv"),
  #          nrows = 1, header = T)[,-1] %>%
  #   reshape2::melt() %>%
  #   mutate(length_cm = as.numeric(stringr::str_replace(variable, "X", "")))
  #
  # f_rate0$lenbin <- cut(f_rate0$length_cm, breaks = seq(-0.001, max(f_rate0$length_cm),0.5), right = TRUE)
  #
  # laa$lenbin <- cut(laa$length_cm, breaks = seq(-0.001, max(f_rate0$length_cm),0.5), right = TRUE)
  #
  # f_rate_age <- merge(f_rate0, laa, by = 'lenbin') %>%
  #   select(age, mean_length = length_cm.y, length_bin = length_cm.x, f_rate = value)

  ## mortality (year x age) ----

  ## don't have age-specific values so keep same
  if(fractional_coverage_use==1){
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
  } ## end if fractional_coverage_use == 1
  ## maturity (year x age) ----
  if(fractional_coverage_use==1){
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
      dplyr::select(-year)

    write.table(maturity_data,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-wham_maturity.csv'),
                row.names = FALSE)

    ## putative maturity curve
    maturity_pars <- optim(par = c(log(19),5),
                           fn = maturity_optim_fn,
                           md = maturity_data,
                           amax = max_age_pop)$par
    maturity_curve <- logistic(x = 1:max_age_pop, x0=maturity_pars[1], k = maturity_pars[2])
    write.table(maturity_curve,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-maturity_curve.csv'),
                row.names = FALSE)

  } ## end if fractional_coverage_use == 1
  ## Survey data ----
  ## sample and build index inputs  (year, index as biom/numbers, cv, vector of ages in numbers/biomass, inputN for comps)
  ## need numbers ('abundance') for ages, biomass for indices
  # units_use <- ifelse(units == 'numbers','abundance','biomass') ## how the files are labeled

  #* survey selectivity ----
  survey_selex <<- if(is.na(srv_selex)) {
    cbind(age = 1:max_age_pop, slx = rep(1,max_age_pop))
  } else if(srv_selex == 'mat'){
    maturity_curve <- read.csv(paste0(wham.dir,"/",file_suffix,'-maturity_curve.csv'))
    cbind(age = 1:max_age_pop, slx = maturity_curve$x)
  } else {
    cbind(age = 1:max_age_pop, slx =   1/(1+exp(-log(19)*((1:max_age_pop)-srv_selex)/(max_age_pop-srv_selex))))
  }

  #* filepaths for abundance, biomass ----
  if(fractional_coverage_use==1){
    age_spatial_path <- list.files(dirtmp,
                                   pattern = paste0('spatial_abundancebyAge-',spname),
                                   recursive = T,
                                   full.names = TRUE)[repID]
    biom_spatial_path <- list.files(dirtmp,
                                    pattern = paste0('spatial_biomassbyAge-',spname),
                                    recursive = T,
                                    full.names = TRUE)[repID]

    # repID2 <-  as.numeric(stringr::str_extract(age_spatial_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

    abundance0 <- ncvar_get(nc_open(age_spatial_path),"abundance") ## this might need to switch with units_use
    # this is age (26) lat (25) lon (52) timestep (24/year 2010-2099)
    ## redefine the timesteps into real year, week, and not NA
    ## and then only take out the July populations
    abundance <- reshape2::melt(abundance0) %>%
      filter(!is.na(value) & Var3 <= max_age_pop) %>% ## drop land
      mutate(year = 2010 + (Var4 - 1) %/% 24,
             month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
      dplyr::select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
      filter(month == 7 & year %in% yrs_use) %>% ## July survey
      group_by(year, age, lat, long, month) %>%
      summarise(value = mean(value)) %>% ## average over the month
      ungroup() %>%
      dplyr::select(-month)

    write.csv(abundance,
              paste0(wham.dir,"/",file_suffix,"-abundance.csv"),
              row.names = FALSE)

    biomass0 <- ncvar_get(nc_open(biom_spatial_path),"biomass")
    biomass <- reshape2::melt(biomass0) %>%
      filter(!is.na(value) & Var3 <= max_age_pop) %>% ## drop land
      mutate(year = 2010 + (Var4 - 1) %/% 24,
             month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
      dplyr::select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
      filter(month == 7 & year %in% yrs_use) %>% ## July survey
      group_by(year, age, lat, long, month) %>%
      summarise(value = mean(value)) %>% ## average over the month
      ungroup() %>%
      dplyr::select(-month)

    write.csv(biomass,
              paste0(wham.dir,"/",file_suffix,"-true_biomass_ysa.csv"),
              row.names = FALSE)

    spatial_biomass <- biomass %>%
      summarise(tot_val = sum(value),   .by = c(year, lat, long)) %>%
      # summarise(abundance_rescale_year =  rescale(tot_val, to=c(0,1), na.rm = T),  .by = c(year))
      mutate(abundance_rescale =  rescale(tot_val, to=c(0,1), na.rm = T),
             replicate = repID2,
             scenario = scen,
             species = spname)

    write.csv(spatial_biomass,
              paste0(wham.dir,"/",file_suffix,"-true_biomass_ys.csv"),
              row.names = FALSE)

    #* true biomass, SSB, no space ----
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
              paste0(wham.dir,"/",file_suffix,"-true_biomass_y.csv"),
              row.names = FALSE)
  } else{
    abundance <- read.csv(paste0(wham.dir,"/",file_suffix,"-abundance.csv"))
    biomass <- read.csv(paste0(wham.dir,"/",file_suffix,"-true_biomass_ysa.csv")) ## annual bio
    spatial_biomass <- read.csv(paste0(wham.dir,"/",file_suffix,"-true_biomass_ys.csv")) ## collapsed to year lat long for maps
    true_biomass <- read.csv(paste0(wham.dir,"/",file_suffix,"-true_biomass_y.csv")) ## collapsed to year lat long for maps
  } ## end if fractional_coverage_use != 1
  # define maximum age above which all entries are NA
  max_age_survey <- abundance %>%
    group_by(age) %>%
    summarise(all_zero = all(value == 0)) %>%
    filter(!all_zero) %>%
    summarise(max_age = max(age))
  max_age_survey <- as.numeric(max_age_survey)
  #* run survey ----
  # Initialize an empty list to store the results
  results_age <- results_index <- results_index_gam <- list()
  survey_array <- read.csv(here::here('outputs','wham_runs',
                                      paste0('2024-05-08-survey_array_',
                                             fractional_coverage_use,'.csv')))

  # For each timestep
  for (timestep in unique(abundance$year)) {
    if(timestep %%2!=0) next ## only even years
    cat(timestep,"\n")
    # pull out the survey stations for this year
    selected_cells <- survey_array[survey_array$year == timestep, ]

    ## age comps (numbers)
    # Filter the data for the current timestep
    timestep_data_age <- abundance[abundance$year == timestep, ]
    timestep_data_biom <- biomass[biomass$year == timestep, ]

    ## survey indices
    # MODEL_BASED: generate the spatialized survey data; pass to GAM once all years ready
    survey_biomass_gam <- semi_join(timestep_data_biom, selected_cells, by = c("lat", "long")) %>%
      merge(., survey_selex, by = 'age') %>%
      group_by(lat, long) %>%
      summarise(station_abund = ifelse(is.na(obs_error),
                                       sum(value*slx),
                                       rnorm(1,sum(value*slx),obs_error*mean(value)))) %>% ## sum over all ages
      ungroup() %>%
      mutate(year = timestep,
             replicate = repID2,
             fc = fractional_coverage_use,
             scenario = scenario,
             species = spname)

    results_index_gam[[paste(timestep)]] <- survey_biomass_gam

    # DESIGN-BASED: expand the mean and sd of the abundance for the selected cells
    survey_biomass <- survey_biomass_gam %>%
      summarise(
        abund_mean_per_cell = mean(station_abund),
        abund_mean = abund_mean_per_cell*total_area,
        term1 = sd(station_abund, na.rm = T)*total_area/sqrt(nrow(selected_cells)),
        term2 = sqrt((total_area-nrow(selected_cells))/(total_area-1)),
        abund_se = term1*term2, ## Spencer method with finite pop correction term
        abund_sd = sqrt(var(station_abund, na.rm = T)*total_area^2 / nrow(selected_cells)), ## Oyafuso method
        #abund_cv = abund_sd/abund_mean  ## Oyafuso method
        abund_cv = ifelse(abund_se/abund_mean < 0.025, 0.025, abund_se/abund_mean)
      ) %>%
      dplyr::select(-term1, -term2) %>%
      mutate(year = timestep,
             replicate = repID2,
             fc = fractional_coverage_use,
             scenario = scenario,
             species = spname)

    results_index[[paste(timestep)]] <- survey_biomass

    ## perform multinomial sampling of ages at each station, given selex
    for (i in 1:nrow(selected_cells)) {
      # Filter the data for the current cell (all ages present in population)
      cell_data <- timestep_data_age[timestep_data_age$long == selected_cells$long[i] &
                                       timestep_data_age$lat == selected_cells$lat[i], ]

      if(all(is.na(cell_data$value))) next ## skip if no data
      # cat(i,"\n")
      # Perform the age sampling and store the results
      results_age[[paste(timestep, selected_cells$long[i], selected_cells$lat[i], sep = "_")]] <-
        sample_ages(cell_data, timestep, long = selected_cells$long[i], lat = selected_cells$lat[i],
                    max_age_pop)

    } ## end stations loop for age comps
  } ## end timesteps loop for survey data

  #* combine & save agecomps ----
  # Combine the results into a single data frame
  results_df_age_spatial <- do.call(rbind, results_age) %>%
    filter(age <= max_age_survey)

  # Aggregate the agecomp data by timestep and age, summing the count (collapse space)
  results_df_age <- results_df_age_spatial %>%
    group_by(timestep, age) %>%
    summarise(count = sum(count)) %>%
    dplyr::select(year = timestep, age, count) %>%
    ungroup()

  results_df_age %>%
    mutate( replicate = repID2,
            fc = fractional_coverage_use,
            scenario = scen,
            species = spname) %>%
    write.csv(.,
              paste0(wham.dir,"/",file_suffix,'-',fractional_coverage_use,
                     '-survey_obs_agecomp-numbers.csv'),
              row.names = FALSE)

  #* combine & save indices ----
  #** model-based (gam) ----
  results_df_index_spatial <- do.call(rbind, results_index_gam)

  #** design-based ----
  results_df_index <- do.call(rbind, results_index)
  results_df_index %>%
    mutate(type = 'design-based',
           abund_mean_rescale= rescale(abund_mean, to = c(0,1)),
           replicate = repID2,
           fc = fractional_coverage_use,
           scenario = scen,
           species = spname) %>%
    write.csv(.,
              paste0(wham.dir,"/",file_suffix,'-',fractional_coverage_use,
                     '-survey_obs_biomass.csv'),
              row.names = FALSE)


  # rescale, reshape to WHAM format and save
  # fill missing years with -999
  survey_results <- results_df_index %>%
    mutate(abund_mean = round(abund_mean/units_scalar),
           abund_cv = round(abund_cv,3)) %>%
    dplyr::select(year, abund_mean, abund_cv) %>%
    tidyr::complete(year = 2010:2099,
                    fill = list(abund_mean = -999,
                                abund_sd  = -999,
                                abund_cv  = -999)) %>%
    merge(., results_df_age %>%
            mutate(count = round(count/units_scalar,4)) %>%
            tidyr::complete(year= 2010:2099,
                            age = 1:max_age_pop,
                            fill = list(count = -999) ) %>%
            tidyr::pivot_wider(., names_from = age, values_from = count), by = 'year') %>%
    mutate(inputN = 100) %>%
    arrange(year)

  # save the results
  write.table(survey_results,
              sep = ' ',
              paste0(wham.dir,"/",file_suffix,'-',fractional_coverage_use,
                     '-wham_survey.csv'),
              row.names = FALSE)


  ## Yield data ----
  ## strip and format catches (want yr x Age with totals on end outputs are spp x Age x timestep)
  if(fractional_coverage_use==1){
    yield_files <- list.files(dirtmp, pattern = 'yieldDistribByAge*', recursive = T, full.names = TRUE)
    yield_path <-   list.files(dirtmp,
                               pattern = paste0('ns_yieldDistribByAge*'),
                               recursive = T,
                               full.names = TRUE)[repID]

    repID2 <-  as.numeric(stringr::str_extract(yield_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

    yield0 <- ncvar_get(nc_open(yield_path),"biomass")

    yield1 <- reshape2::melt(yield0) %>%
      mutate(year = 2010 + (Var3 - 1) %/% 24,
             month = ((Var3 - 1) %% 24) %/% 2 + 1,
             age = ifelse(Var2 >= max_age_pop,max_age_pop,Var2)) %>%
      filter(Var1 == sppIdx) %>%
      dplyr::select(year, age, value) %>%
      group_by(year,age) %>%
      summarise(value = sum(value)) %>%
      ungroup()

    # define maximum age above which all entries are NA or zero
    # max_age_catch <- yield1 %>%
    #   group_by(age) %>%
    #   summarise(all_zero = all(value == 0)) %>%
    #   filter(!all_zero) %>%
    #   summarise(max_age = max(age))
    # max_age_catch <- ifelse(as.numeric(max_age_catch)>max_age_pop,max_age_pop,as.numeric(max_age_catch))

    yield1 %>%
      filter(age <= max_age_pop) %>%
      ## truncate age-zeros and max ages
      mutate(value = case_when(age == 1 ~ 0,
                               # age >= max_age_catch  ~ -999,
                               age <= max_age_pop ~ round(value))) %>%
      tidyr::pivot_wider(names_from = age, values_from = value) %>%
      dplyr::select(-year) %>%
      mutate(total = rowSums(.)+999) %>%
      # save the results
      write.table(.,
                  sep = ' ',
                  paste0(wham.dir,"/",file_suffix,'-wham_catch_at_age.csv'),
                  row.names = FALSE)

    ## WAA matrices ----
    #*   population WAA ----
    ## the WAA used to calculate SSB is given by the meanSizeDistribByAge csvs and the allometric w-L parameters in the model
    ## we take this at midyear to match asap3$dat$fracyr_spawn, the fraction of year elapsed before SSB calculation

    waa <- read.csv(header = T, skip = 1,
                    list.files(dirtmp, pattern = 'meanSizeDistribByAge*',
                               recursive = T, full.names = TRUE)[repID]) %>%
      reshape2::melt(id = c('Time','Age')) %>%
      filter(variable == spname & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
      mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
      mutate(rescaled_time =  scales::rescale(Time, to=c(0,1), na.rm = TRUE),   ## start to end year
             compressed_time = cut(rescaled_time, breaks = 26, labels = FALSE), ## ensure same 26 breaks each year
             .by = c(year)) %>%
      mutate(avg_rescaled_time = mean(rescaled_time), .by = compressed_time) %>% ## account for some rounding issues
      filter(avg_rescaled_time > 0.39 & avg_rescaled_time < 0.44) %>% ## filter to spawntime
      group_by(year,Age) %>%
      summarise(mean_size_cm = mean(value)) %>% ## average size-at-age over year
      ungroup() %>%
      mutate(mean_weight_kg  = round(lw_pars$condition*mean_size_cm^lw_pars$allometric/1000,3)) %>% ## average weight at age across year
      mutate(asymp_weight_kg = max(mean_weight_kg),.by = 'year') %>%
      dplyr::select(-mean_size_cm)  %>%
      tidyr::complete(year = 2010:2099, Age = 1:max_age_pop) %>%
      group_by(year) %>%
      tidyr::fill(mean_weight_kg , .direction = "down")  %>%
      ungroup() %>%
      tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = mean_weight_kg) %>%
      dplyr::select(-year)

    write.table(waa,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-wham_waa_ssb.csv'),
                row.names = FALSE)


    #*   catch WAA ----
    #*   Gave this some thought. TL;DR we don't have NAA/NAL in catch
    #*   and therefore can't infer WAA in the catch explicitly.
    #*   Only have biomass-at-age, total numbers, and mean size.
    #*   Avoid using population NAA and monkeying with selex.
    #*   FOR NOW, assume that the catch WAA matches the population


    waa <- read.csv(header = T, skip = 1,
                    list.files(dirtmp, pattern = 'meanSizeDistribByAge*',
                               recursive = T, full.names = TRUE)[repID]) %>%
      reshape2::melt(id = c('Time','Age')) %>%
      filter(variable == spname & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
      mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
      # mutate(rescaled_time =  scales::rescale(Time, to=c(0,1), na.rm = TRUE),   ## start to end year
      #        compressed_time = cut(rescaled_time, breaks = 26, labels = FALSE), ## ensure same 26 breaks each year
      #        .by = c(year)) %>%
      # mutate(avg_rescaled_time = mean(rescaled_time), .by = compressed_time) %>% ## account for some rounding issues
      # filter(avg_rescaled_time > 0.39 & avg_rescaled_time < 0.44) %>% ## filter to spawntime
      group_by(year,Age) %>%
      summarise(mean_size_cm = mean(value)) %>% ## average size-at-age over year
      ungroup() %>%
      mutate(mean_weight_kg  = round(lw_pars$condition*mean_size_cm^lw_pars$allometric/1000,3)) %>% ## average weight at age across year
      mutate(asymp_weight_kg = max(mean_weight_kg),.by = 'year') %>%
      dplyr::select(-mean_size_cm)  %>%
      tidyr::complete(year = 2010:2099, Age = 1:max_age_pop) %>%
      group_by(year) %>%
      tidyr::fill(mean_weight_kg , .direction = "down")  %>%
      ungroup() %>%
      tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = mean_weight_kg) %>%
      dplyr::select(-year)
    write.table(waa,
                sep = ' ',
                paste0(wham.dir,"/",file_suffix,'-wham_waa_catch.csv'),
                row.names = FALSE)

    ## Summary figures and data ----
    #* Catch Figures ----

    ggplot(yield1, aes(x = year, y = value, fill = as.factor(age))) +
      theme(legend.position = 'none')+
      geom_area(alpha =0.4, color = scenLabs2[scenario,'Pal'])+
      labs(x = 'Year', y = 'Yield (kg)', color = 'Year') +
      scale_fill_manual(values =  monochromeR::generate_palette(scenLabs2[scenario,'Pal'],
                                                                modification = "go_lighter",
                                                                n_colours = length(unique(yield1$age)),
                                                                view_palette = FALSE))
    ggsave(last_plot(), file = paste0(wham.dir,"/",file_suffix,"-catch_at_age.png"),
           width = 6, height = 6, unit = 'in', dpi = 400)
  } ## end if fractional_coverage_use == 1
  #* Survey figures ----
  #** maps of true biomass thru time----
  map <-  spatial_biomass %>%
    filter(year %in% floor(seq(2020,max(biomass$year),length.out = 4)))  %>%
    ggplot(data = ., aes(x = lat, y = long,
                         fill = abundance_rescale))+
    theme_void()+
    geom_raster()+
    geom_point(data = filter(survey_array, year %in% floor(seq(2020,max(biomass$year),length.out = 4))),
               fill = NA,
               color = 'red', shape = 4, size = 0.1)+
    # scale_fill_viridis_c(na.value = NA)+
    scale_fill_gradient2(mid = "#FFF5D1",
                         low = "#efeee7",
                         high = scenLabs2$Pal[scenario])+
    theme(strip.text = element_text(size = 25),
          strip.text.y = element_blank(),
          legend.position = 'none')+
    facet_wrap(~year, ncol = 2)

  #** survey index data----
  true_b <- ggplot(true_biomass, aes(x = year, y =abund_mean_rescale)) +
    geom_line(color = scenLabs2[scenario,'Pal'])+
    scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 10))+
    labs(x = 'Year', y = paste0('True Abundance ',
                                ifelse(units == 'numbers','(millions)','(tons)')))

  index <- ggplot(results_df_index %>%
                    mutate(abund_mean_rescale = rescale(abund_mean, to = c(0,1))),
                  aes(x = year, y = abund_mean)) +
    geom_point(color = scenLabs2[scenario,'Pal'])+
    geom_errorbar(aes(ymin = abund_mean - abund_mean*abund_cv,
                      ymax = abund_mean + abund_mean*abund_cv), width = 0, color = scenLabs2[scenario,'Pal']) +
    geom_line(data = true_biomass, color = 'grey70', aes(y = total_biomass))+ ## will only show if in range
    # geom_errorbar(aes(ymin = abund_mean_rescale - abund_cv,
    #ymax = abund_mean_rescale + abund_cv), width = 0, color = scenLabs2[scenario,'Pal']) +
    scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 10))+
    scale_y_continuous(limits = c(0, 1.2*max(results_df_index$abund_mean)), expand = c(0,0))+
    labs(x = 'Year', y = 'Biomass',
         title = paste0('Survey Index, coverage = ', fractional_coverage_use))

  #** survey age comps----
  comps <- results_df_age %>%
    mutate( scount = sum(count), .by = c(year)) %>%
    mutate(frequency = count/scount)  %>%
    ggplot(.) +
    geom_area(alpha =0.4,
              aes(x = age, y = frequency,
                  group = factor(year),
                  fill = factor(year), color = factor(year)))+
    # geom_bar(aes(x = age, y = frequency, group = factor(year), color = factor(year)),stat = 'identity')+
    scale_x_continuous(breaks = seq(0, max(results_df_age$age), by = 2),expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    labs(x = 'Age', y = 'Frequency', color = 'Year', title = 'Survey Age Comps') +
    scale_fill_manual(values =  monochromeR::generate_palette(scenLabs2[scenario,'Pal'],
                                                              modification = "go_lighter",
                                                              n_colours = length(unique(results_df_age$year)),
                                                              view_palette = FALSE)) +
    scale_color_manual(values =  monochromeR::generate_palette(scenLabs2[scenario,'Pal'],
                                                               modification = "go_lighter",
                                                               n_colours = length(unique(results_df_age$year)),
                                                               view_palette = FALSE)) +
    theme(legend.position = 'none', axis.text.y = element_blank()) #+
  # facet_wrap(~year)


  png(file =  paste0(wham.dir,"/",file_suffix, '-',fractional_coverage_use,
                     '-survey_data.png'),
      height = 5, width = 12, unit = 'in',res = 520)
  Rmisc::multiplot(map, index, comps, cols = 3)
  # Rmisc::multiplot(map, index, cols = 2)
  dev.off()


  cat('Built data and summary figures for',
      paste0(spname,
             ' ',scen,' replicate ',
             repID2),"\n")

}
