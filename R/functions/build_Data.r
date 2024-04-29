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
                     yrs_use = 2010:2099,
                     srv_selex = 11, ## whether or not survey sampled with age-based selex (logistic with A50=11)
                     obs_error = 0.2, ## whether or not survey sampled with observation error
                     units = 'biomass', ## units for survey observations
                     units_scalar = 1){

  ## where the raw evOsmose outputs are stored
  dirtmp <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/",scenLabs2[scenario,2],'/output')

  ## where the WHAM outputs are to be stored
  wham.dir <- here::here('outputs','wham_runs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],"-",Sys.Date()))
  if(!dir.exists(wham.dir)) dir.create(wham.dir)
  # setwd(wham.dir)

  ## load parameters for this species ----
  lw_pars <- read.csv(here::here('outputs','wham_runs','length2weight.csv')) %>%  filter(species == sppLabs2[sppIdx,2])
  max_age_pop <- read.csv(here::here('outputs','wham_runs','max_age.csv')) %>%  filter(species == sppLabs2[sppIdx,2])
  vonBpars <- read.csv(here::here('outputs','wham_runs','vonBpars.csv')) %>%  filter(species == sppLabs2[sppIdx,2])

  #* populate length-at-age vector
  laa <- data.frame(length_cm = NA, age = NA, len_bin = NA)
  for(a in 1:max_age_pop$value){
    laa[a,'age'] <- a
    laa[a,'length_cm'] <- vonBpars$lInf*(1-exp(-vonBpars$K*(a-vonBpars$t0)))
  }


  ## fishing mortality x age ----
  ## OSMOSE outputs the matrix of F by size by timestep, though it is invariant thru time, so just take first row
  ## Don't have Frate_age exactly, rather Frate_size so this also needs to convert on a species-basis
  ## it is only saved ONCE under cc-evo
  f_rate0 <- read.csv(paste0(dirname(dirname(dirtmp)),'/cc_evo',
           '/fishing_rate',"/mortality.fishing.rate.byDt.bySize.file.sp",sppIdx-1,".csv"),
           nrows = 1, header = T)[,-1] %>%
    reshape2::melt() %>%
    mutate(length_cm = as.numeric(stringr::str_replace(variable, "X", "")))

  f_rate0$lenbin <- cut(f_rate0$length_cm, breaks = seq(-0.001, max(f_rate0$length_cm),0.5), right = TRUE)

  laa$lenbin <- cut(laa$length_cm, breaks = seq(-0.001, max(f_rate0$length_cm),0.5), right = TRUE)

  f_rate_age <- merge(f_rate0, laa, by = 'lenbin') %>%
    select(age, mean_length = length_cm.y, length_bin = length_cm.x, f_rate = value)

  # Merge df1 and df2 based on LENGTH_bin
  merged_df <- merge(df1, df2, by.x = "LENGTH", by.y = "LENGTH_bin")


  ## mortality (year x age) ----

  ## don't have age-specific values so keep same

  mort_path <- paste0(dirtmp, '/mortality/', 'ns_mortalityRate-',sppLabs2[sppIdx,2],"_Simu",repID,".csv")
  mort_csv <- read.csv(mort_path, skip = 3, header = F)[,c(1,12,18)] %>% ## timestep, Frecruits, Mrecruits
    mutate(year = floor(as.numeric(stringr::str_split_fixed(V1, "\\.", 1)) -70+2010)) %>%
    group_by(year) %>%
    summarise(Frecruits = round(sum(V12),4), Mrecruits = round(sum(V18),4)) %>% ungroup()

  matrix(rep(mort_csv$Mrecruits, length(1:26)), byrow=FALSE, ncol = length(1:26)) %>%
    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_mortality.csv'),
                row.names = FALSE)

  ## maturity (year x age) ----

  mat_path <- paste0(dirtmp, '/ageIndicators/', 'ns_maturityDistribByAge',"_Simu",repID,".csv")
  read.csv(mat_path, skip = 1, header = T) %>%
    reshape2::melt(id = c('Time','Age')) %>%
    filter(variable == sppLabs2[sppIdx,2] & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
    mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
    group_by(year,Age) %>%
    summarise(value = round(mean(value),3)) %>% ## average maturity across year
    ungroup() %>%
    tidyr::complete(year = 2010:2099, Age = 1:26,
                    fill = list(value = 1)) %>%
    tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = value) %>%
    select(-year) %>%
    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_maturity.csv'),
                row.names = FALSE)








  ## Survey data ----
  ## sample and build index inputs  (year, index as biom/numbers, cv, vector of ages in numbers/biomass, inputN for comps)
  ## need numbers ('abundance') for ages, biomass for indices
  # units_use <- ifelse(units == 'numbers','abundance','biomass') ## how the files are labeled

  survey_selex <- if(is.null(srv_selex)) {rep(1,26)} else { 1/(1+exp(-log(19)*((1:26)-srv_selex)/(15-srv_selex)))}
  survey_selex <- cbind(age = 1:26, slx = survey_selex)

  age_spatial_path <- list.files(dirtmp,
                                 pattern = paste0('spatial_abundancebyAge-',sppLabs2[sppIdx,2]),
                                 recursive = T,
                                 full.names = TRUE)[repID]
  biom_spatial_path <- list.files(dirtmp,
                                  pattern = paste0('spatial_biomassbyAge-',sppLabs2[sppIdx,2]),
                                  recursive = T,
                                  full.names = TRUE)[repID]

  repID2 <-  as.numeric(stringr::str_extract(age_spatial_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

  abundance0 <- ncvar_get(nc_open(age_spatial_path),"abundance") ## this might need to switch with units_use
  ## this is age (26) lat (25) lon (52) timestep (24/year 2010-2099)
  ## redefine the timesteps into real year, week, and not NA
  ## and then only take out the July populations
  abundance <- reshape2::melt(abundance0) %>%
    filter(!is.na(value)) %>% ## drop land
    mutate(year = 2010 + (Var4 - 1) %/% 24,
           month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7 & year %in% yrs_use) %>%
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
  max_age_survey <- as.numeric(max_age)

  total_area <- length(unique(abundance$lat))*length(unique(abundance$long)) ## total number of cells

  biomass0 <- ncvar_get(nc_open(biom_spatial_path),"biomass") ## this might need to switch with units_use
  biomass <- reshape2::melt(biomass0) %>%
    filter(!is.na(value)) %>% ## drop land
    mutate(year = 2010 + (Var4 - 1) %/% 24,
           month = ((Var4 - 1) %% 24) %/% 2 + 1) %>%
    select(lat=Var1, long=Var2, age=Var3, year, value, month) %>%
    filter(month == 7 & year %in% yrs_use) %>%
    group_by(year, age, lat, long, month) %>%
    summarise(value = mean(value)) %>% ## average over the month
    ungroup() %>%
    select(-month)


  # Initialize an empty list to store the results
  results_age <- results_index <- list()

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
    # expand the mean and sd of the abundance for the selected cells
    survey_biomass <- semi_join(timestep_data_biom, selected_cells, by = c("lat", "long")) %>%
      filter(!is.na(value)) %>% ## remove NAs
      merge(., survey_selex, by = 'age') %>%
      group_by(lat, long) %>%
      summarise(station_abund = ifelse(is.null(obs_error),
                                       sum(value*slx),
                                       rnorm(1,sum(value*slx),obs_error*mean(value)))) %>% ## sum over all ages
      ungroup() %>%
      summarise(
        abund_mean = mean(station_abund)*total_area,
        term1 = sd(station_abund, na.rm = T)*total_area/sqrt(nrow(selected_cells)),
        term2 = sqrt((total_area-nrow(selected_cells))/(total_area-1)),
        abund_se = term1*term2, ## Spencer method with finite pop correction term
        abund_sd = sqrt(var(station_abund, na.rm = T)*total_area^2 / nrow(selected_cells)), ## Oyafuso method
        #abund_cv = abund_sd/abund_mean  ## Oyafuso method
        abund_cv = abund_se/abund_mean
      ) %>%
      select(-term1, -term2) %>%
      mutate(year = timestep, replicate = repID2, scenario = scenario, species = sppLabs2[sppIdx,2])

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
        sample_ages(cell_data, timestep, long = selected_cells$long[i], lat = selected_cells$lat[i])

    } ## end stations loop for age comps
  } ## end timesteps loop for survey data

  # Combine the results into a single data frame
  results_df_age_spatial <- do.call(rbind, results_age) %>%
    filter(age <= max_age_survey)
  # mutate(count = ifelse(age <= max_age, count, -999)) ## blank data for unused ages

  # Aggregate the agecomp data by timestep and age, summing the count (collapse space)
  results_df_age <- results_df_age_spatial %>%
    group_by(timestep, age) %>%
    summarise(count = sum(count)) %>%
    select(year = timestep, age, count) %>%
    ungroup()

  results_df_index <- do.call(rbind, results_index)

  # rescale, reshape to WHAM format and save
  # fill missing years with -999
  survey_results <- results_df_index %>%
    mutate(abund_mean = round(abund_mean/units_scalar),
           abund_cv = round(abund_cv,3)) %>%
    select(year, abund_mean, abund_cv) %>%
    tidyr::complete(year = 2010:2099,
                    fill = list(abund_mean = -999,
                                abund_sd  = -999,
                                abund_cv  = -999)) %>%
    merge(., results_df_age %>%
            mutate(count = round(count/units_scalar,4)) %>%
            tidyr::complete(year= 2010:2099,
                            age = 1:26,
                            fill = list(count = -999) ) %>%
            tidyr::pivot_wider(., names_from = age, values_from = count), by = 'year') %>%
    mutate(inputN = 50) %>%
    arrange(year)

  # save the results
  write.table(survey_results,
              sep = ' ',
              paste0(wham.dir,"/",sppLabs2[sppIdx,2],
                     '-rep',repID2,'-',Sys.Date(),
                     '-wham_survey.csv'),
              row.names = FALSE)


  ## Yield data ----
  ## strip and format catches (want yr x Age, outputs are spp x Age x timestep)
  yield_files <- list.files(dirtmp, pattern = 'yieldDistribByAge*', recursive = T, full.names = TRUE)
  yield_path <-   list.files(dirtmp,
                             pattern = paste0('ns_yieldDistribByAge*'),
                             recursive = T,
                             full.names = TRUE)[repID]

  repID2 <-  as.numeric(stringr::str_extract(yield_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

  yield0 <- ncvar_get(nc_open(yield_path),"biomass")

  yield1 <- reshape2::melt(yield0) %>%
    mutate(year = 2010 + (Var3 - 1) %/% 24,
           month = ((Var3 - 1) %% 24) %/% 2 + 1) %>%
    filter(Var1 == sppIdx) %>%
    select(year, age = Var2, value) %>%
    group_by(year,age) %>%
    summarise(value = sum(value)) %>%
    ungroup()

  # define maximum age above which all entries are NA or zero
  # this might need to be 26 or the population max age, no matter what
  max_age_catch <- yield1 %>%
    group_by(age) %>%
    summarise(all_zero = all(value == 0)) %>%
    filter(!all_zero) %>%
    summarise(max_age = max(age))
  max_age_catch <- as.numeric(max_age_catch)

  yield1 %>%
    ## truncate age-zeros and max ages
    mutate(value = case_when(age == 1 ~ -999,
                             age >= max_age_catch  ~ -999,
                             age < max_age_catch ~ round(value))) %>%
    tidyr::pivot_wider(names_from = age, values_from = value) %>%
    select(-year) %>%

    # save the results
    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",sppLabs2[sppIdx,2],
                       '-rep',repID2,'-',Sys.Date(),
                       '-wham_catch_at_age.csv'),
                row.names = FALSE)

  ## WAA matrices ----


  #*   population WAA ----
  ## the WAA used to calculate SSB is given by the meanSizeDistribByAge csvs and the allometric w-L parameters in the model
  read.csv(header = T, skip = 1,
           list.files(dirtmp, pattern = 'meanSizeDistribByAge*',
                      recursive = T, full.names = TRUE)[repID]) %>%
    reshape2::melt(id = c('Time','Age')) %>%
    filter(variable == sppLabs2[sppIdx,2] & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
    mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
    group_by(year,Age) %>%
    summarise(mean_size_cm = mean(value)) %>% ## average size-at-age over year
    ungroup() %>%
    mutate(  mean_weight_kg  = round(lw_pars$condition*mean_size_cm^lw_pars$allometric/1000,3)) %>% ## average weight at age across year
    mutate(asymp_weight_kg = max(mean_weight_kg),.by = 'year') %>%
    select(-mean_size_cm) %>%
    tidyr::complete(year = 2010:2099, Age = 1:26,
                    fill = list(value = .$asymp_weight_kg)) %>%
    tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = mean_weight_kg) %>%
    select(-year) %>% View()
  write.table(.,
              sep = ' ',
              paste0(wham.dir,"/",sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],
                     '-wham_waa_ssb.csv'),
              row.names = FALSE)



  #*   catch WAA ----
  #*   this needs to be interpolated because we don't have a measure of NAA in catch,
  #*   only biomass at age, total numbers, and mean size
  #*   estimate NAA in catch via Frate_age * NAA of the population, then divide into biomass-at-age
  #*   to approximate the weight-at-age of a single landed fish in a given year
  #*
  #*   Sanity check that the derived WAA seems to return the landed biomass...
  true_abundance <- abundance %>%
    group_by(year,age) %>%
    summarise(pop_naa = sum(value,na.rm = T))


  ## Summary figures ----
  #* survey figures ----
  ## maps of true biomass thru time
  map <- biomass %>%
    #group_by(lat, long) %>%
    #group_by(year, lat, long) %>%
    summarise(tot_val = sum(value), .by = c(year, lat, long)) %>%
    ungroup() %>%
    mutate(abundance_rescale =  rescale(tot_val, to=c(0,1), na.rm = T)) %>%

    filter(year %in% floor(seq(2020,max(biomass$year),length.out = 4)))  %>%
    ggplot(data = ., aes(x = lat, y = long,
                         fill = abundance_rescale))+
    theme_void()+
    geom_raster()+
    geom_point(data = filter(survey_array, year %in% floor(seq(2020,max(biomass$year),length.out = 4))),
               fill = NA,
               color = 'red', shape = 4, size = 0.1)+
    scale_fill_viridis_c(na.value = NA)+
    theme(strip.text = element_text(size = 25),
          strip.text.y = element_blank(),
          legend.position = 'none')+
    facet_wrap(~year, ncol = 2)

  ## survey index data
  true_biomass <- biomass %>%
    group_by(year) %>%
    summarise(abund_mean=sum(value,na.rm = T)) %>%
    ungroup() %>%
    mutate(abund_mean_rescale= rescale(abund_mean, to = c(0,1)))

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
    #geom_line(data = true_abund, color = 'grey3')+
    #geom_errorbar(aes(ymin = abund_mean_rescale - abund_cv,
    #ymax = abund_mean_rescale + abund_cv), width = 0, color = scenLabs2[scenario,'Pal']) +
    scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 10))+
    scale_y_continuous(limits = c(0, 1.2*max(results_df_index$abund_mean)), expand = c(0,0))+
    labs(x = 'Year', y = 'Biomass', title = 'Survey Index')

  ## survey age comps
  comps <- results_df_age %>%
    group_by(year) %>%
    mutate(frequency = count / sum(count)) %>%
    ungroup() %>%
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
  # facet_grid(year~., scales = 'free') +
  # theme(panel.spacing = unit(0.2, "lines"),)



  png(file = paste0(wham.dir,"/",sppLabs2[sppIdx,2],
                    '-rep',repID2,'-',Sys.Date(),
                    '-input_data.png'),
      height = 5, width = 12, unit = 'in',res = 520)
  Rmisc::multiplot(map, index, comps, cols = 3)
  # Rmisc::multiplot(map, index, cols = 2)
  dev.off()


  cat('Built data and summary figures for',paste0(sppLabs2[sppIdx,2],
                                                  ' ',scenLabs2[scenario,2],' replicate ',
                                                  repID2))

}
