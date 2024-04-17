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
                     repuse = 1,
                     yrs_use = 2010:2099,
                     srv_selex = 8, ## whether or not there is survey selex (logistic with A50=8)
                     obs_error = 0.2, ## whether or not survey biomass has observation error
                     units = 'numbers',
                     units_scalar = 1e6){

  ## where the raw evOsmose outputs are stored
  dirtmp <- paste0("F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/",scenLabs2[scenario,2],'/output')

  ## where the WHAM outputs are to be stored
  wham.dir <- here::here('outputs','wham_runs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],"-",Sys.Date()))
  if(!dir.exists(wham.dir)) dir.create(wham.dir)
  # setwd(wham.dir)

  ## mortality (year x age) ----
  ## don't have age-specific values so keep same
  readLines(mort_path,n=2,skip = 1)
  mort_path <- paste0(dirtmp, '/mortality/', 'ns_mortalityRate-',sppLabs2[sppIdx,2],"_Simu",repuse,".csv")
  mort_csv <- read.csv(mort_path, skip = 3, header = F)[,c(1,12,18)] %>% ## timestep, Frecruits, Mrecruits
    mutate(year = floor(as.numeric(stringr::str_split_fixed(V1, "\\.", 1)) -70+2010)) %>%
    group_by(year) %>%
    summarise(Frecruits = round(sum(V12),4), Mrecruits = round(sum(V18),4))

  matrix(rep(mort_csv$Mrecruits, length(1:26)), byrow=FALSE, ncol = length(1:26)) %>%
    write.table(.,
                sep = ' ',
                paste0(wham.dir,"/",sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_mortality.csv'),
                row.names = FALSE)

  ## maturity (year x age) ----
  mat_path <- paste0(dirtmp, '/ageIndicators/', 'ns_maturityDistribByAge',"_Simu",repuse,".csv")
  read.csv(mat_path, skip = 1, header = T) %>%
    reshape2::melt(id = c('Time','Age')) %>%
    filter(variable == sppLabs2[sppIdx,2] & !is.na(value) & Age > 0)  %>% ## get species- and age-specific values
    mutate(year = floor(as.numeric(stringr::str_split_fixed(Time, "\\.", 1)) -70+2010)) %>%
    group_by(year,Age) %>%
    summarise(value = mean(value)) %>% ## average maturity across year
    tidyr::complete(Age = 1:26,
                    fill = list(value = 1)) %>%
    tidyr::pivot_wider(., id_cols = year, names_from = Age, values_from = value) %>%
    write.table(.,
                sep = ' ',
                here::here('output','wham_runs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_maturity.csv')),
                row.names = FALSE)


  ## WAA matrices (catch, discards (unused), ssb)


  ## strip and format catches (yr x Age)
  yield_files <- list.files(dirtmp, pattern = 'ns_yield*', recursive = T, full.names = TRUE)

  ## Survey data ----
  ## sample and build index inputs  (year, index as biom/numbers, cv, vector of ages in numbers/biomass, inputN for comps)
  ## need numbers ('abundance') for ages, biomass for indices
  # units_use <- ifelse(units == 'numbers','abundance','biomass') ## how the files are labeled

  survey_selex <- if(is.null(srv_selex)) {rep(1,26)} else { 1/(1+exp(-log(19)*((1:26)-srv_selex)/(15-srv_selex)))}
  survey_selex <- cbind(age = 1:26, slx = survey_selex)

  age_spatial_path <- list.files(dirtmp,
                                 pattern = paste0('spatial_abundancebyAge-',sppLabs2[sppIdx,2]),
                                 recursive = T,
                                 full.names = TRUE)[repuse]
  biom_spatial_path <- list.files(dirtmp,
                                 pattern = paste0('spatial_biomassbyAge-',sppLabs2[sppIdx,2]),
                                 recursive = T,
                                 full.names = TRUE)[repuse]

  repID <-  as.numeric(stringr::str_extract(age_spatial_path, "(?<=Simu)\\d+(?=\\.nc)")) ## might not match replicate input

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
  max_age <- abundance %>%
    group_by(age) %>%
    #summarise(mean(value, na.RM = TRUE)) %>% tail(10)
    summarise(all_zero = all(value == 0)) %>%
    filter(!is.na(all_zero)) %>%
    summarise(max_age = max(age))
  max_age <- as.numeric(max_age)

  total_area <- length(unique(abundance$lat))*length(unique(abundance$long))

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
        abund_sd = sqrt(var(station_abund, na.rm = T)*total_area^2 / nrow(selected_cells)),
        abund_cv = abund_sd/abund_mean ) %>%
      mutate(year = timestep, replicate = repID, scenario = scenario, species = sppLabs2[sppIdx,2])

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
    filter(age <= max_age)
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
            here::here('output','wham_runs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',scenLabs2[scenario,2],'-wham_survey.csv')),
            row.names = FALSE)


  ## summary figures
  ## maps of true biomass thru time
  map <- biomass %>%
    group_by(year, lat, long) %>%
    summarise(tot_val = sum(value)) %>%
    mutate(abundance_rescale =  rescale(tot_val, to=c(0,1), na.rm = T)) %>%
    ungroup() %>%
    filter(year %in% floor(seq(2020,max(biomass$year),length.out = 4)))  %>%
    ggplot(data = ., aes(x = lat, y = long,
                         fill = abundance_rescale))+
    theme_void()+
    geom_raster()+
    geom_point(data = filter(survey_array, year %in% floor(seq(2020,max(biomass$year),length.out = 4))),
          fill = NA,
               color = 'red', shape = 4)+
    scale_fill_viridis_c(na.value = NA)+
    theme(strip.text = element_text(size = 25),
          strip.text.y = element_blank(),
          legend.position = 'none')+
    facet_wrap(~year, ncol = 2)

  ## survey index data
  true_abund <- biomass %>%
    group_by(year) %>%
    summarise(abund_mean=sum(value,na.rm = T)) %>%
    ungroup() %>%
    mutate(abund_mean_rescale= rescale(abund_mean, to = c(0,1)))

  true_b <- ggplot(true_abund, aes(x = year, y =abund_mean_rescale)) +
    geom_line(color = scenLabs2[scenario,'Pal'])+
    scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 10))+
    labs(x = 'Year', y = paste0('True Abundance ',
                                ifelse(units == 'numbers','(millions)','(tons)')))


  index <- ggplot(results_df_index %>%
                    mutate(abund_mean_rescale = rescale(abund_mean, to = c(0,1))),
                  aes(x = year, y = abund_mean_rescale)) +
    geom_point(color = scenLabs2[scenario,'Pal'])+
    geom_line(data = true_abund, color = 'grey3')+
    geom_errorbar(aes(ymin = abund_mean_rescale - abund_cv,
                      ymax = abund_mean_rescale + abund_cv), width = 0, color = scenLabs2[scenario,'Pal']) +
    scale_x_continuous(breaks = seq(min(yrs_use), max(yrs_use), by = 10))+
    labs(x = 'Year', y = 'Biomass (rescaled)', title = 'Survey Index')

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
    scale_x_continuous(breaks = seq(0, max(results_df_age$age), by = 2))+
    labs(x = 'Age', y = 'Frequency', color = 'Year', title = 'Surveyed Age Comps') +
    scale_fill_manual(values =  monochromeR::generate_palette(scenLabs2[scenario,'Pal'],
                                                               modification = "go_lighter",
                                                               n_colours = length(unique(results_df_age$year)),
                                                               view_palette = FALSE)) +
    scale_color_manual(values =  monochromeR::generate_palette(scenLabs2[scenario,'Pal'],
                                                               modification = "go_lighter",
                                                               n_colours = length(unique(results_df_age$year)),
                                                               view_palette = FALSE)) +
    theme(legend.position = 'none') #+
    # facet_grid(year~., scales = 'free') +
    # theme(panel.spacing = unit(0.2, "lines"),)



  png(file = here::here('figs',paste0(sppLabs2[sppIdx,2],'-rep',repID,'-',
                                      scenLabs2[scenario,2],'-abundance-',units,'-',Sys.Date(),'.png')),
      height = 8, width = 12, unit = 'in',res = 520)
  Rmisc::multiplot(map, index, comps, cols = 3)
  # Rmisc::multiplot(map, index, cols = 2)

  dev.off()


  cat('Built data and summary figures for',paste0(sppLabs2[sppIdx,2],
                                                  ' ',scenLabs2[scenario,2],' replicate ',
                                                  repID))

}
