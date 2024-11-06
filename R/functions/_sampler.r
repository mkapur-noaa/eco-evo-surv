## sampler functions
## standalone functions to sample individuals from a population
## consider NA handling (precludes a fixed set of cells for all spp/yea)

## establish fixed sampling design for each year
## ensure that the survey is truly multispecies, sampling at the same location for all species in a given year

build_survey_array <- function(n_years = length(2010:2099), fractional_coverage ){
  all_cells <- read.csv(here::here('data','northsea_latlong.csv')) ## only the marine bit
  sample_cells <- function() {
    sample_nn <- nrow(all_cells) * fractional_coverage
    sampled_cells <- sample_n(all_cells, sample_nn)
    return(sampled_cells)
  }
  # Replicate the sampling for each year
  survey_array0 <- replicate(n_years, sample_cells(), simplify = FALSE)
  # Add the year to each sampled data frame and bind them all together
  survey_array0 <- do.call(rbind, lapply(seq_along(survey_array0), function(i){cbind(year = 2010 + i - 1, survey_array0[[i]])}))
  #return(survey_array0)
  write.table(survey_array0,
              sep = ',',
              here::here('outputs','wham_runs',paste0(Sys.Date(),'-survey_array_',
                                                      fractional_coverage,'.csv')),
              row.names = FALSE)

}

# A function that samples a station for abundance (biomass or numbers) and age comps
sample_ages <- function(df = cell_data, timestep=2010, long, lat,max_age_pop) {
  # Exclude NAs
  # df <- df[!is.na(df$value), ]

  # If there are no positive values, return a data frame with zero counts for all age groups
  if (all(df$value <= 0) | nrow(df) == 0){
   return(data.frame(timestep = timestep, long = long, lat = lat,
   age = 1:max_age_pop, count =0))
  }

  ##
  # Calculate the proportions of each age group
  props <- df$value/ sum(df$value)

  # Create a data frame with zero counts for all age groups
  result <- data.frame(timestep = timestep, long = long, lat = lat,
                       age = 1:max_age_pop, count = integer(max(df$age)))

  # Sample 500 individuals using a multinomial distribution & apply selex
  # Fill in the sampled counts for the age groups that were present in the input data
  result[result$age %in% df$age, "count"] <- rmultinom(1, 500, props*survey_selex[,2])

  # Return the number of individuals in each age group
  return(result)
  #return(data.frame(timestep = timestep, long = long, lat = lat, age = df$age, count = result))
}
