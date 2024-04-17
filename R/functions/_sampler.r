## sampler functions
## standalone functions to sample individuals from a population
## consider NA handling (precludes a fixed set of cells for all spp/yea)

## establish fixed sampling design for each year
## ensure that the survey is truly multispecies, sampling at the same location for all species in a given year

build_survey_array <- function(n_years = length(2010:2099), fractional_coverage = 0.5){
  #Create a data frame with all combinations of lat and long
  all_cells <- expand.grid(lat = 1:25, long = 1:52)
  # Function to sample 50% of the cells
  sample_cells <- function() {
    sample_n <- nrow(all_cells) * fractional_coverage
    sampled_cells <- sample_n(all_cells, sample_n)
    #return(sampled_cells)
  }
  # Replicate the sampling for each year
  survey_array0 <- replicate(n_years, sample_cells(), simplify = FALSE)
  # Add the year to each sampled data frame and bind them all together
  survey_array0 <- do.call(rbind, lapply(seq_along(survey_array0), function(i){cbind(year = 2010 + i - 1, survey_array0[[i]])}))
  #return(survey_array0)
}

# A function that samples a station for abundance (biomass or numbers) and age comps
sample_ages <- function(df, timestep, long, lat) {
  # Exclude NAs
  df <- df[!is.na(df$value), ]

  # If there are no positive values, return a data frame with zero counts for all age groups
  if (all(df$value <= 0) | nrow(df) == 0){
   return(data.frame(timestep = timestep, long = long, lat = lat,
   age = 1:26, count =0))
  }
  # Calculate the proportions of each age group
  props <- df$value / sum(df$value)

  # Sample 500 individuals using a multinomial distribution
  sampled <- rmultinom(1, 500, props)

  # Create a data frame with zero counts for all age groups
  result <- data.frame(timestep = timestep, long = long, lat = lat,
                       age = 1:max(df$age), count = integer(max(df$age)))

  # Fill in the sampled counts for the age groups that were present in the input data
  result[result$age %in% df$age, "count"] <- sampled

  # Return the number of individuals in each age group
  return(result)
  #return(data.frame(timestep = timestep, long = long, lat = lat, age = df$age, count = result))
}

