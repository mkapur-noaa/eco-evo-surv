#################################################
# Authors:
#   Maia Kapur
#   Megsie Siple
# Description:
# Compare a simple design-based estimator to a model-based one that includes spatiotemporal processes.
#
# From Maia: 'It would be nice to have a simple illustration of whether or not accounting for spatio-temporal processes in a standardization framework dramatically change(s)'
#################################################

library(tidyverse)
library(mgcv)

raw <- read.csv("data/sampler_raw.csv") # From Maia Feb 6 2024
# Data description
# lat: 52 vertical latitudinal cells
# lon: 25 horizontal longitudinal cells
# species/species name: integer and string ID for the species. For now I've mainly been inspecting IDs 8, 11, 0, and 3: Atlantic Cod, Whiting, Atl. Herring, European Sprat
# timestep: 1-72, where each number is a bimonthly measure, so three years total. Might want to add a "year" column if you're experimenting with generating an annualized estimate.
# true_biomass: this is actually the sampled/observed biomass, which would get summarized to create an annual value

# Sample randomly from the grid at one timestep ---------------------------
# We don't have data at every timestep for every species
table(raw$species_name, raw$timestep)

# Focus on one species for now
single_sp <- raw |>
  filter(species_name == "European Sprat") # could be anything!
unique(single_sp$timestep)

survstep <- 25 # The timestep when the survey happens

single_sp |>
  filter(timestep == survstep) |>
  ggplot(aes(x = lon, y = lat, fill = log(true_biomass))) +
  geom_tile() +
  coord_equal()

# For a simple random sample, the cells are already subsetted to the 'survey' locations (50% of cells, simple random sample). For now we will assume the survey visits a subset of cells but observes the ecosystem perfectly.

# NOTE: Will need to fix this; we should be using the full grid to predict to. This estimate will only pertain to 50% of the area.
grid <- raw |>
  distinct(lon, lat)

grid |>
  ggplot(aes(x = lon, y = lat)) +
  geom_tile(fill = "red") +
  coord_equal()

get_gam_index <- function(dat, survey_timestep = 12, grid = grid) {
  # dat is the data for a single species - here we use true values from a SRS that assumes every observation is of the true state.
  # survey_timestep is an index between 1 and 72 defining which bimonthly time point has data collected for it.
  # grid is the survey grid (the full area that we want to extrapolate to)
  if (survey_timestep < 1 | survey_timestep > 24) {
    stop("Survey timestep (when the survey happens) must be between 1 and 24, assuming data on true biomass are indexed every two weeks.")
  }

  x <- unique(dat$timestep) / survey_timestep
  y <- x[which(x == floor(x))]
  survidx <- y * survey_timestep # vector of timesteps where you have data
  if (length(survidx) == 0) {
    stop("No years with survey data.")
  }
  survdat <- dat[which(dat$timestep %in% survidx), ]
  survdat$year <- survdat$timestep / survey_timestep

  mod <- gam(
    formula = true_biomass ~ as.factor(year) + # temporal
      s(lon, lat, bs = c("ts")) + # spatial
      s(lon, lat, bs = c("ts"), by = as.factor(year), id = 1), # spatiotemporal
    family = tw(link = "log"),
    data = survdat
  )

  # Replicate grid for each year we have data
  grid_list <- list()

  for (y in 1:length(unique(survdat$year))) {
    df <- grid
    df$year <- unique(survdat$year)[y]
    grid_list[[y]] <- df
  }

  longgrid <- bind_rows(grid_list)
  row.names(longgrid) <- NULL


  pred_gam <- predict(mod,
    type = "response",
    newdata = longgrid
  ) # This takes a long time.
  pred_gam_df <- cbind(longgrid, pred_gam)
  pred_gam_df$area <- 1 # default cell area

  gam_idx_mt <- pred_gam_df |>
    dplyr::group_by(year) |>
    summarize(gam_total_biomass = sum(pred_gam))

  return(gam_idx_mt)
}
