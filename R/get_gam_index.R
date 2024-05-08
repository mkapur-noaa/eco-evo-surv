#################################################
# Authors:
#   Maia Kapur
#   Megsie Siple
# Description:
# Compare a simple design-based estimator to a model-based one that includes spatiotemporal processes.
#
# From Maia: 'It would be nice to have a simple illustration of whether or not accounting for spatio-temporal processes in a standardization framework dramatically change(s)'
#################################################

# library(dplyr)
# library(ggplot2)
# library(mgcv)
# library(gratia)
# library(here)
#
# raw <- read.csv(here::here("data","demo_data", "sampler_raw.csv")) # From Maia Feb 6 2024
#
# # Data description
# # lat: 52 vertical latitudinal cells
# # lon: 25 horizontal longitudinal cells
# # species/species name: integer and string ID for the species. For now I've mainly been inspecting IDs 8, 11, 0, and 3: Atlantic Cod, Whiting, Atl. Herring, European Sprat
# # timestep: 1-72, where each number is a bimonthly measure, so three years total. Might want to add a "year" column if you're experimenting with generating an annualized estimate.
# # true_biomass: this is actually the sampled/observed biomass, which would get summarized to create an annual value
#
# # Sample randomly from the grid at one timestep ---------------------------
# # We don't have data at every timestep for every species
# table(raw$species_name, raw$timestep)
#
# # Focus on one species for now
# single_sp <- raw |>
#   filter(species == 3) # could be anything!
# unique(single_sp$timestep)
#
# survstep <- 25 # The timestep when the survey happens
#
# single_sp |>
#   dplyr::filter(timestep == survstep) |>
#   ggplot(aes(x = lon, y = lat, fill = log(true_biomass))) +
#   geom_tile() +
#   coord_equal()
#
# # For a simple random sample, the cells are already subsetted to the 'survey' locations (50% of cells, simple random sample). For now we will assume the survey visits a subset of cells but observes the ecosystem perfectly.
#
# # NOTE: Will need to fix this; we should be using the full grid to predict to. This estimate will only pertain to 50% of the area.
# grid <- raw |>
#   distinct(lon, lat)
#
# grid |>
#   ggplot(aes(x = lon, y = lat)) +
#   geom_tile(fill = "red") +
#   coord_equal()

cv <- function(x) {
  y <- sd(x) / mean(x)
  return(y)
}

## dat is the survey_obs_biomass (what the survey saw)
get_gam_index <- function(survdat = dat,  ## full time series of spatialized survey obs
                          surv_ratio = nrow(survey_array[survey_array$year == 2010, ])/total_area, ## proportional coverage of survey
                          survey_timestep = 7,
                          grid = grid,
                          sims_gratia = 10) {
  # dat is the data for a single species - here we use true values from a SRS that assumes every observation is of the true state.
  # survey_timestep is an index between 1 and 72 defining which bimonthly time point has data collected for it.
  # grid is the survey grid (the full area that we want to extrapolate to)
  # if (survey_timestep < 1 | survey_timestep > 24) {
  #   stop("Survey timestep (when the survey happens) must be between 1 and 24, assuming data on true biomass are indexed every two weeks.")
  # }

  # nyr <- max(dat$timestep) / 24
  # survsteps <- survey_timestep + seq(0, (nyr - 1) * 24, by = 24) # which timesteps have a survey in them
  #
  # if (length(survsteps) == 0) {
  #   stop("No years with survey data.")
  # }
  # survdat <- dat[which(dat$timestep %in% survsteps), ]
  #
  # # survdat$year <- NA
  # for (i in 1:nrow(survdat)) {
  #   survdat$year[i] <- which(survsteps == survdat$timestep[i])
  # }

  mod <- gam(
    formula = station_abund ~ as.factor(year) + # temporal
      s(long, lat, bs = c("ts")) + # spatial
      s(long, lat, bs = c("ts"), by = as.factor(year), id = 1), # spatiotemporal
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

  longgrid <- dplyr::bind_rows(grid_list)
  row.names(longgrid) <- NULL

  # use gratia to draw fitted values from the posterior of the model
  sims <- gratia::fitted_samples(mod,
                                 n = sims_gratia, data = longgrid,
                                 scale = "response", seed = 123
  )
  sims$year <- longgrid$year[sims$row]
  sims$area <- 1
  sims$biomass <- sims$fitted * sims$area # expand from density to biomass, given area - this is not needed when all areas = 1 but I'm leaving it in for thoroughness

  gam_idx_mt <- sims |>
    dplyr::group_by(year, draw) |>
    dplyr::summarise_at("biomass", list(biomass = sum)) |>
    dplyr::group_by(year) |>
    dplyr::summarize_at("biomass", list(est = mean, sd = sd, CV = cv)) |>
    dplyr::ungroup()

  return(list(idx = gam_idx_mt, survey_data = survdat))
}

sprat_idx <- get_gam_index(dat = single_sp, survey_timestep = 12, grid = grid, sims_gratia = 1000)

# NOTE: You'll want to multiply these indices by two to get the total estimated index for a given year, because the predictions are only for the 50% of the area that was surveyed (if that makes sense). Basically, since I don't have access to the full grid for the area that the survey is representing, I am only predicting the total biomass of the surveyed area, here.

db <- sprat_idx$survey_data |>
  group_by(year) |>
  summarize(idx = 1/surv_ratio*sum(true_biomass), #These say "true" but don't be fooled! They're only calculated from a subset of the data!
            sd = sd(true_biomass),
            cv = cv(true_biomass))

mb <- sprat_idx$idx

plot(idx/1000 ~ year, data = db,
     xlab = "year", ylab = "index", pch = 19,
     ylim = c(min(db$idx/1000),
              max(db$idx/1000)*1.5))
lines(idx/1000 ~ year, data = db)
segments(x0 = db$year,
         x1 = db$year,
         y0 = db$idx/1000 - db$sd/1000,
         y1 = db$idx/1000 + db$sd/1000)

points(est/1000 ~ year, data = mb, col = "red", pch = 19)
lines(est/1000 ~ year, data = mb, col = "red")
segments(x0 = db$year,
         x1 = db$year,
         y0 = mb$est/1000 - mb$sd/1000,
         y1 = mb$est/1000 + mb$sd/1000,
         col = "red")

legend("topright", legend = c("design-based", "model-based"), col = c("black", "red"), pch = c(19, 19), lty = c(1, 1))
