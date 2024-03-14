#################################################
# Authors:
#   Maia Kapur
#   Megsie Siple
# Description:
# Compare a simple design-based estimator to a model-based one that includes spatiotemporal processes.
#
# From Maia: 'It would be nice to have a simple illustration of whether or not accounting for spatio-temporal processes in a standardization framework dramatically change(s)'
#################################################

raw <- read.csv("data/sampler_raw.csv") # From Maia Feb 6 2024

# lat: 52 vertical latitudinal cells
# lon: 25 horizontal longitudinal cells
# species/species name: integer and string ID for the species. For now I've mainly been inspecting IDs 8, 11, 0, and 3: Atlantic Cod, Whiting, Atl. Herring, European Sprat
# timestep: 1-72, where each number is a bimonthly measure, so three years total. Might want to add a "year" column if you're experimenting with generating an annualized estimate.
# true_biomass: this is actually the sampled/observed biomass, which would get summarized to create an annual value


# Sample randomly from the grid at one timestep ---------------------------

