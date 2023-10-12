## Survey simultor for Ev-OSMOSE outputs
## Biomass sampler based on on code from Z Oyafuso (FIMS 2023)
## M Kapur Fall 2023 maia.kapur@noaa.gov
## syntax
## a = age
## x = sex
## m = morph
## s = species (1 = cod, 2 = herring, 3 = dab, 4 = herring)
## t = spatial cell
## z = survey strata (clump of cells)
## y = year

## we ignore selex here

library(here)
library(dplyr)
library(ggplot2)
theme_set(theme_minimal())

## POPULATION SPECIFICATIONS ----
n_a <- 20
n_x <- 2
n_m <- 2
n_s <- 4
n_y <- 20
spp_list <- c('cod','herring','dab','sole')
n_iters <- 25 ## number of stochastic replicates


## load spatial extents
## need barents sea with ID'd lat-long etc
## for mockup use gridded domain
n_t <- 676 ## total number of grid cells in domain
domain <- matrix(NA, nrow = sqrt(n_t), ncol = sqrt(n_t))

grid_cells <- matrix(NA, ncol = 3, nrow = n_t); colnames(grid_cells) <- c('strata','x','y')
idx = 1
for(i in 1:sqrt(n_t)){
for(j in 1:sqrt(n_t)){

grid_cells[idx,'x'] <- i
grid_cells[idx,'y'] <- j
idx = idx+1
}
}
grid_cells[,'strata'] <- 1:n_t



cell_area <- 1 ## km2 per cell
total_area <- cell_area*n_t


## SURVEY SPECIFICATIONS ----
target_n <- 5 ## sampling effort (number of hauls)

## area of each cell. uniform for osmose, but here the assumption is you are chunking the
## total domain into a certain number of strata and want to know the spatial extent
Ah <- tapply(X = cell_area, INDEX = solution, FUN = sum)

## make up a single array slice for one species
## age, sex, morph, strata, year, iteration (replicate)
## we really just need the biomass in each species x cell x year x experiment
biomass_styi <- array(NA, dim = c(n_s,n_t,n_y,n_iters),
dimnames = list(spp_list,NULL,NULL,NULL))

## generate random probabilities for the sepecies
for(iter in 1:n_iters){
for(year in 1:n_y){
    for(ispp in 1:n_s){
## for each of 20 years dump 10k individuals randomly thru the domain
## probs vary by spp, one prob for each cell
        pop_size <- runif(1, 500,10000)
        x <- runif(n_t)
        y <- x/sum(x)
        biomass_styi[ispp,,year,iter] <- rmultinom(n = 1, size = pop_size, prob = y)
}
}
}

 
reshape2::melt(biomass_styi) %>%
dplyr::rename('species' = Var1, 'strata' = Var2, 
'year' = Var3, 'iteration' = Var4, 'biomass' = value) %>% 
merge(.,grid_cells, by = 'strata') -> pop_bio

 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate SRS ----
## Simulate Simple Random Designs at varying sampling efforts   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
index_srs <- cv_srs <- rb_index_srs <- true_index <- array(NA, 
dim = c(n_s, n_y,n_iters),
dimnames = list(spp_list,1:n_y,NULL))

set.seed(234)
for(year in 1:n_y){
  for (iter in 1:n_iters) { ## Loop over iterations -- start

    true_index[,year,iter] <-  rowSums(biomass_styi[,,year,iter]) 
    
    ## Randomly sample grid indices within domain
    srs_sample_idx <- sample(x = 1:n_t, size = target_n)
    ## extract population size of all ages, sexes, morphs at sampled cells (in year 1)
    ## returns a row for each species, columns are locations
    srs_sample <- biomass_styi[,srs_sample_idx,year,iter] 
    
    ## Calculate survey statistics by spp
    srs_tau <- total_area * rowMeans(srs_sample)
    srs_sd <- sqrt(apply(X = srs_sample, MARGIN = 1, FUN = var) *
                     total_area^2 / target_n)
    
    ## Record Survey Performance Metrics
    index_srs[ ,year, iter] <- srs_tau 
    rb_index_srs[ ,year, iter] <- 
      100 * (srs_tau - true_index[,year,iter]) /  true_index[,year,iter]
    cv_srs[, year,iter] <-  srs_sd / srs_tau

## make a master DF with year, species, true, obs, error %, true cv, obs cv

  }  ## Loop over iterations -- end
  } ## loop over years -- end




## PLOTS

## population abundance
ggplot(subset(pop_bio, year == 5), aes(x = x, y = y, fill = biomass)) +
geom_tile() +
facet_wrap(~species)



## compare time series (obs & error)
bind_rows(reshape2::melt(true_index) %>%mutate(src = 'True Index'),
reshape2::melt(index_srs) %>% mutate(src = 'SRS Index')) %>% 
dplyr::rename('species' = Var1, 
'year' = Var2,'biomass' = value) %>%
## add a small offset to year for plotting
mutate(year = ifelse(src == 'True Index', year+5e-1,year)) %>%
group_by(species, year, src) %>%
reframe(enframe(quantile(biomass, c(0.25, 0.5, 0.75)), "quantile", "biomass")) %>% 
tidyr::pivot_wider(id_cols = c(species,year,src),names_from = quantile, values_from = biomass) -> indices

ggplot(indices, aes(x = year, color = src)) +
geom_point(aes(y = `50%`)) +
geom_errorbar(width = 0, aes(ymin = `25%`, ymax = `75%`)) +
scale_color_manual(values = c('#023047','#ffb703')) +
labs(x = 'Year', y = 'Survey Abundance', color = '') +
facet_wrap(~species)

ggsave(last_plot(), file = here('figs',paste0('sampler_mockup_ntows=',target_n,'.png')))

## relative error time series

ggplot(indices, aes(x = year, color = src)) +
geom_point(aes(y = `50%`)) +
geom_errorbar(width = 0, aes(ymin = `25%`, ymax = `75%`)) +
scale_color_manual(values = c('#023047','#ffb703')) +
labs(x = 'Year', y = 'Survey Abundance', color = '') +
facet_wrap(~species)
