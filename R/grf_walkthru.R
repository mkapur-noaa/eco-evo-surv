## Illustrative example using gaussian random fields
## of how spatial heterogeneity can induce bias in samples.

library(geoR)
library(here)
library(dplyr)
library(ggplot2)
library(stringr)
set.seed(731)

## BUILD SIMULATED POPULATIONS ----
## four plots along spatial-morpho axis
## cov pars are sigma-sq (distance) and phi
## the 'reg' grid has extent 0, 1

pars_sims <- matrix(c(0,0,0.5,0,5,0,5,1),ncol = 2, byrow = TRUE) ## sigma, phi by scenario
border_mat <- matrix(c(5,0,0,5,0,0,5,5),ncol = 2, byrow = FALSE) ## domain


## initial examples (for viz)
for(isim in 1:4){
tmp <- grf(1e3, grid = 'reg', 
cov.pars = pars_sims[isim,],
borders = border_mat)
assign(value = tmp, x =paste0('sim',isim))
}

## create df of results
lapply(list(sim1,sim2,sim3,sim4), 
FUN = function(x){
data.frame(x$coords,x$data) %>%
mutate(X1 = round(X1,2), X2 = round(X2,2), cellid = seq(length(sim1$data)))
}) %>% 
bind_rows %>%
mutate(id = paste0('sim',rep(1:4,each = length(sim1$data)))) -> sims
sims$X1 <- factor(sims$X1)
sims$X2 <- factor(sims$X2)
#save(sims, file = here('data','grf_sims.rdata'))



## ITERATE SAMPLES ----
## Generate n_rep realizations of the spatial/morpho phenomenae parameterized above.
grf_fun <- function(isim){
x <- grf(1e3, grid = 'reg', 
cov.pars = pars_sims[isim,],
borders = border_mat)
return(cbind(x$coords,x$data))
}

big_sims <- NULL
n_reps = 100

for(isim in 1:4){
tmp <- do.call( rbind, replicate(n_reps, grf_fun(isim), simplify=FALSE ) ) %>%
data.frame() %>%
mutate(X1 = factor(round(X1,2)), X2 = factor(round(X2,2)), 
simid = paste0('sim',isim),
cellid = rep(seq(length(sim1$data)),n_reps), 
repid = rep(1:n_reps, each = length(sim1$data))) %>%
select(X1,X2,x.data = X3,cellid,simid, repid)
big_sims <- rbind(big_sims, tmp); rm(tmp)
}

#save(big_sims, file = here('data','grf_big_sims.rdata'))

## SAMPLING ROUTINE ----
library(scales)
## For a given tow number, sample the 'biomass' in the cells and return 
## the true and observed mean & cv

## show sampling routine
n_stations <- 25
stations <- sample(seq(length(sim1$data)),n_stations) ## 100 random sampling sites within domain (don't change)
total_area <- length(sim1$data) ## the number of cells; assume each cell = 1

big_sampler <- matrix(NA, nrow = 4*n_reps, ncol = 10) ## replicates x sims x variable
idx = 1
for(isim in 1:4){ ## loop four scenarios
for(irep in 1:n_reps){ ## loop n_rep realizations

sim_true <- subset(big_sims, simid == paste0('sim',isim) & repid == irep)['x.data']
sim_true <- rescale(sim_true$x.data, to=c(0,1)) ## do away with negs
true_mean <- sum(sim_true)
true_sd <- sqrt(var(sim_true))
true_cv <- true_sd/true_mean

sim_sample <- subset(big_sims, simid == paste0('sim',isim) & cellid %in% stations & repid == irep)['x.data']
sim_sample <- rescale(sim_sample$x.data, to=c(0,1)) ## do away with negs
sim_mean <- mean(sim_sample) * total_area ## observed index
sim_sd <- sqrt(var(sim_sample)*total_area^2 / n_stations)
sim_cv <- sim_sd/sim_mean

rel_err_mean <- 100*(sim_mean-true_mean)/true_mean
rel_err_cv <- 100*(sim_cv-true_cv)/true_cv

big_sampler[idx,1] <- isim
big_sampler[idx,2] <- irep
big_sampler[idx,3:10] <- c(true_mean,true_sd,true_cv,
sim_mean,sim_sd,sim_cv, rel_err_mean,rel_err_cv)
idx = idx +1

}}

big_sampler_df <-data.frame(big_sampler)
names(big_sampler_df) <- c('simid', 'replicate', paste0('true_',c('mean','sd','cv')),
paste0('sim_',c('mean','sd','cv')), 'rel_err_mean','rel_err_cv')

#save(big_sampler_df, file = here('data','grf_big_sampler.rdata'))


## FIGURES

## Biomass plots ----
ggplot(sims, aes(x = X1, y= X2, fill = x.data, color =x.data)) +
geom_raster()+
coord_equal() +
scale_color_viridis_c()+
scale_fill_viridis_c()+
facet_wrap(~id, ncol = 2) +
theme_void() +
scale_y_discrete(expand = c(0,0))+
scale_x_discrete(expand = c(0,0))+
theme(legend.position = 'none', strip.text = element_blank())

ggsave(last_plot(), file = here('figs','grf_example.png'),
width = 8, height = 8, unit = 'in', dpi = 520)

## biomass plots with survey sites ----
ggplot(sims, aes(x = X1, y= X2, fill = x.data, color =x.data)) +
geom_raster()+
coord_equal() +
scale_color_viridis_c()+
scale_fill_viridis_c()+
facet_wrap(~id, ncol = 2) +
theme_void() +
geom_point(data=subset(sims, cellid %in% stations), 
pch = 4, color = 'blue', size = 3) +
scale_y_discrete(expand = c(0,0))+
scale_x_discrete(expand = c(0,0))+
theme(legend.position = 'none', strip.text = element_blank())

ggsave(last_plot(), file = here('figs','grf_example_25stations.png'),
width = 8, height = 8, unit = 'in', dpi = 520)


## Sampling error ----
big_sampler_df %>%
select(replicate, simid, sim_mean, true_mean, true_cv, sim_cv) %>%
reshape2::melt(id = c('replicate', 'simid')) %>%
mutate(
src = ifelse(grepl('sim',variable),'sim','true'),
variable= ifelse(grepl('cv',str_sub(variable,-4)),'cv',str_sub(variable,-4))) -> big_sampler_df_reshape

ggplot(data = subset(big_sampler_df_reshape,variable == 'mean'), 
aes(color = src))+
geom_boxplot(aes(x = factor(simid), y = value)) +
scale_color_manual(values = c('blue','grey22'), labels = c('Survey','Population') ) +
scale_x_discrete(labels = c('Low Morpho, Low Spatial', 'Low Morpho, High Spatial',
'High Morpho, Low Spatial','High Morpho, High Spatial'))+
labs(y = 'Mean Biomass in Domain', color = '')+
theme(axis.title.x = element_blank())

ggsave(last_plot(),
file = here('figs','mean_biomass_replicates_25stations.png'),
width = 6, height = 6, unit = 'in')


big_sampler_df_reshape %>% 
tidyr::pivot_wider(id_cols = c('replicate','simid','src'), 
values_from = 'value', names_from = 'variable') %>% 
mutate(lci = mean-mean*cv, uci = mean+mean*cv) -> big_sampler_df_reshape2


big_sampler_df_reshape2 %>%
group_by(simid, src) %>%
summarise(mean(mean),mean(cv))

big_sampler_df_reshape2 %>%
filter(replicate == 4) %>%
ggplot(., aes(x = factor(simid), color = src)) +
scale_color_manual(values = c('blue','grey22'), labels = c('Survey','Population') ) +
scale_x_discrete(labels = c('Low Morpho, Low Spatial', 'Low Morpho, High Spatial',
'High Morpho, Low Spatial','High Morpho, High Spatial'))+
geom_point(aes(y = mean)) +
labs(y = 'Biomass in Domain', color = '')+
geom_errorbar(aes(ymin = lci, ymax = uci), width = 0)
theme(axis.title.x = element_blank())

ggsave(last_plot(),
file = here('figs','surv_biomass_replicates_25stations.png'),
width = 6, height = 6, unit = 'in')

