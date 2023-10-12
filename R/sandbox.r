require(dplyr)
require(ggplot2)
library(BAT)
require(here)
## these are dependencies for mixfishsim that have been archived
#install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz", repos=NULL)
#install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.14.tar.gz", repos=NULL)
# devtools::install_github('pdolder/mixfishsim')
## there were bugs so i forked it
devtools::install_local("C://Users/maia.kapur/Work/projects/eco-evo-surv/MixFishSim")
library(MixFishSim)

## working with the mixfishsim vignette  https://github.com/pdolder/MixFishSim/blob/master/vignettes/Simple_MixFishSim_Example.Rmd
sim <- init_sim(nrows = 10, ncols = 10, n_years = 10, n_tows_day = 4,
		n_days_wk_fished = 5, n_fleets = 2, n_vessels = 20, n_species =
			2, move_freq = 2)

class(sim)
sim$idx
names(sim$brk.idx)

habby <- create_hab(sim_init = sim, 
		  spp.ctrl = list(
				  "spp.1" = list('nu' = 1/0.015,
						'var' = 1,
						'scale' = 1,
						'Aniso' = matrix(nc = 2, c(1.5, 3, -3, 4))),
				  "spp.2" = list('nu' = 1/0.05,
						'var'  = 2,
						'scale' = 12,
						'Aniso' = matrix(nc = 2, c(1, 2, -1, 2)))
				  ),
		  spawn_areas = list(
				     "spp1" = list(
						   'area1' = c(2,5,2,5),
						   'area2' = c(6,8,6,8)
						   ),
				     "spp2" = list(
						   'area1' = c(5,6,6,6)
						   )),
				     spwn_mult = 10, 
				     plot.dist = FALSE)

print(habby)

## Plot the unadjusted habitat fields
plot_habitat(habby$hab)

## Plot the adjusted habitat fields
plot_habitat(habby$spwn_hab)


Pop <- init_pop(sim_init = sim, 
        Bio = c("spp1" = 1e5, "spp2" = 1e5),
		hab = habby[["hab"]], 
        start_cell = c(5,5),
		lambda = c("spp1" = 0.1, "spp2" = 0.1),
		init_move_steps = 20,
		rec_params = list("spp1" = c("model" = "BH", "a" = 54, "b" = 2, "cv" = 0.7),
				  "spp2" = c("model" = "BH", "a" = 27, "b" = 4,"cv" = 0.3)),
		rec_wk = list("spp1" = 3:6, "spp2" = 4:8),
		spwn_wk = list("spp1" = 4:8, "spp2" = 4:8),
		M = c("spp1" = 0.2, "spp2" = 0.2),
		K = c("spp1" = 0.3, "spp2" = 0.3),
        wt = c("spp1" = 1, "spp2" = 1), ## these default to the wrong dimlength
        wtm1 = c("spp1" = 1, "spp2" = 1) ## these default to the wrong dimlength
		)

names(Pop)

Pop$dem_params

image(Pop$Start_pop[[1]], main = "spp1 starting biomass")
image(Pop$Start_pop[[2]], main = "spp2 starting biomass")

## fish fleet (storage)
fleets <- init_fleet(sim_init = sim, 
            VPT = list("spp1" = 4, "spp2" = 3),
		     Qs = list("fleet 1" = c("spp1" = 1e-5, "spp2" = 3e-5),
			       "fleet 2" = c("spp1" = 5e-5, "spp2" = 1e-5)
			       ),
		     fuelC = list("fleet1" = 3, "fleet 2" = 8),
		     step_params = list("fleet 1" = c("rate" = 3, "B1" = 1, "B2" = 2, "B3" = 3),
					"fleet 2" = c("rate" = 3, "B1" = 2, "B2" = 4, "B3" = 4)
					),				
		     past_knowledge = TRUE,
		     past_year_month = TRUE,
		     past_trip = TRUE,
		     threshold = 0.7
		     )

test_step(step_params = fleets$fleet_params[[1]]$step_params, rev.max = 1e2)
test_step(step_params = fleets$fleet_params[[2]]$step_params, rev.max = 1e2)

## Movement Model given Temperature
moveCov <- init_moveCov(sim_init = sim, steps = 52,
			spp_tol = list("spp1" = list("mu" = 12, "va" = 8),
				       "spp2" = list("mu" = 15, "va" = 7)
				       )
			)

plot(norm_fun(x = 0:25, mu = 12, va = 8)/max(norm_fun(0:25, 12, 8)), 
     type = "l", xlab = "Temperature", ylab = "Tolerance", lwd = 2)
lines(norm_fun(x = 0:25, mu = 15, va = 7)/ max(norm_fun(0:25, 15, 7)),
      type = "l", col = "blue", lwd = 2)
legend(x = 2, y = 0.9, legend = c("spp1", "spp2"), lwd = 2, col = c("black", "blue"))

plot_spatiotemp_hab(hab = habby, moveCov = moveCov, 
spwn_wk = list("spp1" = 4:8, "spp2" = 4:8), 
plot.file = here::here('figs'))

## initialize the survey
surv <- init_survey(sim_init = sim, 
design = 'fixed_station', 
n_stations = 50,
start_day = 1, stations_per_day = 5, 
Qs= rep(1,2))

## the simulator looks for 'Qs.spp1'
names(surv$survey_settings)[4:5] <- paste0('Qs.spp',seq(2))

fleet_params = fleets[["fleet_params"]][[1]]
fleet_catches = fleets[["fleet_catches"]][[1]][[1]] ## emtpy matrix
sp_fleet_catches = fleets[["sp_fleet_catches"]][[1]][[1]]

## first fleet, first vessel, tow t
catches[[1]][['fleets_catches']][[1]][1:5,]
catches[[1]][['sp_fleets_catches']][[1]][1:5,]

reshape2::melt(catches[[1]][['sp_fleets_catches']]) %>%
filter(Var2 == 'spp1' & !is.na(value)) %>% 
sum(.$value)
reshape2::melt(catches[[2]][['sp_fleets_catches']]) %>%
filter(Var2 == 'spp1' & !is.na(value)) %>% 
sum(.$value)

as.list(rowSums(sapply(catches[[1]][['fleets_catches']], unlist)))

str(catches[[1]][['fleets_catches']]) 
sum_fleets_catches(sim_init = sim_init, 
                fleets_log = catches)

sum_fleet_catches(fleet_log = catches[[1]], sim_init = sim_init)

sum_fleet_catches <- function (sim_init = sim, fleet_log = NULL) {
    n_spp <- sim_init[["idx"]][["n.spp"]]
    y <- lapply(seq(n_spp), function(x) {
        z <- lapply(fleet_log[["fleets_catches"]], "[[", x)
        return(z)
    })
    flt_catches <- lapply(y, function(x) {
        Reduce("+", x)
    })
    names(flt_catches) <- paste("spp", seq(n_spp), sep = "")
    return(flt_catches)
}


## run the simulation
res <- run_sim(sim_init = sim,
	       pop_init = Pop,
	       move_cov = moveCov,
	       fleets_init = fleets,
	       hab_init = habby,
	       save_pop_bio = TRUE,
	       survey = surv,
	       closure = NULL)



## PLOTS           
p1 <- plot_pop_summary(results = res,
		 timestep = "annual",
		 save = FALSE
		 )
p1

p2 <- plot_daily_fdyn(res)
p2

## Fishery

logs <- combine_logs(res[["fleets_catches"]])

p3 <- plot_vessel_move(sim_init = sim, logs = logs, fleet_no = 1, vessel_no = 5,
		 year_trip = 5, trip_no = 10)
p3

p4 <- plot_realised_stepF(logs = logs, fleet_no = 1, vessel_no = 1)
p4    


## playing with survey sampler
load(here::here('data','optimization_data.rdata'))
source("https://raw.githubusercontent.com/zoyafuso-NOAA/chukchi_survey_evaluation/main/code/01_analysis/simulate_survey.R")
## DEPRECATED----

## mockup study design
grid_dim <- 10
sptl_grid <- matrix(0, nrow = grid_dim, ncol = grid_dim)

par(mfrow = c(2 ,2))

for(i in 1:4){
	dist_use <- c('aggregated','random','uniform','gradient')[i]
    comm <- BAT::sim.spatial(100, 1, distribution = dist_use)
	plot(comm$x, comm$y, main = dist_use, xlim = c(0,1), ylim = c(0,1))
}

 