## eco-evo-surv settings
## establish packages, plotting settings, and directories


#0=AtlanticHerring 1=AtlanticMackerel 2=Sandeel 3=EuropeanSprat 4=NorwayPout
# 5=EuropeanPlaice 6=CommonSole 7=Saithe 8=AtlanticCod 9=Haddock 10=HorseMackerel
##11=Whiting 12=CommonDab 13=GreyGurnard 14=Hake 15=Shrimp

## packages ----
#https://pjbartlein.github.io/REarthSysSci/netCDF.html#reading-restructuring-and-writing-netcdf-files-in-r
#install.packages('ncdf4')
library(ncdf4)
library(dplyr)
library(here)
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library(testthat)
library(wham)
library(mgcv)
library(gratia)
# remotes::install_github("cararthompson/monochromeR"
library(monochromeR)
library(RColorBrewer)
library(dismo)

theme_set(ggsidekick::theme_sleek()+
            theme(
                  legend.background = element_blank(),
                  plot.background = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
options(dplyr.summarise.inform = FALSE)

## definitions and labels ----
scenLabs <-  c('No Climate Change & No Evolution',
'No Climate Change & Evolution',
'Climate Change & No Evolution',
'Climate Change & Evolution')
names(scenLabs) <- c('noCC_noEvo','noCC_Evo','CC_noEvo','CC_Evo')


fcLabs <- c('perfect information', '50% Coverage', '30% Coverage',
            'Simple Survey','Simple Survey + Selex')
names(fcLabs)  <- c(1,0.5,0.3,0.15,0.15001)


sppLabs <-c("Atlantic Herring" , "Atlantic Mackerel" ,"Sandeel" ,"European Sprat","Norway Pout" ,
"European Plaice","Common Sole","Saithe","Atlantic Cod","Haddock","Horse Mackerel" ,
"Whiting","Common Dab","Grey Gurnard" ,"Hake","Shrimp")
#names(sppLabs) <- as.character(0:15)
names(sppLabs) <- gsub(" ","",sppLabs)

sppLabs2 <- data.frame(Var1 = sppLabs,
                  Var2=names(sppLabs),
                  Var3 = 0:(length(sppLabs)-1),
                  Var4 = FALSE)
sppLabs2$Var5[sppLabs2$Var3%in%c(0,3,7,8)] <- c('MidForage','BaseForage','Predator1','Predator2') ## mask names for spp
sppLabs2$Var4[sppLabs2$Var3%in%c(0,3,7,8)] <- TRUE ## indicate spp of interest
sppLabs[sppLabs2$Var3%in%c(0,3,7,8)] <- sppLabs2$Var5[sppLabs2$Var3%in%c(0,3,7,8)] ## mask names for spp


## color palettes ----
scenPal <- c('#a00000','#d8a6a6','#9fc8c8','#298c8c') ## four scenario colors: reds have climate change, blues do not

scenLabs2 <- data.frame(cbind(Var1 = scenLabs, Var2=names(scenLabs), Var3 = 1:4, Pal = rev(scenPal)))
scenLabs2$Pal2 <- c('#1D6060','#1D6060', '#6A0000','#6A0000')

sppPal <- gplots::rich.colors(n = 13)

## function from Cole to calculate actual q used given bounds (logit)
q_f <- function(x,a=0,b=1000) a+(b-a)/(1+exp(-x))
## get q from the outputs (from wham plotting code; equivalent to above)
# se = t(matrix(as.list(mod_use$sdrep, "Std. Error")$logit_q, nrow = NCOL(mod_use$rep$logit_q_mat),
#               ncol = NROW(mod_use$rep$logit_q_mat)))
# logit_q_lo = mod_use$rep$logit_q_mat - qnorm(0.975)*se
# logit_q_hi = mod_use$rep$logit_q_mat + qnorm(0.975)*se
# q = mean(t(mod_use$input$data$q_lower + (mod_use$input$data$q_upper - mod_use$input$data$q_lower)/(1+exp(-t(mod_use$rep$logit_q_mat)))))

## build survey arrays at many coverage types for lookup
## this ensures that for a given coverage fraction, the exact same design
## is used for every species, scenario, and replicate
# for(i in c(seq(0.05,1,0.05))){
#   cat(i,"\n")
#   build_survey_array(fractional_coverage = i)
# }

# https://stackoverflow.com/questions/48297440/list-files-recursive-up-to-a-certain-level-in-r
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


## for use in optimizers
logistic <- function(x,x0,k){
  return(1/(1+exp(-k*(x-x0))))
}
# cbind(age = 1:max_age_pop, slx =   1/(1+exp(-log(19)*((1:max_age_pop)-srv_selex)/(max_age_pop-srv_selex))))

logistic2 <- function(x,x0,amax){
  return( 1/(1+exp(-log(19)*((x)-x0)/(amax-x0))))
}

maturity_optim_fn <- function(pars = c(x0,k), md = maturity_data, amax = max_age_pop){
  x0_use = pars[1];k_use = pars[2]
  mat_temp <- logistic(x = 1:amax, x0=x0_use, k = k_use)
  obj <- 0
  for(i in 1:90){ ## terminal year has an NA
    obj <- obj+sum((md[i,] - mat_temp)^2)
  }
  return(obj)
}

## index for alphabetized reps
ordered_reps0 <- data.frame(lab = paste0('rep',0:27)) %>%
  arrange(lab) %>%
  mutate(idx0 = as.numeric(gsub('rep','',lab)), ## true replicate ID
         idx = 0:27, ## ordered list
         idx1 = idx+1) %>% ## index-friendly
  arrange(idx0)

ordered_reps <- ordered_reps0 %>%
  dplyr::select(idx1) %>%
  t()

# maturity_optim_fn <- function(pars){
#   x0_use = pars[1];
#   mat_temp <- logistic2(x = 1:max_age_pop, x0=x0_use, amax = max_age_pop)
#   obj <- 0
#   for(i in 1:90){ ## terminal year has an NA
#     obj <- obj+sum((maturity_data[i,] - mat_temp)^2)
#   }
#   return(obj)
# }
## build and save parameter lookup tables ----

# parm_path <- "F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/cc_evo/ns_param-species.csv" ## stored once for all runs, not time-invariant


#* putative maturity curves


#* max age (lifespan)
# read.csv(parm_path, header = F) %>%
#   filter(grepl('lifespan.sp',V1)) %>%
#   filter(!grepl('#species',V1)) %>% ## drop commented-out ones
#   mutate(sppIdxTEMP = as.numeric(stringr::str_extract(V1, "\\d+$"))) %>%
#   merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
#   select(species0 =  sppIdxTEMP, species = Var2, value = V2) %>%
#   write.csv(here::here('outputs','wham_runs','max_age.csv'), row.names = FALSE)


#* length-at-age von B pars ----
# read.csv(parm_path, header = F) %>%
#   filter(., grepl(paste(c('t0','lInf','K'), collapse="|"), V1)) %>%
#  mutate(sppIdxTEMP = stringr::str_extract(V1, "\\d+$"),
#   variable = stringr::str_extract(V1, "(?<=\\.)\\w+(?=\\.)")) %>%
# merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
#  select(species0 =  sppIdxTEMP, species = Var2, variable, value = V2) %>%
#  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
#  write.csv(here::here('outputs','wham_runs','vonBpars.csv'), row.names = FALSE)


#* length-weight allometric pars ----
# read.csv(parm_path, header = F) %>%
# filter(grepl('length2weight',V1)) %>%
#  mutate(sppIdxTEMP = stringr::str_extract(V1, "\\d+$"),
#   variable = stringr::str_extract(V1, "(?<=length2weight\\.)\\w+(?=\\.)")) %>%
# merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
#  select(species0 =  sppIdxTEMP, species = Var2, variable, value = V2) %>%
#  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
#  write.csv(here::here('outputs','wham_runs','length2weight.csv'), row.names = FALSE)
