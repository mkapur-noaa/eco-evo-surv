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
# remotes::install_github("cararthompson/monochromeR"
library(monochromeR)

theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

## definitions and labels ----
scenLabs <-  c('No Climate Change & No Evolution',
'No Climate Change & Evolution',
'Climate Change & No Evolution',
'Climate Change & Evolution')
names(scenLabs) <- c('noCC_noEvo','noCC_Evo','CC_noEvo','CC_Evo')



sppLabs <-c("Atlantic Herring" , "Atlantic Mackerel" ,"Sandeel" ,"European Sprat","Norway Pout" ,
"European Plaice","Common Sole","Saithe","Atlantic Cod","Haddock","Horse Mackerel" ,
"Whiting","Common Dab","Grey Gurnard" ,"Hake","Shrimp")
#names(sppLabs) <- as.character(0:15)
names(sppLabs) <- gsub(" ","",sppLabs)
sppLabs2 <- cbind(Var1 = sppLabs, Var2=names(sppLabs), Var3 = 0:(length(sppLabs)-1))


## color palettes ----
scenPal <- c('#a00000','#d8a6a6','#9fc8c8','#298c8c') ## four scenario colors: reds have climate change, blues do not
scenLabs2 <- cbind(Var1 = scenLabs, Var2=names(scenLabs), Var3 = 1:4, Pal = rev(scenPal))
sppPal <- gplots::rich.colors(n = 13)

survey_array <- build_survey_array(fractional_coverage = 0.3)


## build and save parameter lookup tables ----

# parm_path <- "F:/Ev-osmose/Ev-OSMOSE outputs_15April2024/cc_evo/ns_param-species.csv" ## stored once for all runs, not time-invariant

#* max age (lifespan)
read.csv(parm_path, header = F) %>%
  filter(grepl('lifespan.sp',V1)) %>%
  filter(!grepl('#species',V1)) %>% ## drop commented-out ones
  mutate(sppIdxTEMP = as.numeric(stringr::str_extract(V1, "\\d+$"))) %>%
  merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
  select(species0 =  sppIdxTEMP, species = Var2, value = V2) %>%
  write.csv(here::here('outputs','wham_runs','max_age.csv'), row.names = FALSE)


#* length-at-age von B pars ----
# read.csv(parm_path, header = F) %>%
# filter(grepl('length2weight',V1)) %>%
#  mutate(sppIdxTEMP = stringr::str_extract(V1, "\\d+$"),
#   variable = stringr::str_extract(V1, "(?<=length2weight\\.)\\w+(?=\\.)")) %>%
# merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
#  select(species0 =  sppIdxTEMP, species = Var2, variable, value = V2) %>%
#  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
#  write.csv(here::here('outputs','wham_runs','length2weight.csv'), row.names = FALSE)


#* length-weight allometric pars ----
# read.csv(parm_path, header = F) %>%
# filter(grepl('length2weight',V1)) %>%
#  mutate(sppIdxTEMP = stringr::str_extract(V1, "\\d+$"),
#   variable = stringr::str_extract(V1, "(?<=length2weight\\.)\\w+(?=\\.)")) %>%
# merge(., sppLabs2, by.x = 'sppIdxTEMP', by.y = 'Var3') %>%
#  select(species0 =  sppIdxTEMP, species = Var2, variable, value = V2) %>%
#  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
#  write.csv(here::here('outputs','wham_runs','length2weight.csv'), row.names = FALSE)
