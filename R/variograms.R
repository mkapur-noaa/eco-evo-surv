library(readr)
library(geoR)
# http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/c173c273_lec4.pdf
ftr <- list.files(files_to_run,
           pattern ="*true_biomass_ys\\b.csv",
           recursive = TRUE,
           full.names = TRUE)

acod_geodata <- read.csv(ftr[1]) %>% filter(year == 2040)
b <-as.geodata(acod_geodata[c('lat','long','tot_val')])
plot(b);points(b, cex.min=1, cex.max=3, col="gray")
# coords_use  =  as.matrix(cbind(acod_geodata$lat,acod_geodata$long), ncol = 2)
v1 <- geoR::variog( b)
plot(v1)


get_variogram_stats <- function(x){
  ## assuming x is the true_biomass_ys dataframe
  coords_use  =  as.matrix(cbind(x$lat,x$long), ncol = 2)
  v1 <- geoR::variog( coords = coords_use, data =x$tot_val)


  return(v1$v) ## variogram values at distance
}
