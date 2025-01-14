#* input WAA at maturity across scenario----
strip_waa <- function(x) {
  tmpwaa <- read.table(x, sep = ' ', header = TRUE)
  names(tmpwaa) <- 1:ncol(tmpwaa)
  tmpwaa$year <- 2010:(2009 + nrow(tmpwaa))

  filen2 = basename(x)
  # filen2 <- basename(dirname(dirname(file_suffix)))
  tmpwaa$scenario <- strsplit(filen2, '-')[[1]][2]
  tmpwaa$rep <- strsplit(filen2, '-')[[1]][3]
  tmpwaa$species <-   strsplit(filen2, '-')[[1]][1]
  waa <-
    reshape2::melt(tmpwaa, id = c('year', 'rep', 'scenario', 'species'))
  names(waa)[5:6] <- c('age', 'weight_kg')
  waa$age <- as.numeric(waa$age)
  return(waa)
}
