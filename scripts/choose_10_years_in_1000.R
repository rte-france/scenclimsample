library(tidyverse)
library(lubridate)
library(data.table)

source("R/get_score.R")
source("R/energy.R")

set.seed(2022)

nb_year <- 200
dt_series <- data.table(mcYear = rep(seq(nb_year), each = 8760),
                        horodate = rep(seq(ymd_h("2001-01-01 00"), ymd_h("2001-12-31 23"), by = "hour"), nb_year),
                         value_1 = rnorm(8760*nb_year),
                         value_2 = rnorm(8760*nb_year))

l_dt_series_sorted <- map(c("value_1", "value_2"), function(ccol){
  dt_dumb <- dt_series[,.SD,.SDcols = c("mcYear", ccol)]
  setnames(dt_dumb, ccol, "value")
  setorder(dt_dumb, value)
  setindex(dt_dumb, mcYear)
  })
names(l_dt_series_sorted) <- c("value_1", "value_2")

nrow_complet <- nrow(dt_series)
quantile_vec_allyears <- seq(1/nrow_complet*nb_year, 1-1/nrow_complet*nb_year, by = 1/nrow_complet*nb_year)

l_dt_series_quantiles <- map(l_dt_series_sorted, function(DT) DT[,.(value = quantile(value, quantile_vec_allyears, type = 3),
                                                       quantile_level = quantile_vec_allyears)])


dt_series_mean <- dt_series[,lapply(.SD, mean), .SDcols = -"mcYear"]

l_mat_A <- map(c("value_1", "value_2"), function(ccol)
  build_euclid_dist_matrix(dt_series[,.SD,.SDcols = c("mcYear", ccol)], block_id = "mcYear")
  )
names(l_mat_A) <- c("value_1", "value_2")


# l_mat_A <- list(zz = array(rnorm(1000*1000), c(1000, 1000)), zz2 = array(rnorm(1000*1000), c(1000, 1000)))

## liste des combinaisons à tester
Ntry <- 1e3
nb_year <- max(dt_series$mcYear)
nb_year_select <- 10

# combinaison aléatoire
array_sel_random <- sapply(seq(Ntry), function(x) sort(sample.int(nb_year, nb_year_select, replace = FALSE)))
colnames(array_sel_random) <- paste0("random_", seq(ncol(array_sel_random)))


system.time({dt_score_mean <- get_dt_score_mean(array_sel = array_sel_random, dt_series, dt_series_mean)})
system.time({dt_score_energy <- get_dt_score_energy(array_sel = array_sel_random, l_mat_A)})
system.time({dt_score_quantile <- get_dt_score_quantile(array_sel = array_sel_random, l_dt_series_sorted, l_dt_series_quantiles)})
