library(tidyverse)
library(lubridate)
library(data.table)

# evtools::load_all()
source("R/get_score.R")
source("R/energy.R")

set.seed(2022)

nb_year <- 200
dt_series <- data.table(block_id = rep(seq(nb_year), each = 8760),
                        horodate = rep(seq(ymd_h("2001-01-01 00"), ymd_h("2001-12-31 23"), by = "hour"), nb_year),
                         value_1 = rnorm(8760*nb_year),
                         value_2 = rnorm(8760*nb_year))

l_mat_A <- map(c("value_1", "value_2"), function(ccol){

    dt_series_tmp <- dcast(dt_series[,.SD,.SDcols = c("block_id", "horodate", ccol)],
          date(horodate) + block_id ~ hour(horodate), value.var = ccol)

  build_euclid_dist_matrix(dt_series_tmp[,.SD,.SDcols = -("horodate")], block_id = "block_id")
  })
names(l_mat_A) <- c("value_1", "value_2")


# l_mat_A <- list(zz = array(rnorm(1000*1000), c(1000, 1000)), zz2 = array(rnorm(1000*1000), c(1000, 1000)))

## liste des combinaisons à tester
Ntry <- 1e3
nb_year <- max(dt_series$block_id)
nb_year_by_sample <- 10

# random samples
array_sample_random <- sapply(seq(Ntry), function(x) sort(sample.int(nb_year, nb_year_by_sample, replace = FALSE)))
colnames(array_sample_random) <- paste0("random_", seq(ncol(array_sample_random)))




system.time({dt_score_avg <- get_dt_score_avg(array_sample = array_sample_random, dt_series = dt_series[,.SD,.SDcols = -"horodate"])})
system.time({dt_score_energy <- get_dt_score_energy(array_sample = array_sample_random, l_mat_A = l_mat_A)})
system.time({dt_score_quantile <- get_dt_score_quantile(array_sample = array_sample_random, dt_series = dt_series, col_vec = c("value_1", "value_2"))})
