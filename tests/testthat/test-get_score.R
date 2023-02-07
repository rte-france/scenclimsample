# Copyright (C) 2022-2023, RTE (http://www.rte-france.com)
# SPDX-License-Identifier: MPL-2.0

## Test moyenne par colonne
dt_series <- data.table(block_id = rep(seq(2), each=3),
                        value01 = c(rep(0, 3), seq(3)))

array_sample <- array(c(1, 2), c(1, 2))
colnames(array_sample) <- c("sampletest1", "sampletest2")
dt_score_avg_ref <- data.table(value01 = c(1, 1), name_sample = c("sampletest1", "sampletest2"))

test_that("get_dt_score_avg works", {
  expect_equal(get_dt_score_avg(array_sample, dt_series), dt_score_avg_ref)
})


## Test get_dt_score_quantile_one_col
dt_series2 <- data.table(block_id = rep(seq(2), each=1001),
                        value01 = c(rep(0, 1001), seq(1001)))
dt_score_quantile_ref <- structure(list(mae_q_value01 = c(250, 250.5),
                                        erreur_q99_value01 = c(979.02, 495.01),
                                        mae_pointe_value01 = c(500, 250.5),
                                        name_sample = c("sampletest1", "sampletest2")),
                                   row.names = c(NA, -2L), class = c("data.table", "data.frame"))

test_that("get_dt_score_quantile works", {
    expect_equal(get_dt_score_quantile(array_sample, dt_series2, col_vec = c("value01"), upper_quantile_boundary = 0.5),
                 dt_score_quantile_ref)
})

## Test get_dt_score_energy
dt_series3 <- data.table(block_id = rep(seq(2), each=1001),
                         value01 = c(rep(0, 1001), seq(1001)))

l_mat_A <- list(dumbA = build_euclid_dist_matrix(dt_series3))

dt_score_energy_ref <- structure(list(dumbA = c(83.5417082917083, 83.5417082917082),
               name_sample = c("sampletest1", "sampletest2")),
          class = c("data.table", "data.frame"), row.names = c(NA, -2L))

test_that("get_dt_score_energy works", {
  expect_equal(get_dt_score_energy(array_sample, l_mat_A),
               dt_score_energy_ref)
})


