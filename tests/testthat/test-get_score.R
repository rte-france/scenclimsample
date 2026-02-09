# Copyright (C) 2022-2023, RTE (http://www.rte-france.com)
# SPDX-License-Identifier: MPL-2.0

array_sample <- array(c(1, 2), c(1, 2))
colnames(array_sample) <- c("sampletest1", "sampletest2")

test_that("get_dt_score_avg works : single column", {
  dt_series <- data.table(block_id = rep(seq(2), each=3),
                          value01 = c(rep(0, 3), seq(3)))

  dt_score_avg_ref <- data.table(value01 = c(1, 1), name_sample = c("sampletest1", "sampletest2"))

  expect_equal(get_dt_score_avg(array_sample, dt_series), dt_score_avg_ref)
})

test_that("get_dt_score_avg works : multiple column", {
  dt_series <- data.table(block_id = rep(seq(2), each=3),
                          value01 = c(rep(0, 3), seq(3)),
                          value02 = c(rep(0, 3), seq(3)*10))

  dt_score_avg_ref <- data.table(value01 = c(1, 1), value02 = c(10, 10), name_sample = c("sampletest1", "sampletest2"))

  expect_equal(get_dt_score_avg(array_sample, dt_series), dt_score_avg_ref)
})


## Test get_dt_score_quantile

dt_series2 <- data.table(block_id = rep(seq(2), each=1000),
                        value01 = c(rep(0, 1000), seq(1000)))
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


test_that("build_l_block_cdf returns correct structure and values", {
  dt <- data.table(block_id = rep(1:2, each=4), v1 = 1:8, v2 = 8:1)
  l_ref <- build_l_block_cdf(dt, col_vec = c("v1", "v2"))
  expect_type(l_ref, "list")
  expect_equal(names(l_ref), c("v1", "v2"))
  expect_equal(l_ref$v1[,1], l_ref$v2[,2])

  n_all <- 8
  nb_blocks <- 2
  quantile_vec_all <- seq(1/(n_all)*nb_blocks/2, 1-1/(n_all)*nb_blocks/2, by=1/(n_all)*nb_blocks)
  quantile_vec_all_value <- quantile(dt$v1, quantile_vec_all)
  ecdf_v1 <- ecdf(dt$v1)

  for (ii in seq(4)){
    for(bb in seq(2)){
      ecdf_v1_block_b <- ecdf_v1(dt[block_id == bb, v1])
      expect_equal(l_ref$v1[ii,bb], mean(ecdf_v1_block_b <= quantile_vec_all[ii]))
    }
  }

})


test_that("get_dt_score_CDF_parallel works on toy example", {
  dt <- data.table(block_id = rep(1:2, each=4), v1 = 1:8, v2 = 8:1)
  l_ref <- build_l_block_cdf(dt, col_vec = "v1")
  array_sample <- array(1:2, c(1,2)); colnames(array_sample) <- c("sampletest1", "sampletest2")
  res <- get_dt_score_CDF(array_sample, l_ref)
  expect_s3_class(res, "data.table")
  expect_equal(res$name_sample, c("sampletest1","sampletest2"))
})

