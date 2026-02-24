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
