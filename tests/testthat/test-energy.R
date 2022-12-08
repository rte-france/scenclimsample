mat_init <- t(array(c(3, 4, 5, 2, 1, 1, 1, 0), c(4, 2)))
mat_avg <- t(array(c(3.5, 3.5, 1, 0.5), c(2, 2)))


test_that("matrixblockmean works", {
  expect_equal(matrixblockmean(mat_init, 1, 2), mat_avg)
})



dt_series <- data.table(block_id = rep(seq(2), each=3),
                        value00 = c(rep(0, 3), seq(3)))

aa <- expand.grid(a = seq(3), b = seq(3))
dist22 <- mean(abs(aa$a - aa$b))
mat_A_ref <- array(t(c(0, 2, 2, dist22)), c(2, 2))

build_euclid_dist_matrix(dt_series)

test_that("matrixblockmean works", {
  expect_equal(build_euclid_dist_matrix(dt_series), mat_A_ref)
})
