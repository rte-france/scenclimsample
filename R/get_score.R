#' Compute average mean absolute error for all columns
#'
#' @param array_sample array of samples to test, one sample by column
#' @param dt_series data.table of series
#'
#' @return a data.table of scores, one row per sample
#'
#' @export
#' @import data.table
#'
get_dt_score_avg <- function(array_sample, dt_series){

  # Average by block and average over all blocks
  dt_avg_over_all_blocks <- dt_series[, colMeans(.SD), .SDcols = -"block_id"]
  dt_avg <- dt_series[, colMeans(.SD), by = "block_id"]
  setorder(dt_avg, block_id)

  # loop over samples
  dt_avg_sample <- purrr::map_dfr(seq(ncol(array_sample)), function(ii){
    dt_avg[array_sample[, ii], colMeans(.SD), .SDcols = -"block_id"]
  })
  dt_avg_sample$name_sample = colnames(array_sample)
  setDT(dt_avg_sample)

  # scoring
  dt_score_avg <- abs(dt_avg_over_all_blocks[rep(1, nrow(dt_avg_sample))] - dt_avg_sample[,.SD, .SDcols = -"name_sample"])
  dt_score_avg[,name_sample := dt_avg_sample$name_sample]
  dt_score_avg
}

#' Compute energy scores
#'
#' @param array_sample array of samples to test, one sample by column
#' @param l_mat_A list of A matrices for the energy scores
#'
#' @return a data.table of scores, one row per sample
#'
#' @export
#' @import data.table
#'
get_dt_score_energy <- function(array_sample, l_mat_A){

  # Loop over matrices A
  dt_score_energy <- purrr::map_dfc(l_mat_A, function(mat_A){

    Nyear_total <- nrow(mat_A)
    w0 <- array(1/Nyear_total, c(Nyear_total,1))

    # precomputation
    mat_A_w0 <- mat_A %*% w0
    w0_mat_A_w0_half <- 0.5*t(w0) %*% mat_A %*% w0

    # loop over samples
    vec_score_energy <- purrr::map_dbl(seq(ncol(array_sample)), function(ii){
      ws <- array(0, c(Nyear_total,1))
      ws[array_sample[, ii],1] <- 1/Nyear_total

      # energy score
      t(ws) %*% mat_A_w0 - w0_mat_A_w0_half - 0.5*t(ws) %*% mat_A %*% ws
    })

    data.frame(t(vec_score_energy))
  })

  setDT(dt_score_energy)
  setnames(dt_score_energy, names(l_mat_A))
  dt_score_energy[,name_sample := colnames(array_sample)]
}


#' Compute quantile scores for a single column
#'
#' @param array_sample array of samples to test, one sample by column
#' @param dt_series_one_col data.table with two columns block_id and value
#' @param upper_quantile_boundary num boundary for specific score on the upper quantiles
#'
#' @return a data.table of scores, one row per sample
#'
#' @export
#' @import data.table
#'
get_dt_score_quantile_one_col <- function(array_sample, dt_series_one_col, upper_quantile_boundary = 0.998){

  # dt_quantiles computed with all blocks
  nrow_all_blocks <- nrow(dt_series_one_col)
  nb_blocks <- length(unique(dt_series_one_col$block_id))
  quantile_vec_all_blocks <- seq(1/nrow_all_blocks*nb_blocks, 1-1/nrow_all_blocks*nb_blocks, by = 1/nrow_all_blocks*nb_blocks)

  dt_quantiles <- dt_series_one_col[,.(value = stats::quantile(value, quantile_vec_all_blocks, type = 3),
                                       quantile_level = quantile_vec_all_blocks)]

  # sort dt_series
  dt_series_sorted <- copy(dt_series_one_col)
  setnames(dt_series_sorted, ccol, "value")
  setorder(dt_series_sorted, value)
  setindex(dt_series_sorted, block_id)

  # known quantile positions in the sample
  nrow_sample <- nrow(dt_series_sorted[.(array_sample[,1]), on = "block_id"])
  nb_blocks_sample <- nrow(array_sample)
  quantile_vec_sample_blocks <- seq(1/nrow_sample*nb_blocks_sample, 1-1/nrow_sample*nb_blocks_sample, by = 1/nrow_sample*nb_blocks_sample)
  quantile_position_index <- round(quantile_vec_sample_blocks*nrow_sample)
  quantile_position_upper_index <- which(quantile_vec_sample_blocks>upper_quantile_boundary)

  # loop over samples
  dt_score_quantile <- purrr::map_dfr(seq(ncol(array_sample)), function(ii){

    # estimated quantiles of the sample
    value_estim <- sort(dt_series_sorted[.(array_sample[,ii]), on = "block_id"]$value)[quantile_position_index]
    diff_vec <- dt_quantiles$value - value_estim

    # scores of the estimated quantiles
    data.table(mae_q = mean(abs(diff_vec)),
               erreur_q99 = stats::quantile(abs(diff_vec), 0.99),
               mae_pointe = mean(abs(diff_vec[quantile_position_upper_index])))
  })

}



#' Compute and concatenate quantile scores for multiple columns
#'
#' @param array_sample array of samples to test, one sample by column
#' @param dt_series data.table of series
#' @param col_vec a character vector with columns names for quantile score computation
#' @param ... Other arguments passed to internal get_dt_score_quantile_one_col
#'
#' @return a data.table of scores, one row per sample
#'
#' @export
#' @import data.table
#'
get_dt_score_quantile <- function(array_sample, dt_series, col_vec, ...){

  dt_score_quantiles <- purrr::reduce(purrr::map(col_vec, function(ccol){

    dt_series_one_col <- dt_series[,.SD, .SDcol = c("block_id", ccol)]
    setnames(dt_series_one_col, ccol, "value")
    dt_score_quantiles_one_col <- get_dt_score_quantile_one_col(array_sample, dt_series_one_col, ...)
    names(dt_score_quantiles_one_col) <- paste0(names(dt_score_quantiles_one_col), "_", ccol)
    dt_score_quantiles_one_col
  }),
  cbind)
  dt_score_quantiles[,name_sample := colnames(array_sample)]

}

