# Copyright (C) 2022-2023, RTE (http://www.rte-france.com)
# SPDX-License-Identifier: MPL-2.0

#' Compute average mean absolute error for all columns
#'
#' @param array_sample array of samples to test, one sample by column
#' @param dt_series data.table of series
#' @param col_vec a character vector with columns names for quantile score computation
#'
#' @return a data.table of scores, one row per sample
#'
#' @export
#' @import data.table
#'
get_dt_score_avg <- function(array_sample, dt_series, col_vec = NULL){

  if (is.null(col_vec)){
    col_vec <- grep("^block_id$", names(dt_series), value = TRUE, invert = TRUE)
  }

  dt_avg_over_all_blocks <- dt_series[, lapply(.SD, mean), .SDcols = col_vec]

  dt_avg <- dt_series[, lapply(.SD, mean), by = block_id]
  setorder(dt_avg, block_id)

  mat_avg <- as.matrix(dt_avg[, ..col_vec])
  vec_avg_all <- as.numeric(dt_avg_over_all_blocks[1, ..col_vec])

  mat_score <- sapply(seq_len(ncol(array_sample)), function(ii) {
    abs(colMeans(mat_avg[array_sample[, ii], , drop = FALSE]) - vec_avg_all)
  })

  if(length(col_vec)>1){
    dt_score_avg <- as.data.table(t(mat_score))
  } else {
    dt_score_avg <- as.data.table(list(V1 = mat_score))
  }

  setnames(dt_score_avg, col_vec)
  dt_score_avg[, name_sample := colnames(array_sample)]
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
  dt_score_energy <- purrr::imap_dfc(l_mat_A, function(mat_A, mat_A_name){

    Nyear_total <- nrow(mat_A)
    Nyear_sample <- nrow(array_sample)
    w0 <- array(1/Nyear_total, c(Nyear_total,1))

    # precomputation
    mat_A_w0 <- mat_A %*% w0
    w0_mat_A_w0_half <- 0.5*t(w0) %*% mat_A %*% w0

    # loop over samples
    vec_score_energy <- purrr::map_dbl(seq(ncol(array_sample)), function(ii){

      # # energy score V1
      # ws <- array(0, c(Nyear_total,1))
      # ws[array_sample[, ii],1] <- 1/Nyear_sample
      # t(ws) %*% mat_A_w0 - w0_mat_A_w0_half - 0.5*t(ws) %*% mat_A %*% ws

      # # energy score V2
      sum(mat_A_w0[array_sample[, ii]])/Nyear_sample - w0_mat_A_w0_half - 0.5*sum(mat_A[array_sample[, ii], array_sample[, ii]])/(Nyear_sample)^2

      # # energy score V3
      # ws_sp <- Matrix::Matrix(ws, sparse = T)
      # as.numeric(Matrix::t(ws_sp) %*% mat_A_w0 - w0_mat_A_w0_half - 0.5*Matrix::t(ws_sp) %*% mat_A %*% ws)

    })

    list_output <- list()
    list_output[[mat_A_name]] <- vec_score_energy
    data.frame(list_output)
  })

  setDT(dt_score_energy)
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

  # sort dt_series
  dt_series_sorted <- copy(dt_series_one_col)
  setorder(dt_series_sorted, value)
  setindex(dt_series_sorted, block_id)

  # quantiles computed with all blocks
  nrow_all_blocks <- nrow(dt_series_one_col)
  nb_blocks <- length(unique(dt_series_one_col$block_id))

  nrow_sample <- nrow(dt_series_sorted[.(array_sample[,1]), on = "block_id"])
  nb_blocks_sample <- nrow(array_sample)

  nb_quantiles <- floor(nrow_all_blocks/nb_blocks)
  quantile_vec <- seq(round(nb_blocks/2), nrow_all_blocks-nb_blocks/2, by = nb_blocks)/nrow_all_blocks
  quantile_position_upper_index <- which(quantile_vec>upper_quantile_boundary)

  quantile_position_index_all_blocks <- quantile(seq(nrow_all_blocks), quantile_vec, type = 1)
  quantile_position_index <- quantile(seq(nrow_sample), quantile_vec, type = 1)

  value_obs <- dt_series_sorted$value[quantile_position_index_all_blocks]

  # loop over samples
  dt_score_quantile <- purrr::map_dfr(seq(ncol(array_sample)), function(ii){

    # estimated quantiles of the sample
    value_estim <- sort(dt_series_sorted[.(array_sample[,ii]), on = "block_id"]$value)[quantile_position_index]
    diff_vec <- value_obs - value_estim

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
get_dt_score_quantile <- function(array_sample, dt_series, col_vec = NULL, ...){

  if (is.null(col_vec)){
    col_vec <- grep("^block_id$", names(dt_series), value = TRUE, invert = TRUE)
  }

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

#' Compute and concatenate quantile scores for multiple columns, fast version
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
get_dt_score_quantile_fast <- function(array_sample, dt_series, col_vec = NULL, upper_quantile_boundary = 0.998) {
  if (is.null(col_vec)) {
    col_vec <- setdiff(names(dt_series), "block_id")
  }

  if(ncol(array_sample)>0){

    outlist <- furrr::future_map(col_vec, function(ccol) {
      value_vec <- dt_series[[ccol]]
      dtval <- data.table(block_id = dt_series$block_id, value = value_vec)
      setkey(dtval, block_id)
      value_sorted <- sort(value_vec)
      n_all <- length(value_sorted)
      nb_blocks <- uniqueN(dtval$block_id)
      quantile_vec_all <- seq(1/(n_all)*nb_blocks/2, 1-1/(n_all)*nb_blocks/2, by=1/(n_all)*nb_blocks)
      quantile_idx_all <- round(quantile_vec_all * n_all)
      quantile_idx_all <- pmax(1, pmin(n_all, quantile_idx_all))
      value_obs <- value_sorted[quantile_idx_all]
      n_sample <- ncol(array_sample)

      idx <- array_sample[, 1]
      sample_vals <- dtval[J(idx), value]
      # sample_sorted <- sort(sample_vals)
      n_samp <- length(sample_vals)
      quantile_vec_samp <- seq(1/(n_samp)*length(idx)/2, 1-1/(n_samp)*length(idx)/2, by=1/(n_samp)*length(idx))
      quantile_idx_samp0 <- quantile_vec_samp * n_samp
      quantile_idx_samp <- round(quantile_idx_samp0 + 1e-10) ## rounding 5 at even digit
      quantile_idx_samp <- pmax(1, pmin(n_samp, quantile_idx_samp))

      upper_idx <- which(quantile_vec_samp > upper_quantile_boundary)


      results <- vector("list", n_sample)
      for (ii in seq_len(n_sample)) {
        idx <- array_sample[, ii]
        sample_vals <- dtval[J(idx), value]
        sample_sorted <- sort(sample_vals)
        # n_samp <- length(sample_sorted)

        value_estim <- sample_sorted[quantile_idx_samp]
        diff_vec <- value_obs - value_estim
        abs_diff_vec <- abs(diff_vec)
        results[[ii]] <- list(
          mae_q = mean(abs_diff_vec),
          erreur_q99 = quantile(abs_diff_vec, 0.99),
          mae_pointe = mean(abs(diff_vec[upper_idx]))
        )
      }

      dt_out <- rbindlist(results)
      setnames(dt_out, c("mae_q", "erreur_q99", "mae_pointe"),
               paste0(c("mae_q_", "erreur_q99_", "mae_pointe_"), ccol))
      dt_out
    })
    dt_score_quantile <- do.call(cbind, outlist)
    dt_score_quantile[, name_sample := colnames(array_sample)]
  } else {
    dt_score_quantile <- data.table()
  }
  dt_score_quantile[]
}
