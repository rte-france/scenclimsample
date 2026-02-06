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



#' build_l_block_cdf
#'
#' Compute blockwise empirical CDFs for specified columns in a data.table.
#'
#' @param dt_series data.table with a 'block_id' column and numeric columns.
#' @param col_vec character vector of column names to compute CDFs for. If NULL, all except 'block_id'.
#'
#' @return named list of matrices where each matrix is blockwise CDFs (rows: quantile grid, cols: blocks).
#' @examples
#' dt <- data.table(block_id = rep(1:3, each = 10), val1 = rnorm(30), val2 = rnorm(30))
#' build_l_block_cdf(dt, col_vec = c("val1", "val2"))
#' @export
#' @import data.table
#'
build_l_block_cdf <- function(dt_series, col_vec = NULL){

  if (is.null(col_vec)){
    col_vec <- grep("^block_id$", names(dt_series), value = TRUE, invert = TRUE)
  }

  n_blocks <- uniqueN(dt_series$block_id)
  qvec <- get_quantile_idx_all(n_all = nrow(dt_series), n_blocks)/nrow(dt_series)

  l_out <- purrr::map(col_vec, function(ccol){

    dt_tmp <- copy(dt_series[,.SD, .SDcols = c("block_id", ccol)])
    dt_tmp[,idx := seq(.N), by = .(block_id)]

    X <- dcast(dt_tmp, idx ~ block_id, value.var = ccol)
    X[["idx"]] <- NULL
    X <- as.matrix(X)
    q_grid <- quantile(dt_series[[ccol]], probs = qvec)

    # CDF values
    block_cdf <- matrix(0, nrow = length(qvec), ncol = n_blocks)

    for (b in seq_len(n_blocks)) {
      block_cdf[, b] <- colMeans(outer(X[, b], q_grid, `<=`))
    }
    block_cdf
  })

  names(l_out) <- col_vec
  l_out

}

#' get_dt_score_CDF_parallel
#'
#' Compute three CDF-based scores for each sample subset and each variable (mae_cdf, erreur_ks99, mae_cdfpointe).
#'
#' @param array_sample Integer matrix. Each column represents a subset/sample (block indices).
#' @param l_block_cdf Named list of blockwise CDF matrices (output from build_l_block_cdf).
#' @param upper_quantile_boundary Numeric. Quantile cutoff for point-wise CDF score (default: 0.998).
#'
#' @return data.table of CDF scores for all samples and all specified variables.
#' @examples
#' l_cdf <- build_l_block_cdf(dt, col_vec = c("val1"))
#' arr <- matrix(1:3, nrow = 3, ncol = 1)
#' get_dt_score_CDF_parallel(arr, l_cdf)
#' @export
#' @import data.table
#'
get_dt_score_CDF_parallel <- function(array_sample, l_block_cdf, upper_quantile_boundary = 0.998){

  n_sample <- ncol(array_sample)
  n_grid <- nrow(l_block_cdf[[1]])

  seq_grid <- seq(1/n_grid/2, 1-1/n_grid/2, by=1/n_grid)
  upper_idx <- which(seq_grid >= min(upper_quantile_boundary, max(seq_grid)))

  outlist <- furrr::future_imap(l_block_cdf, function(block_cdf, ccol) {

    full_cdf <- rowMeans(block_cdf)

    results <- vector("list", n_sample)
    for (ii in seq_len(n_sample)) {
      idx <- array_sample[, ii]
      subset_cdf <- rowMeans(block_cdf[, idx, drop = FALSE])
      abs_diff_vec <- abs(full_cdf - subset_cdf)


      results[[ii]] <- list(
        mae_cdf = mean(abs_diff_vec),
        erreur_ks99 = quantile(abs_diff_vec, 0.99),
        mae_cdfpointe = mean(abs_diff_vec[upper_idx])
      )

    }
    dt_out <- rbindlist(results)
    setnames(dt_out, c("mae_cdf", "erreur_ks99", "mae_cdfpointe"),
             paste0(c("mae_cdf_", "erreur_ks99_", "mae_cdfpointe_"), ccol))
    dt_out
  })

  names(outlist) <- NULL
  dt_score_cdf <- do.call(cbind, outlist)
  dt_score_cdf[, name_sample := colnames(array_sample)]
  dt_score_cdf[]

}


get_quantile_idx_all <- function(n_all, nb_blocks){

  quantile_vec_all <- seq(1/(n_all)*nb_blocks/2, 1-1/(n_all)*nb_blocks/2, by=1/(n_all)*nb_blocks)
  quantile_idx_all <- round(quantile_vec_all * n_all + 1e-10)
  pmax(1, pmin(n_all, quantile_idx_all))

}
