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
get_dt_score_avg <- function(array_sample, dt_series, col_vec = NULL, ref_block_id = NULL){

  if (is.null(col_vec)){
    col_vec <- grep("^block_id$", names(dt_series), value = TRUE, invert = TRUE)
  }

  # Average by block and average over all blocks
  if (is.null(ref_block_id)){
    dt_avg_over_all_blocks <- dt_series[, lapply(.SD, mean), .SDcols = col_vec]
  } else {
    dt_avg_over_all_blocks <- dt_series[block_id %in% ref_block_id, lapply(.SD, mean), .SDcols = col_vec]
  }

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
get_dt_score_energy <- function(array_sample, l_mat_A, ref_block_id = NULL){

  # Loop over matrices A
  dt_score_energy <- furrr::future_imap_dfc(l_mat_A, function(mat_A, mat_A_name){

    Nyear_total <- nrow(mat_A)
    Nyear_sample <- nrow(array_sample)

    if(is.null(ref_block_id)){
      w0 <- array(1/Nyear_total, c(Nyear_total,1))
    } else {
      w0 <- array(0, c(nrow(l_mat_A[[1]]),1))
      w0[ref_block_id,1] <- 1/length(ref_block_id)
    }

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
get_dt_score_quantile <- function(array_sample, dt_series, col_vec = NULL, upper_quantile_boundary = 0.998, ref_block_id = NULL) {
  if (is.null(col_vec)) {
    col_vec <- setdiff(names(dt_series), "block_id")
  }

  if (is.null(ref_block_id)){
    ref_block_id <- sort(unique(dt_series$block_id))
  }

  if(ncol(array_sample)>0){

    outlist <- furrr::future_map(col_vec, function(ccol) {

      dtval <- dt_series[,.SD, .SDcols = c("block_id", ccol)]
      setnames(dtval, ccol, "value")
      setkey(dtval, block_id)

      dtval_ref <- dtval[block_id %in% ref_block_id]
      setkey(dtval_ref, block_id)

      value_sorted <- sort(dtval_ref$value)
      n_all <- length(value_sorted)
      nb_blocks <- uniqueN(dtval_ref$block_id)
      quantile_vec_all <- seq(1/(n_all)*nb_blocks/2, 1-1/(n_all)*nb_blocks/2, by=1/(n_all)*nb_blocks)
      quantile_idx_all <- round(quantile_vec_all * n_all)
      quantile_idx_all <- pmax(1, pmin(n_all, quantile_idx_all))
      value_obs <- value_sorted[quantile_idx_all]

      idx <- array_sample[, 1]
      length_sample <- nrow(dtval[J(idx)])
      quantile_vec_samp <- seq(1/(length_sample)*length(idx)/2, 1-1/(length_sample)*length(idx)/2, by=1/(length_sample)*length(idx))
      quantile_idx_samp0 <- quantile_vec_samp * length_sample
      quantile_idx_samp <- round(quantile_idx_samp0 + 1e-10) ## rounding 5 at even digit
      quantile_idx_samp <- pmax(1, pmin(length_sample, quantile_idx_samp))

      upper_idx <- which(quantile_vec_samp > upper_quantile_boundary)

      n_sample <- ncol(array_sample)
      results <- vector("list", n_sample)
      for (ii in seq_len(n_sample)) {
        idx <- array_sample[, ii]
        sample_vals <- dtval[J(idx), value]
        sample_sorted <- sort(sample_vals)

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
get_dt_score_CDF <- function(array_sample, l_block_cdf, upper_quantile_boundary = 0.998, ref_block_id = NULL){

  n_sample <- ncol(array_sample)
  n_grid <- nrow(l_block_cdf[[1]])

  seq_grid <- seq(1/n_grid/2, 1-1/n_grid/2, by=1/n_grid)
  upper_idx <- which(seq_grid >= min(upper_quantile_boundary, max(seq_grid)))

  outlist <- furrr::future_imap(l_block_cdf, function(block_cdf, ccol) {

    if (is.null(ref_block_id)){
      full_cdf <- rowMeans(block_cdf)
    } else {
      full_cdf <- rowMeans(block_cdf[,ref_block_id])
    }

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
