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
