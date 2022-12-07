matrixblockmean <- function(m, blockszi, blockszj){
  out <- matrix(NA, nrow(m)/blockszi, ncol(m)/blockszj)
  for(i in 1:(nrow(m)/blockszi)){
    for(j in 1:(ncol(m)/blockszj)){
      iblock <- (1:blockszi) + blockszi*(i - 1)
      jblock <- (1:blockszj) + blockszj*(j - 1)
      out[i, j] <- mean(m[iblock, jblock])
    }
  }
  out
}

#' @export
build_euclid_dist_matrix <- function(dt_series, block_id){
  unique_blocks <- unique(dt_series[[block_id]])
  nblocks <- length(unique_blocks)
  block_length <- dt_series[,.N, by = block_id][,mean(N)]

  mat_A <- purrr::map(unique_blocks, function(ii_block)
    fields::rdist(x1 = dt_series[.(ii_block),.SD,.SDcols = -c(block_id), on = block_id],
                  x2 = dt_series[,.SD,.SDcols = -c(block_id)]) %>%
      matrixblockmean(block_length, block_length)) %>%
    do.call("rbind", args = .)

  mat_A
}

