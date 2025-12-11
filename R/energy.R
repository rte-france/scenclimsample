# Copyright (C) 2022-2023, RTE (http://www.rte-france.com)
# SPDX-License-Identifier: MPL-2.0

#' Compute block mean of a matrix
#'
#' @param m Matrix
#' @param blockszi row block size
#' @param blockszj col block size
#'
#' @return an averaged matrix
#'
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

#' Build averaged pairwise euclidean distance between blocks
#'
#' @param dt_series a data.table of series with block_id blocks
#'
#' @return a matrix of averaged pairwise euclidean distance between blocks
#'
#' @export
#' @import data.table
#'
build_euclid_dist_matrix <- function(dt_series){
  unique_blocks <- dt_series[,unique(block_id)]
  nblocks <- length(unique_blocks)
  block_length <- dt_series[,.N, by = block_id][,mean(N)]

  mat_A <- purrr::map(unique_blocks, function(ii_block){
    ii_block_below <- seq(ii_block)
    row_mat_A <- fields::rdist(x1 = dt_series[.(ii_block),.SD,.SDcols = -c("block_id"), on = "block_id"],
                               x2 = dt_series[.(ii_block_below),.SD,.SDcols = -c("block_id"), on = "block_id"])
    row_mat_A_mean <- c(matrixblockmean(row_mat_A, block_length, block_length), rep(0, nblocks - ii_block))

    row_mat_A_mean[ii_block] <- row_mat_A_mean[ii_block]/2
    row_mat_A_mean
  })

  mat_A <- do.call("rbind", args = mat_A)
  mat_A <- (mat_A + t(mat_A))

  mat_A
}

