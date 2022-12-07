## calcul des moyennes
get_dt_score_mean <- function(array_sel, dt_complet, dt_complet_mean){

  dt_sel_mean <- map_dfr(seq(ncol(array_sel)), function(ii){
    dt_complet[array_sel[, ii], colMeans(.SD), .SDcols = -"mcYear"]
  })
  dt_sel_mean$name_sel = colnames(array_sel)
  setDT(dt_sel_mean)

  ## score vs. 1000 années
  dt_score_mean <- abs(dt_complet_mean[rep(1, nrow(dt_sel_mean))] - dt_sel_mean[,.SD, .SDcols = -"name_sel"])
  dt_score_mean[,name_sel := dt_sel_mean$name_sel]
  dt_score_mean
}

## calcul des scores d'énergie
get_dt_score_energy <- function(array_sel, l_mat_A){

  dt_score_energy <- imap_dfc(l_mat_A, function(mat_A, mat_A_name){

    Nyear_total <- nrow(mat_A)
    w0 <- array(1/Nyear_total, c(Nyear_total,1))
    mat_A_w0 <- mat_A %*% w0
    w0_mat_A_w0_half <- 0.5*t(w0) %*% mat_A %*% w0

    vec_score_energy <- map_dbl(seq(ncol(array_sel)), function(ii){
      ws <- array(0, c(Nyear_total,1))
      ws[array_sel[, ii],1] <- 1/Nyear_total

      t(ws) %*% mat_A_w0 - w0_mat_A_w0_half - 0.5*t(ws) %*% mat_A %*% ws
    })

    data.table(vec_score_energy) %>% setnames(mat_A_name)
  }) %>%
    setDT()
  dt_score_energy[,name_sel := colnames(array_sel)]
}


## calcul des scores quantiles pour une colonne
get_dt_score_quantile_one_col <- function(array_sel, dt_complet_sorted, dt_complet_quantiles, quantile_pointe = 0.998){

  nrow_sel <- nrow(dt_complet_sorted[.(array_sel[,1]), on = "mcYear"])
  nb_year_sel <- nrow(array_sel)
  quantile_vec_allyears <- seq(1/nrow_sel*nb_year_sel, 1-1/nrow_sel*nb_year_sel, by = 1/nrow_sel*nb_year_sel)
  quantile_position_index <- round(quantile_vec_allyears*nrow_sel)
  quantile_position_pointe_index <- which(quantile_vec_allyears>quantile_pointe)

  dt_score_quantile <- map_dfr(seq(ncol(array_sel)), function(ii){

    value_estim <- sort(dt_complet_sorted[.(array_sel[,ii]), on = "mcYear"]$value)[quantile_position_index]
    diff_vec <- dt_complet_quantiles$value - value_estim

    data.table(mae_q = mean(abs(diff_vec)),
               erreur_q99 = quantile(abs(diff_vec), 0.99),
               mae_pointe = mean(abs(diff_vec[quantile_position_pointe_index])))
  })

}

## concatenation des scores quantiles pour toutes les colonnes
get_dt_score_quantile <- function(array_sel, l_dt_complet_sorted, l_dt_complet_quantiles, ...){

  dt_score_quantiles <- reduce(map(seq(length(l_dt_complet_sorted)), function(ccol){
    dt_score_quantiles_one_col <- get_dt_score_quantile_one_col(array_sel, l_dt_complet_sorted[[ccol]], l_dt_complet_quantiles[[ccol]], ...)
    names(dt_score_quantiles_one_col) <- paste0(names(dt_score_quantiles_one_col), "_", names(l_dt_complet_sorted)[ccol])
    dt_score_quantiles_one_col
  }),
  cbind)
  dt_score_quantiles[,name_sel := colnames(array_sel)]

}
