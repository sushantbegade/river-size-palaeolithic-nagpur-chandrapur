# R/numerics_check.R -- functions for script 22 (EXPLORATORY numerical checks, post-freeze): epsilon of the down-weighted likelihood,
# ridge vs maximum likelihood, and the origin of the bootstrap block grid. Uses the frozen helpers unchanged.
# Needs: stepzero.R (fit_ppm), stepzero_b.R (cluster_ids), stepzero_c.R (fit_ppm_ridge), model_primary.R (ci_row), diagnostics.R (boot_onestep_j).

# refit every k at several presence weights eps (frozen value 1e-6); also the unpenalised ML fit, and the share of the eps term in the objective
eps_check <- function(Xp_list, Xq_list, a, iqr_by_k, eps_vec = c(1e-4, 1e-6, 1e-8)) {
  do.call(rbind, lapply(names(Xp_list), function(k) {
    Xp <- Xp_list[[k]]; Xq <- Xq_list[[k]]; j <- ncol(Xq)
    ml <- fit_ppm(Xp, Xq, a, eps = 1e-6)
    do.call(rbind, lapply(eps_vec, function(e) {
      f <- fit_ppm_ridge(Xp, Xq, a, eps = e)
      data.frame(k = as.integer(k), eps = e, converged = f$converged, beta_iqr = f$coef[j] * iqr_by_k[[k]],
                 beta_iqr_ml = ml$coef[j] * iqr_by_k[[k]],
                 eps_term_share = e * sum(exp(drop(Xp %*% f$coef))) / sum(a * exp(drop(Xq %*% f$coef))))
    }))
  }))
}

# block ids on a grid whose origin is shifted by (ox_km, oy_km); (0, 0) reproduces cluster_ids()
shift_clusters <- function(xy, block_km, ox_km, oy_km) {
  cluster_ids(cbind(xy[, 1] - ox_km * 1000, xy[, 2] - oy_km * 1000), block_km * 1000)
}

# one-step cluster-bootstrap CI for the exposure slope with the block grid shifted
origin_row <- function(f, Xp, xy, block_km, ox_km, oy_km, iqr_per_sd, sd_log, B = 2000L) {
  cl <- shift_clusters(xy, block_km, ox_km, oy_km); j <- length(f$coef)
  ci <- ci_row(boot_onestep_j(f, Xp, cl, j, B), iqr_per_sd, sd_log, f$coef[j])
  data.frame(block_km = block_km, ox_km = ox_km, oy_km = oy_km, n_clusters = length(unique(cl)),
             max_cluster_share = round(max(table(cl)) / length(cl), 3), beta_iqr = ci$beta_iqr, lo = ci$lo, hi = ci$hi,
             excludes_zero = ci$lo > 0 | ci$hi < 0, n_boot_ok = ci$n_boot_ok)
}