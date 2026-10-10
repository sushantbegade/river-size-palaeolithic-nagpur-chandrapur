# R/robustness.R -- functions for script 16 (pre-registered robustness matrix, plus two labelled additions)
# Needs: model_primary.R, diagnostics.R (subset_profile, exposure_var), stepzero_c.R (PRIOR_SD), quadrature.R (nearest_dist, SIZE_K)

# random greedy thinning: keep a random subset in which no two points are closer than d_m (returns row indices)
thin_greedy <- function(xy, d_m) {
  n <- nrow(xy); keep <- integer(0)
  for (i in sample.int(n)) {
    if (!length(keep) || all(sqrt((xy[keep, 1] - xy[i, 1])^2 + (xy[keep, 2] - xy[i, 2])^2) >= d_m)) keep <- c(keep, i)
  }
  sort(keep)
}

# repeated thinning: distribution of the k-th threshold coefficient (per IQR) over R random thinnings
thin_profile <- function(label, sub, d_m, k, des, ex, Xq, a, R = 100L) {
  Xs <- des$apply(sub); z <- ex$apply(sub)
  xy <- as.matrix(sub[, c("x", "y")])
  res <- lapply(seq_len(R), function(r) {
    i <- thin_greedy(xy, d_m)
    f <- fit_real(cbind(int = 1, Xs[i, , drop = FALSE], z = z[i]), Xq, a)
    data.frame(n_kept = length(i), beta_iqr = if (f$bad) NA_real_ else f$coef[length(f$coef)] * ex$iqr_per_sd, ok = !f$bad)
  })
  r <- do.call(rbind, res); g <- r[r$ok, ]
  data.frame(subset = label, thin_km = d_m / 1000, k = k, replicates = R, usable = nrow(g), mean_n_kept = round(mean(r$n_kept), 1),
             median_beta_iqr = round(stats::median(g$beta_iqr), 3), q025 = round(stats::quantile(g$beta_iqr, 0.025, names = FALSE), 3),
             q975 = round(stats::quantile(g$beta_iqr, 0.975, names = FALSE), 3), share_positive = round(mean(g$beta_iqr > 0), 3))
}

# build everything one variant needs and run the profile for the given sample sets
run_variant <- function(label, dom_v, sets_v, ks, exposure_fun, drop_cols = "gm_alluvial", block_km = 40, B = 2000L) {
  des <- design_from_domain(dom_v, drop = drop_cols)
  ex <- stats::setNames(lapply(ks, exposure_fun), ks)
  Xq <- lapply(ks, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)); names(Xq) <- ks
  a <- dom_v$area_km2
  do.call(rbind, lapply(names(sets_v), function(nm) {
    r <- subset_profile(nm, sets_v[[nm]], ks, des, ex, Xq, a, block_km, B)
    cbind(variant = label, r)
  }))
}