# R/model_primary.R -- functions for script 11 (primary model on the real sites).
# Uses the FROZEN estimator and bootstrap helpers unchanged: fit_ppm_ridge (stepzero_c.R), boot_onestep, cluster_ids (stepzero_b.R).
# Needs: config.R, quadrature.R (log_d), stepzero.R (nuisance_design), stepzero_b.R, stepzero_c.R

# raw nuisance covariates for any table with the grid's column names (must match nuisance_design exactly)
nuisance_raw <- function(d) {
  cbind(elev = d$elev_m,
        slope = log(d$slope_deg + 0.1),
        gm_plateau  = as.numeric(d$gm4 == "plateau_hills"),
        gm_alluvial = as.numeric(d$gm4 == "alluvial"),
        crop  = clamp01(d$frac_crop),
        built = sqrt(clamp01(d$frac_built)),
        road  = log_d(d$d_road_major),
        settle = log_d(d$d_settle))
}

# standardise raw covariates with the ANALYSIS-DOMAIN centre and scale (sites are never standardised on their own)
design_from_domain <- function(dom, drop = "gm_alluvial") {
  Xd <- nuisance_design(dom)
  cen <- attr(Xd, "scaled:center"); scl <- attr(Xd, "scaled:scale")
  keep <- setdiff(colnames(Xd), drop)
  dom_std <- Xd[, keep, drop = FALSE]
  list(dom = dom_std, center = cen[keep], scale = scl[keep],
       apply = function(d) {
         R <- nuisance_raw(d)[, keep, drop = FALSE]
         sweep(sweep(R, 2, cen[keep], "-"), 2, scl[keep], "/")
       })
}

# exposure for threshold k: standardised with the domain mean/sd of log distance
exposure_domain <- function(dom, k) {
  v <- log_d(dom[[paste0("d_k", k)]])
  m <- mean(v); s <- stats::sd(v)
  list(z_dom = (v - m) / s, sd_log = s, iqr_per_sd = stats::IQR(v) / s,
       apply = function(d) (log_d(d[[paste0("d_k", k)]]) - m) / s)
}

# ridge fit + Poisson-process AIC on real presence rows
fit_real <- function(Xp, Xq, a) {
  f <- fit_ppm_ridge(Xp, Xq, a)
  ll <- sum(drop(Xp %*% f$coef)) - sum(a * exp(drop(Xq %*% f$coef)))
  # separation heuristic applies to the slopes only: with few sites the intercept (log intensity per km2) is legitimately below -8
  bad <- !f$converged || any(abs(f$coef[-1]) > 8) || !is.finite(ll)
  c(f, list(loglik = ll, aic = -2 * ll + 2 * ncol(Xq), bad = bad))
}

# cluster bootstrap with full refits for several k at once (same resample for every k)
boot_refit <- function(Xq_list, Xp_list, a, cl, B, progress = 50L) {
  ks <- names(Xq_list)
  ucl <- unique(cl); nc <- length(ucl)
  rows <- split(seq_len(nrow(Xp_list[[1]])), factor(cl, levels = ucl))
  beta <- matrix(NA_real_, B, length(ks), dimnames = list(NULL, ks)); aic <- beta
  for (b in seq_len(B)) {
    idx <- unlist(rows[sample.int(nc, nc, TRUE)], use.names = FALSE)
    for (k in ks) {
      f <- fit_real(Xp_list[[k]][idx, , drop = FALSE], Xq_list[[k]], a)
      if (!f$bad) { beta[b, k] <- f$coef[length(f$coef)]; aic[b, k] <- f$aic }
    }
    if (b %% progress == 0L) message("    bootstrap ", b, " / ", B)
  }
  list(beta = beta, aic = aic)
}

ci_row <- function(b_std, iqr_per_sd, sd_log, est_std) {
  q <- stats::quantile(b_std[is.finite(b_std)], c(0.025, 0.975), names = FALSE)
  e <- est_std * iqr_per_sd; lo <- q[1] * iqr_per_sd; hi <- q[2] * iqr_per_sd
  data.frame(beta_iqr = round(e, 3), lo = round(lo, 3), hi = round(hi, 3),
             ratio_q25_vs_q75 = round(exp(-e), 2), ratio_lo = round(exp(-hi), 2), ratio_hi = round(exp(-lo), 2),
             per_doubling = round(exp(-log(2) * est_std / sd_log), 2),
             per_doubling_lo = round(exp(-log(2) * q[2] / sd_log), 2), per_doubling_hi = round(exp(-log(2) * q[1] / sd_log), 2),
             n_boot_ok = sum(is.finite(b_std)))
}