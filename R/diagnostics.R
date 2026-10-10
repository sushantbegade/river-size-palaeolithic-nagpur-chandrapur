# R/diagnostics.R -- functions for script 12 (pre-registered secondary analyses on the real sites)
# Uses frozen helpers unchanged. Needs: model_primary.R (fit_real, boot_refit, ci_row, exposure_domain), stepzero_b.R (cluster_ids)

# exposure for any distance column (e.g. d_lineament), standardised with the domain mean/sd of log distance
exposure_var <- function(dom, var) {
  v <- log_d(dom[[var]]); m <- mean(v); s <- stats::sd(v)
  list(z_dom = (v - m) / s, sd_log = s, iqr_per_sd = stats::IQR(v) / s,
       apply = function(d) (log_d(d[[var]]) - m) / s)
}

# one-step cluster bootstrap for coefficient j (generalises boot_onestep, which uses the last coefficient)
boot_onestep_j <- function(f, Xp, cl, j, B = 2000L) {
  S <- rowsum(Xp, cl); nc <- nrow(S)
  if (nc < 2L) return(rep(NA_real_, B))
  M <- stats::rmultinom(B, nc, rep(1 / nc, nc))
  G <- crossprod(M - 1, S)
  f$coef[j] + drop(G %*% f$vcov[, j])
}

# profile over k for one subset of sites: ridge fit, one-step cluster-bootstrap CI, AIC
subset_profile <- function(label, sub, ks, des, ex, Xq_list, a, block_km, B = 2000L) {
  cl <- cluster_ids(as.matrix(sub[, c("x", "y")]), block_km * 1000)
  Xs <- des$apply(sub)
  do.call(rbind, lapply(as.character(ks), function(k) {
    Xp <- cbind(int = 1, Xs, z = ex[[k]]$apply(sub))
    f <- fit_real(Xp, Xq_list[[k]], a)
    base <- data.frame(subset = label, n_sites = nrow(sub), n_clusters = length(unique(cl)), k = as.integer(k))
    if (f$bad) return(cbind(base, data.frame(beta_iqr = NA, lo = NA, hi = NA, ratio_q25_vs_q75 = NA, ratio_lo = NA, ratio_hi = NA,
                                             per_doubling = NA, per_doubling_lo = NA, per_doubling_hi = NA, n_boot_ok = 0L, aic = NA, fit_ok = FALSE)))
    j <- length(f$coef)
    b <- boot_onestep_j(f, Xp, cl, j, B)
    if (sum(is.finite(b)) < 50L) b <- rep(NA_real_, 2L)
    ci <- if (all(is.na(b))) data.frame(beta_iqr = round(f$coef[j] * ex[[k]]$iqr_per_sd, 3), lo = NA, hi = NA, ratio_q25_vs_q75 = round(exp(-f$coef[j] * ex[[k]]$iqr_per_sd), 2),
                                        ratio_lo = NA, ratio_hi = NA, per_doubling = NA, per_doubling_lo = NA, per_doubling_hi = NA, n_boot_ok = 0L)
    else ci_row(b, ex[[k]]$iqr_per_sd, ex[[k]]$sd_log, f$coef[j])
    cbind(base, ci, aic = round(f$aic, 2), fit_ok = TRUE)
  }))
}

# leave-one-group-out refits at one k (groups = references or 40 km blocks)
leave_one_out <- function(sub, group, des, ex_k, Xq, a) {
  Xs <- des$apply(sub); Xp_all <- cbind(int = 1, Xs, z = ex_k$apply(sub))
  g <- sub[[group]]
  full <- fit_real(Xp_all, Xq, a)
  out <- do.call(rbind, lapply(sort(unique(g)), function(u) {
    keep <- g != u
    f <- fit_real(Xp_all[keep, , drop = FALSE], Xq, a)
    data.frame(left_out = u, n_left_out = sum(!keep), n_remaining = sum(keep),
               beta_iqr = if (f$bad) NA_real_ else round(f$coef[length(f$coef)] * ex_k$iqr_per_sd, 3), fit_ok = !f$bad)
  }))
  out$change_vs_full <- round(out$beta_iqr - full$coef[length(full$coef)] * ex_k$iqr_per_sd, 3)
  out[order(-abs(out$change_vs_full)), ]
}