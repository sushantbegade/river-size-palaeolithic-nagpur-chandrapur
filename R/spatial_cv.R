# R/spatial_cv.R -- functions for script 14 (frozen design: blocked spatial cross-validation of the point-process models)
# Needs: model_primary.R (fit_real, design_from_domain, exposure_domain), diagnostics.R (exposure_var), stepzero_b.R (cluster_ids)

# design matrices (domain cells and sites) for every model compared
cv_design <- function(dom, sub, des, ex, exl, ks) {
  Xn_dom <- des$dom; Xn_sub <- des$apply(sub)
  q <- list(M0_intercept = cbind(int = 1, rep(0, nrow(dom)))[, 1, drop = FALSE],
            M1_nuisance = cbind(int = 1, Xn_dom))
  p <- list(M0_intercept = cbind(int = rep(1, nrow(sub))),
            M1_nuisance = cbind(int = 1, Xn_sub))
  for (k in ks) {
    nm <- paste0("M2_river_k", k)
    q[[nm]] <- cbind(int = 1, Xn_dom, z = ex[[as.character(k)]]$z_dom)
    p[[nm]] <- cbind(int = 1, Xn_sub, z = ex[[as.character(k)]]$apply(sub))
  }
  q[["M3_lineament"]] <- cbind(int = 1, Xn_dom, z = exl$z_dom)
  p[["M3_lineament"]] <- cbind(int = 1, Xn_sub, z = exl$apply(sub))
  list(q = q, p = p)
}

# leave-one-block-out CV. Training uses ONLY cells and sites outside the held-out block (no leakage);
# the held-out score is the Poisson-process log-likelihood on the held-out block: sum log lambda(sites) - integral of lambda.
cv_run <- function(Xq_models, Xp_models, a, cell_blk, site_blk) {
  folds <- sort(unique(cell_blk))
  out <- vector("list", length(folds) * length(Xq_models)); ii <- 0L
  for (f in folds) {
    tc <- cell_blk != f; ts <- site_blk != f
    for (m in names(Xq_models)) {
      Xq <- Xq_models[[m]]; Xp <- Xp_models[[m]]
      fit <- fit_real(Xp[ts, , drop = FALSE], Xq[tc, , drop = FALSE], a[tc])
      ll <- NA_real_
      if (!fit$bad) {
        b <- fit$coef
        ll <- sum(drop(Xp[!ts, , drop = FALSE] %*% b)) - sum(a[!tc] * exp(drop(Xq[!tc, , drop = FALSE] %*% b)))
      }
      ii <- ii + 1L
      out[[ii]] <- data.frame(fold = f, model = m, n_test_sites = sum(!ts), n_test_cells = sum(!tc), loglik = ll, fit_ok = !fit$bad)
    }
  }
  do.call(rbind, out)
}

summarise_cv <- function(folds_df, sub_n) {
  w <- reshape(folds_df[, c("fold", "model", "loglik")], idvar = "fold", timevar = "model", direction = "wide")
  names(w) <- sub("^loglik\\.", "", names(w))
  ok <- stats::complete.cases(w)
  w <- w[ok, , drop = FALSE]
  models <- setdiff(names(w), "fold")
  base <- w[["M1_nuisance"]]
  do.call(rbind, lapply(models, function(m) {
    d <- w[[m]] - base
    data.frame(model = m, n_folds_used = nrow(w), n_folds_dropped = sum(!ok), cv_loglik = round(sum(w[[m]]), 2),
               delta_vs_nuisance = round(sum(d), 2), delta_per_site = round(sum(d) / sub_n, 3),
               folds_improved = sum(d > 1e-9), folds_worse = sum(d < -1e-9),
               se_delta_sum = round(stats::sd(d) * sqrt(nrow(w)), 2))
  }))
}