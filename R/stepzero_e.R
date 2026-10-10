# R/stepzero_e.R -- functions for script 09e (profile recovery: which size threshold k does the AIC profile pick?)
# Needs: stepzero.R, stepzero_b.R, stepzero_c.R

K_CLASS <- function(k) ifelse(k <= 5, "low", ifelse(k <= 7, "mid", "high"))

# ridge fit for every k on ONE pattern; AIC = -2 * Poisson-process loglik at the ridge estimate + 2p
fit_profile <- function(cnt, Xk, a, ik) {
  idx <- rep(seq_along(cnt), cnt)
  do.call(rbind, lapply(names(Xk), function(k) {
    Xq <- Xk[[k]]; Xp <- Xq[idx, , drop = FALSE]
    f <- fit_ppm_ridge(Xp, Xq, a)
    ll <- sum(drop(Xp %*% f$coef)) - sum(a * exp(drop(Xq %*% f$coef)))
    bad <- !f$converged || any(abs(f$coef) > 8) || !is.finite(ll)
    data.frame(k = as.integer(k), aic = if (bad) NA_real_ else -2 * ll + 2 * ncol(Xq),
               beta_iqr = f$coef[ncol(Xq)] * ik[[k]], failed = bad)
  }))
}

run_profile <- function(X, zs, a, xy, ks, k_trues, ratios, ns, sigma, range_m, reps,
                        truth = SIM_TRUTH_C, progress = 10L) {
  Xk <- lapply(ks, function(k) cbind(int = 1, X, z = zs[[as.character(k)]]$z)); names(Xk) <- ks
  ik <- lapply(ks, function(k) zs[[as.character(k)]]$iqr_per_sd); names(ik) <- ks
  eta_nuis <- drop(X %*% truth[colnames(X)])
  out <- list(); ii <- 0L
  for (r in seq_len(reps)) {
    f <- rff_field(sigma, range_m)(xy)
    for (n_t in ns) for (kt in k_trues) for (ra in ratios) {
      bt <- log(ra) / zs[[as.character(kt)]]$iqr_per_sd
      cnt <- simulate_counts(eta_nuis + f + bt * zs[[as.character(kt)]]$z, a, n_t)
      if (sum(cnt) < 8L) next
      pr <- fit_profile(cnt, Xk, a, ik)
      if (anyNA(pr$aic)) next                          # profile needs every k fitted
      ii <- ii + 1L
      out[[ii]] <- data.frame(rep = r, n_target = n_t, k_true = kt, ratio = ra,
                              k_hat = pr$k[which.min(pr$aic)],
                              daic_true = pr$aic[pr$k == kt] - min(pr$aic))
    }
    if (r %% progress == 0L) message("  replicate ", r, " / ", reps)
  }
  do.call(rbind, out)
}

summarise_profile <- function(res) {
  g <- split(res, list(res$n_target, res$k_true, res$ratio), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    data.frame(n_target = d$n_target[1], k_true = d$k_true[1], ratio = d$ratio[1], n_ok = nrow(d),
               p_exact = round(mean(d$k_hat == d$k_true), 3),
               p_within1 = round(mean(abs(d$k_hat - d$k_true) <= 1), 3),
               p_class = round(mean(K_CLASS(d$k_hat) == K_CLASS(d$k_true)), 3),
               mean_abs_err = round(mean(abs(d$k_hat - d$k_true)), 2),
               median_daic_true = round(stats::median(d$daic_true), 2),
               p90_daic_true = round(stats::quantile(d$daic_true, 0.9, names = FALSE), 2))
  }))
}

confusion_profile <- function(res, n_target, ratio, ks, k_trues) {
  d <- res[res$n_target == n_target & res$ratio == ratio, ]
  as.data.frame.matrix(table(k_true = factor(d$k_true, levels = k_trues), k_hat = factor(d$k_hat, levels = ks)))
}