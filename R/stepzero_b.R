# R/stepzero_b.R -- functions for script 09b (power, size and CI coverage with a cluster bootstrap)
# Needs: stepzero.R (fit_ppm, simulate_counts, rff_field, nuisance_design, z_size, SIM_NUIS_TRUTH)

# ---- decision gates (PROPOSED; fixed in docs/deviation_log.md BEFORE the real run) ---------------------
GATE_TYPE1_MAX   <- 0.07   # null rejection rate, cluster-bootstrap percentile CI, pooled over k, per (sigma, n)
GATE_COVER_MIN   <- 0.90   # coverage of the 95% CI, pooled over k and effect sizes, per (sigma, n)
GATE_MDE_MAX     <- 2.0    # minimum detectable intensity ratio per IQR (80% power) required at n = 93
GATE_POWER       <- 0.80

# alluvial collapsed into pediplain (expected ~3 of 93 sites: not estimable; quasi-separation in S09a)
nuisance_design_c <- function(dom) {
  X <- nuisance_design(dom)
  X[, setdiff(colnames(X), "gm_alluvial"), drop = FALSE]
}
SIM_TRUTH_C <- SIM_NUIS_TRUTH[setdiff(names(SIM_NUIS_TRUTH), "gm_alluvial")]

cluster_ids <- function(xy, block_m = 20000) {
  paste(floor(xy[, 1] / block_m), floor(xy[, 2] / block_m))
}

# One-step cluster bootstrap. At the MLE the quadrature part of the score equals the presence part,
# so a resample changes the score by sum_c (m_c - 1) S_c  (S_c = cluster sums of presence covariates)
# and one Newton step gives beta* = beta_hat + V g*. No refit and no quadrature pass: microseconds.
boot_onestep <- function(f, Xp, cl, B = 300L) {
  p <- ncol(Xp)
  S <- rowsum(Xp, cl)
  nc <- nrow(S)
  if (nc < 2L) return(rep(NA_real_, B))
  M <- stats::rmultinom(B, nc, rep(1 / nc, nc))
  G <- crossprod(M - 1, S)
  f$coef[p] + drop(G %*% f$vcov[, p])
}

# Full-refit cluster bootstrap (validation only; slow)
boot_full <- function(Xq, aq, Xp, cl, B = 60L) {
  ucl <- unique(cl); nc <- length(ucl)
  rows <- split(seq_len(nrow(Xp)), factor(cl, levels = ucl))
  vapply(seq_len(B), function(b) {
    pick <- sample.int(nc, nc, TRUE)
    idx <- unlist(rows[pick], use.names = FALSE)
    fit_ppm(Xp[idx, , drop = FALSE], Xq, aq)$coef[ncol(Xq)]
  }, numeric(1))
}

# one fit + bootstrap for one pattern and one k; returns one row (effects on per-IQR log scale)
fit_boot_one <- function(cnt, Xq, a, cl_cell, ik, B = 300L, beta_true_iqr = 0) {
  idx <- rep(seq_along(cnt), cnt)
  Xp <- Xq[idx, , drop = FALSE]
  f <- fit_ppm(Xp, Xq, a)
  j <- ncol(Xq)
  b <- f$coef[j] * ik
  sep <- !f$converged && any(abs(f$coef) > 8)
  if (!f$converged || !is.finite(f$se[j])) {
    return(data.frame(n = length(idx), beta_iqr = b, se_model = NA_real_, se_boot = NA_real_, lo = NA_real_, hi = NA_real_,
                      converged = f$converged, separation = sep, covers = NA, rej_z = NA, rej_ci = NA))
  }
  bs <- boot_onestep(f, Xp, cl_cell[idx], B) * ik
  q <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  se_b <- stats::sd(bs, na.rm = TRUE)
  data.frame(n = length(idx), beta_iqr = b, se_model = f$se[j] * ik, se_boot = se_b, lo = q[1], hi = q[2],
             converged = TRUE, separation = FALSE,
             covers = (q[1] <= beta_true_iqr) & (beta_true_iqr <= q[2]),
             rej_z = 2 * stats::pnorm(-abs(b / se_b)) < 0.05,
             rej_ci = (q[1] > 0) | (q[2] < 0))
}

run_power <- function(X, zs, a, xy, ks, ns, ratios, sigmas, reps, range_m, B = 300L,
                      truth = SIM_TRUTH_C, progress = 5L) {
  cl_cell <- cluster_ids(xy)
  Xk <- lapply(ks, function(k) cbind(int = 1, X, z = zs[[as.character(k)]]$z)); names(Xk) <- ks
  ik <- lapply(ks, function(k) zs[[as.character(k)]]$iqr_per_sd); names(ik) <- ks
  eta_nuis <- drop(X %*% truth[colnames(X)])
  out <- list(); ii <- 0L
  for (r in seq_len(reps)) {
    for (sg in sigmas) {
      f <- if (sg > 0) rff_field(sg, range_m)(xy) else 0
      for (n_t in ns) {
        # null pattern: shared by all k (truth has no exposure effect)
        cnt0 <- simulate_counts(eta_nuis + f, a, n_t)
        if (sum(cnt0) >= 8L) for (k in ks) {
          row <- fit_boot_one(cnt0, Xk[[as.character(k)]], a, cl_cell, ik[[as.character(k)]], B, 0)
          ii <- ii + 1L; out[[ii]] <- cbind(row, rep = r, sigma = sg, n_target = n_t, k = k, ratio = 1)
        }
        for (k in ks) for (ra in ratios[ratios > 1]) {
          bt <- log(ra) / ik[[as.character(k)]]
          cnt <- simulate_counts(eta_nuis + f + bt * zs[[as.character(k)]]$z, a, n_t)
          if (sum(cnt) < 8L) next
          row <- fit_boot_one(cnt, Xk[[as.character(k)]], a, cl_cell, ik[[as.character(k)]], B, log(ra))
          ii <- ii + 1L; out[[ii]] <- cbind(row, rep = r, sigma = sg, n_target = n_t, k = k, ratio = ra)
        }
      }
    }
    if (r %% progress == 0L) message("  replicate ", r, " / ", reps)
  }
  do.call(rbind, out)
}

summarise_power <- function(res) {
  g <- split(res, list(res$sigma, res$n_target, res$k, res$ratio), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    ok <- d$converged & is.finite(d$se_boot)
    e <- d[ok, ]
    data.frame(sigma = d$sigma[1], n_target = d$n_target[1], k = d$k[1], ratio = d$ratio[1],
               n_reps = nrow(d), share_failed = round(1 - nrow(e) / nrow(d), 3),
               mean_n_sites = round(mean(e$n), 1),
               bias_beta = round(mean(e$beta_iqr) - log(d$ratio[1]), 3),
               mean_se_model = round(mean(e$se_model), 3), mean_se_boot = round(mean(e$se_boot), 3),
               sd_beta = round(stats::sd(e$beta_iqr), 3),
               rej_z = round(mean(e$rej_z), 3), rej_ci = round(mean(e$rej_ci), 3),
               coverage = round(mean(e$covers), 3))
  }))
}

# smallest ratio reaching target power (linear interpolation of power on log ratio)
mde_table <- function(ps, target = GATE_POWER, col = "rej_ci") {
  g <- split(ps[ps$ratio > 1, ], list(ps$sigma[ps$ratio > 1], ps$n_target[ps$ratio > 1], ps$k[ps$ratio > 1]), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    d <- d[order(d$ratio), ]
    pw <- d[[col]]; lr <- log(d$ratio)
    mde <- if (all(pw < target)) NA_real_ else if (pw[1] >= target) d$ratio[1] else {
      i <- which(pw >= target)[1]
      exp(lr[i - 1] + (target - pw[i - 1]) * (lr[i] - lr[i - 1]) / (pw[i] - pw[i - 1]))
    }
    data.frame(sigma = d$sigma[1], n_target = d$n_target[1], k = d$k[1],
               mde_ratio = round(mde, 2),
               note = ifelse(is.na(mde), paste0(">", max(d$ratio)), ifelse(pw[1] >= target, paste0("<=", d$ratio[1]), "")))
  }))
}

evaluate_gates <- function(ps, mde) {
  null <- ps[ps$ratio == 1, ]; alt <- ps[ps$ratio > 1, ]
  keys <- unique(ps[, c("sigma", "n_target")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    sg <- keys$sigma[i]; nt <- keys$n_target[i]
    t1 <- mean(null$rej_ci[null$sigma == sg & null$n_target == nt])
    cv <- mean(c(alt$coverage[alt$sigma == sg & alt$n_target == nt], null$coverage[null$sigma == sg & null$n_target == nt]))
    m  <- mde$mde_ratio[mde$sigma == sg & mde$n_target == nt]
    data.frame(sigma = sg, n_target = nt, type1_pooled = round(t1, 3), gate_type1 = t1 <= GATE_TYPE1_MAX,
               coverage_pooled = round(cv, 3), gate_coverage = cv >= GATE_COVER_MIN,
               worst_mde = if (any(is.na(m))) NA_real_ else max(m),
               gate_mde = if (nt >= 90) (!any(is.na(m)) && all(m <= GATE_MDE_MAX)) else NA)
  }))
}