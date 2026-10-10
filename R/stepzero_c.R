# R/stepzero_c.R -- functions for script 09c (ridge-stabilised fit + cluster-block-size comparison)
# Needs: stepzero.R, stepzero_b.R

PRIOR_SD       <- 2.5    # weak Gaussian prior (ridge) on standardised slopes; intercept unpenalised
GATE_FAIL_MAX  <- 0.10   # NEW gate: share of replicates with no usable fit (non-convergence / separation)

# Newton-Raphson for the down-weighted Poisson point-process likelihood with a ridge penalty.
fit_ppm_ridge <- function(Xp, Xq, aq, prior_sd = PRIOR_SD, eps = SIM_EPS, maxit = 80L, tol = 1e-10) {
  p <- ncol(Xq)
  pen <- c(0, rep(1 / prior_sd^2, p - 1L))
  b <- rep(0, p); b[1] <- log(nrow(Xp) / sum(aq))
  ll <- function(b) {
    ep <- drop(Xp %*% b); eq <- drop(Xq %*% b)
    v <- sum(ep) - eps * sum(exp(ep)) - sum(aq * exp(eq)) - 0.5 * sum(pen * b^2)
    if (is.finite(v)) v else -Inf
  }
  cur <- ll(b); conv <- FALSE; it <- 0L
  for (it in seq_len(maxit)) {
    ep <- exp(drop(Xp %*% b)); eq <- exp(drop(Xq %*% b)); wq <- aq * eq
    g <- colSums(Xp) - eps * colSums(Xp * ep) - colSums(Xq * wq) - pen * b
    H <- -(eps * crossprod(Xp * ep, Xp) + crossprod(Xq * wq, Xq)) - diag(pen, p)
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1; new <- -Inf
    for (h in 1:30) { new <- ll(b + s * step); if (new >= cur) break; s <- s / 2 }
    if (!is.finite(new) || new < cur) break
    b <- b + s * step
    d <- abs(new - cur); cur <- new
    if (d < tol * (1 + abs(cur)) && max(abs(s * step)) < 1e-6) { conv <- TRUE; break }
  }
  ep <- exp(drop(Xp %*% b)); eq <- exp(drop(Xq %*% b)); wq <- aq * eq
  H <- -(eps * crossprod(Xp * ep, Xp) + crossprod(Xq * wq, Xq)) - diag(pen, p)
  V <- tryCatch(solve(-H), error = function(e) matrix(NA_real_, p, p))
  list(coef = b, se = sqrt(diag(V)), vcov = V, converged = conv, iter = it)
}

# one ridge fit, one cluster bootstrap per block size (clusters list: name -> cluster id per cell)
fit_boot_multi <- function(cnt, Xq, a, cl_list, ik, B = 300L, beta_true_iqr = 0) {
  idx <- rep(seq_along(cnt), cnt)
  Xp <- Xq[idx, , drop = FALSE]
  f <- fit_ppm_ridge(Xp, Xq, a)
  j <- ncol(Xq)
  b <- f$coef[j] * ik
  bad <- !f$converged || !is.finite(f$se[j]) || any(abs(f$coef) > 8)
  do.call(rbind, lapply(names(cl_list), function(nm) {
    if (bad) return(data.frame(block_km = as.numeric(nm), n = length(idx), beta_iqr = b, se_model = NA_real_, se_boot = NA_real_,
                               lo = NA_real_, hi = NA_real_, failed = TRUE, covers = NA, rej_ci = NA))
    bs <- boot_onestep(f, Xp, cl_list[[nm]][idx], B) * ik
    q <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.frame(block_km = as.numeric(nm), n = length(idx), beta_iqr = b, se_model = f$se[j] * ik,
               se_boot = stats::sd(bs, na.rm = TRUE), lo = q[1], hi = q[2], failed = FALSE,
               covers = (q[1] <= beta_true_iqr) & (beta_true_iqr <= q[2]), rej_ci = (q[1] > 0) | (q[2] < 0))
  }))
}

run_power_c <- function(X, zs, a, xy, ks, ns, ratios, sigmas, blocks_km, reps, range_m, B = 300L,
                        truth = SIM_TRUTH_C, progress = 5L) {
  cl_list <- stats::setNames(lapply(blocks_km, function(b) cluster_ids(xy, b * 1000)), blocks_km)
  Xk <- lapply(ks, function(k) cbind(int = 1, X, z = zs[[as.character(k)]]$z)); names(Xk) <- ks
  ik <- lapply(ks, function(k) zs[[as.character(k)]]$iqr_per_sd); names(ik) <- ks
  eta_nuis <- drop(X %*% truth[colnames(X)])
  out <- list(); ii <- 0L
  for (r in seq_len(reps)) {
    for (sg in sigmas) {
      f <- if (sg > 0) rff_field(sg, range_m)(xy) else 0
      for (n_t in ns) {
        cnt0 <- simulate_counts(eta_nuis + f, a, n_t)
        if (sum(cnt0) >= 8L) for (k in ks) {
          row <- fit_boot_multi(cnt0, Xk[[as.character(k)]], a, cl_list, ik[[as.character(k)]], B, 0)
          ii <- ii + 1L; out[[ii]] <- cbind(row, rep = r, sigma = sg, n_target = n_t, k = k, ratio = 1)
        }
        for (k in ks) for (ra in ratios[ratios > 1]) {
          bt <- log(ra) / ik[[as.character(k)]]
          cnt <- simulate_counts(eta_nuis + f + bt * zs[[as.character(k)]]$z, a, n_t)
          if (sum(cnt) < 8L) next
          row <- fit_boot_multi(cnt, Xk[[as.character(k)]], a, cl_list, ik[[as.character(k)]], B, log(ra))
          ii <- ii + 1L; out[[ii]] <- cbind(row, rep = r, sigma = sg, n_target = n_t, k = k, ratio = ra)
        }
      }
    }
    if (r %% progress == 0L) message("  replicate ", r, " / ", reps)
  }
  do.call(rbind, out)
}

# power counts failed fits as non-rejections (unconditional); calibration uses successful fits only
summarise_power_c <- function(res) {
  g <- split(res, list(res$sigma, res$n_target, res$k, res$ratio, res$block_km), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    e <- d[!d$failed, ]
    data.frame(sigma = d$sigma[1], n_target = d$n_target[1], k = d$k[1], ratio = d$ratio[1], block_km = d$block_km[1],
               n_reps = nrow(d), share_failed = round(mean(d$failed), 3),
               bias_beta = round(mean(e$beta_iqr) - log(d$ratio[1]), 3),
               sd_beta = round(stats::sd(e$beta_iqr), 3),
               se_model_med = round(stats::median(e$se_model), 3), se_boot_med = round(stats::median(e$se_boot), 3),
               se_boot_over_sd = round(stats::median(e$se_boot) / stats::sd(e$beta_iqr), 2),
               rej_uncond = round(sum(e$rej_ci) / nrow(d), 3),
               coverage = round(mean(e$covers), 3))
  }))
}

mde_table_c <- function(ps, target = GATE_POWER) {
  a <- ps[ps$ratio > 1, ]
  g <- split(a, list(a$sigma, a$n_target, a$k, a$block_km), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    d <- d[order(d$ratio), ]
    pw <- d$rej_uncond; lr <- log(d$ratio)
    mde <- if (all(pw < target)) NA_real_ else if (pw[1] >= target) d$ratio[1] else {
      i <- which(pw >= target)[1]
      exp(lr[i - 1] + (target - pw[i - 1]) * (lr[i] - lr[i - 1]) / (pw[i] - pw[i - 1]))
    }
    data.frame(sigma = d$sigma[1], n_target = d$n_target[1], k = d$k[1], block_km = d$block_km[1],
               mde_ratio = round(mde, 2), note = ifelse(is.na(mde), paste0(">", max(d$ratio)), ""))
  }))
}

evaluate_gates_c <- function(ps, mde) {
  null <- ps[ps$ratio == 1, ]; alt <- ps[ps$ratio > 1, ]
  keys <- unique(ps[, c("sigma", "n_target", "block_km")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    sg <- keys$sigma[i]; nt <- keys$n_target[i]; bk <- keys$block_km[i]
    sel <- function(d) d$sigma == sg & d$n_target == nt & d$block_km == bk
    t1 <- mean(null$rej_uncond[sel(null)])
    cv <- mean(c(alt$coverage[sel(alt)], null$coverage[sel(null)]))
    fl <- mean(ps$share_failed[sel(ps)])
    m  <- mde$mde_ratio[mde$sigma == sg & mde$n_target == nt & mde$block_km == bk]
    data.frame(sigma = sg, n_target = nt, block_km = bk,
               type1 = round(t1, 3), gate_type1 = t1 <= GATE_TYPE1_MAX,
               coverage = round(cv, 3), gate_coverage = cv >= GATE_COVER_MIN,
               share_failed = round(fl, 3), gate_fail = fl <= GATE_FAIL_MAX,
               worst_mde = if (any(is.na(m))) NA_real_ else max(m),
               gate_mde = if (nt >= 90) (!any(is.na(m)) && all(m <= GATE_MDE_MAX)) else NA)
  }))
}