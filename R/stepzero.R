# R/stepzero.R -- functions for script 09a (Step 0: engine + simulation; no real site is used)
# Needs: config.R (CRS_ANALYSIS), quadrature.R (log_d, SIZE_K)

SIM_EPS <- 1e-6   # presence weight of the down-weighted Poisson likelihood

# Assumed nuisance truth on the standardised scale (ASSUMPTIONS, not estimates; sensitivity axis in S09b)
SIM_NUIS_TRUTH <- c(elev = 0.10, slope = -0.10, gm_plateau = -0.40, gm_alluvial = 0.30,
                    crop = 0.30, built = -0.20, road = -0.20, settle = -0.10)

# focal means can leave tiny negative / >1 values (floating point): clamp fractions to [0, 1]
clamp01 <- function(x) pmin(1, pmax(0, x))

# standardised nuisance design (no intercept)
nuisance_design <- function(dom) {
  X <- cbind(elev = dom$elev_m,
             slope = log(dom$slope_deg + 0.1),
             gm_plateau  = as.numeric(dom$gm4 == "plateau_hills"),
             gm_alluvial = as.numeric(dom$gm4 == "alluvial"),
             crop  = clamp01(dom$frac_crop),
             built = sqrt(clamp01(dom$frac_built)),
             road  = log_d(dom$d_road_major),
             settle = log_d(dom$d_settle))
  sds <- apply(X, 2, stats::sd)
  if (any(!is.finite(sds) | sds == 0)) stop("Constant or NA nuisance column: ", paste(names(sds)[!is.finite(sds) | sds == 0], collapse = ", "), call. = FALSE)
  scale(X)
}

# standardised log-distance for size threshold k; iqr_per_sd converts SD-scale slope to per-IQR scale
z_size <- function(dom, k) {
  z <- log_d(dom[[paste0("d_k", k)]])
  s <- stats::sd(z)
  list(z = (z - mean(z)) / s, iqr_per_sd = stats::IQR(z) / s)
}

# Down-weighted Poisson point-process fit by Newton-Raphson with step halving.
# Presence rows Xp (weight eps, response 1/eps); quadrature rows Xq with weights aq (cell area, km2).
# Column 1 of both matrices must be the intercept.
fit_ppm <- function(Xp, Xq, aq, eps = SIM_EPS, maxit = 60L, tol = 1e-10) {
  p <- ncol(Xq)
  b <- rep(0, p); b[1] <- log(nrow(Xp) / sum(aq))
  ll <- function(b) {
    ep <- drop(Xp %*% b); eq <- drop(Xq %*% b)
    v <- sum(ep) - eps * sum(exp(ep)) - sum(aq * exp(eq))
    if (is.finite(v)) v else -Inf
  }
  cur <- ll(b); conv <- FALSE; it <- 0L
  for (it in seq_len(maxit)) {
    ep <- exp(drop(Xp %*% b)); eq <- exp(drop(Xq %*% b)); wq <- aq * eq
    g <- colSums(Xp) - eps * colSums(Xp * ep) - colSums(Xq * wq)
    H <- -(eps * crossprod(Xp * ep, Xp) + crossprod(Xq * wq, Xq))
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1; new <- -Inf
    for (h in 1:30) { new <- ll(b + s * step); if (new >= cur) break; s <- s / 2 }
    if (!is.finite(new) || new < cur) break
    b <- b + s * step
    d <- abs(new - cur); cur <- new
    if (d < tol * (1 + abs(cur)) && max(abs(s * step)) < 1e-6) { conv <- TRUE; break }
  }
  eq <- exp(drop(Xq %*% b)); ep <- exp(drop(Xp %*% b)); wq <- aq * eq
  H <- -(eps * crossprod(Xp * ep, Xp) + crossprod(Xq * wq, Xq))
  V <- tryCatch(solve(-H), error = function(e) matrix(NA_real_, p, p))
  ll_pp <- sum(drop(Xp %*% b)) - sum(wq)
  list(coef = b, se = sqrt(diag(V)), vcov = V, loglik_pp = ll_pp, aic = -2 * ll_pp + 2 * p,
       converged = conv, iter = it)
}

# random-Fourier-feature approximation of a Gaussian field with kernel sigma^2 * exp(-d / range_m)
# (Matern nu = 1/2; spectral density = bivariate t with 1 df). Returns a function(xy) -> field values.
rff_field <- function(sigma, range_m, M = 400L) {
  W <- cbind(stats::rnorm(M), stats::rnorm(M)) / sqrt(stats::rchisq(M, 1)) / range_m
  ph <- stats::runif(M, 0, 2 * pi)
  function(xy, chunk = 20000L) {
    n <- nrow(xy); out <- numeric(n)
    for (s in seq(1L, n, by = chunk)) {
      i <- s:min(n, s + chunk - 1L)
      out[i] <- rowSums(cos(sweep(xy[i, , drop = FALSE] %*% t(W), 2, ph, "+")))
    }
    sigma * sqrt(2 / M) * out
  }
}

# counts per cell with expected total n_target
simulate_counts <- function(eta, a, n_target) {
  lam <- exp(eta - max(eta))
  lam <- lam * n_target / sum(a * lam)
  stats::rpois(length(a), a * lam)
}

# one replicate: returns data.frame(scenario, k, ...) fitted for all k
fit_all_k <- function(counts, Xk, a, ik, eps = SIM_EPS) {
  idx <- rep(seq_along(counts), counts)
  do.call(rbind, lapply(names(Xk), function(k) {
    Xq <- Xk[[k]]
    f <- fit_ppm(Xq[idx, , drop = FALSE], Xq, a, eps)
    j <- ncol(Xq)
    b <- f$coef[j] * ik[[k]]; se <- f$se[j] * ik[[k]]
    data.frame(k = as.integer(k), n = length(idx), beta_iqr = b, se_iqr = se,
               z = b / se, p = 2 * stats::pnorm(-abs(b / se)), aic = f$aic, converged = f$converged,
               separation = !f$converged && any(abs(f$coef) > 8))
  }))
}

run_sim <- function(X, zs, a, xy, ks, sigmas, ratio, k_true, n_target, reps, range_m, truth = SIM_NUIS_TRUTH,
                    progress = 10L) {
  Xk <- lapply(ks, function(k) cbind(int = 1, X, z = zs[[as.character(k)]]$z)); names(Xk) <- ks
  ik <- lapply(ks, function(k) zs[[as.character(k)]]$iqr_per_sd); names(ik) <- ks
  eta_nuis <- drop(X %*% truth[colnames(X)])
  b_alt_std <- log(ratio) / zs[[as.character(k_true)]]$iqr_per_sd
  zt <- zs[[as.character(k_true)]]$z
  out <- vector("list", reps * length(sigmas))
  ii <- 0L
  for (r in seq_len(reps)) {
    for (sg in sigmas) {
      f <- if (sg > 0) rff_field(sg, range_m)(xy) else 0
      for (sc in c("null", "alt")) {
        eta <- eta_nuis + f + if (sc == "alt") b_alt_std * zt else 0
        cnt <- simulate_counts(eta, a, n_target)
        if (sum(cnt) < 10L) next
        res <- fit_all_k(cnt, Xk, a, ik)
        res$rep <- r; res$sigma <- sg; res$scenario <- sc
        ii <- ii + 1L; out[[ii]] <- res
      }
    }
    if (r %% progress == 0L) message("  replicate ", r, " / ", reps)
  }
  do.call(rbind, out[seq_len(ii)])
}

summarise_sim <- function(sim, alpha = 0.05, true_beta_iqr = NULL) {
  g <- split(sim, list(sim$scenario, sim$sigma, sim$k), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    ok <- d$converged & is.finite(d$se_iqr)
    d <- d[ok, ]
    data.frame(scenario = d$scenario[1], sigma = d$sigma[1], k = d$k[1], n_ok = nrow(d),
               mean_n_sites = round(mean(d$n), 1),
               mean_beta_iqr = round(mean(d$beta_iqr), 3), sd_beta_iqr = round(stats::sd(d$beta_iqr), 3),
               mean_se_iqr = round(mean(d$se_iqr), 3),
               se_ratio = round(mean(d$se_iqr) / stats::sd(d$beta_iqr), 2),
               reject_rate = round(mean(d$p < alpha), 3))
  }))
}

khat_table <- function(sim) {
  a <- sim[sim$scenario == "alt" & sim$converged, ]
  a$key <- paste(a$sigma, a$rep)
  g <- split(a, a$key)
  g <- Filter(function(d) nrow(d) == length(unique(sim$k)), g)
  kh <- do.call(rbind, lapply(g, function(d) data.frame(sigma = d$sigma[1], khat = d$k[which.min(d$aic)])))
  as.data.frame(table(sigma = kh$sigma, khat = kh$khat), responseName = "n_reps")
}

# empirical correlation of the approximate field against exp(-d/range)
field_check <- function(xy, sigma, range_m, nfield = 300L, npair = 3000L, M = 400L) {
  n <- nrow(xy)
  i <- sample.int(n, npair, TRUE); j <- sample.int(n, npair, TRUE)
  d <- sqrt((xy[i, 1] - xy[j, 1])^2 + (xy[i, 2] - xy[j, 2])^2)
  Fi <- matrix(0, nfield, npair); Fj <- Fi
  for (r in seq_len(nfield)) {
    f <- rff_field(sigma, range_m, M)
    Fi[r, ] <- f(xy[i, , drop = FALSE]); Fj[r, ] <- f(xy[j, , drop = FALSE])
  }
  cr <- vapply(seq_len(npair), function(q) stats::cor(Fi[, q], Fj[, q]), numeric(1))
  br <- c(0, 0.25, 0.5, 1, 2) * range_m
  br <- c(br, Inf)
  bin <- cut(d, br, include.lowest = TRUE)
  out <- do.call(rbind, lapply(levels(bin), function(l) {
    s <- bin == l
    data.frame(bin = l, n_pairs = sum(s), mean_dist_km = round(mean(d[s]) / 1000, 1),
               emp_corr = round(mean(cr[s], na.rm = TRUE), 3),
               target_corr = round(mean(exp(-d[s] / range_m)), 3))
  }))
  list(by_bin = out, marginal_sd = round(stats::sd(c(Fi)), 3), target_sd = sigma)
}