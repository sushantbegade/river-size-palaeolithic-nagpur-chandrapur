# R/survey_bias.R -- functions for script 15 (survey-effort sensitivity: reference-aware clusters and a reference-conditional likelihood)
# Needs: model_primary.R (fit_real, ci_row), diagnostics.R (boot_onestep_j), stepzero_b.R (cluster_ids), stepzero_c.R (PRIOR_SD)

# survey footprint per reference = hull of ALL its sites (any type; >= min_n sites) or the points themselves, buffered, clipped to the study
ref_footprints <- function(sites_sf, study, buf_m = 10000, min_n = 5L) {
  by <- split(seq_len(nrow(sites_sf)), sites_sf$ref_clean)
  fp <- lapply(by, function(ix) {
    u <- sf::st_union(sf::st_geometry(sites_sf)[ix])
    h <- if (length(ix) >= min_n) sf::st_convex_hull(u) else u
    sf::st_geometry(sf::st_intersection(sf::st_buffer(h, buf_m), study))[[1]]
  })
  stats::setNames(sf::st_sfc(fp, crs = sf::st_crs(sites_sf)), names(by))
}

# indices of domain cells inside each footprint
footprint_cells <- function(footprints, cell_xy, crs) {
  pts <- sf::st_as_sf(as.data.frame(cell_xy), coords = c(1, 2), crs = crs)
  hit <- sf::st_intersects(footprints, pts)
  stats::setNames(lapply(hit, as.integer), names(footprints))
}

# Conditional-on-reference point-process likelihood: each reference's sites are compared only with cells of its own
# footprint (a reference-specific intensity level cancels). beta has no intercept. Ridge prior on all slopes.
fit_cond_ref <- function(Xq, a, cells, Xs, w = NULL, prior_sd = PRIOR_SD, maxit = 80L, tol = 1e-10) {
  R <- length(cells); p <- ncol(Xq)
  if (is.null(w)) w <- rep(1, R)
  nr <- vapply(Xs, nrow, integer(1))
  sx <- Reduce(`+`, Map(function(X, wr) wr * colSums(X), Xs, w))
  pen <- rep(1 / prior_sd^2, p)
  pieces <- function(b) {
    g <- sx - pen * b; H <- -diag(pen, p); ll <- -0.5 * sum(pen * b^2)
    for (r in seq_len(R)) {
      if (w[r] == 0) next
      Xr <- Xq[cells[[r]], , drop = FALSE]; eta <- drop(Xr %*% b)
      m <- max(eta); e <- a[cells[[r]]] * exp(eta - m); s <- sum(e); P <- e / s
      mu <- colSums(Xr * P)
      ll <- ll + w[r] * (sum(drop(Xs[[r]] %*% b)) - nr[r] * (m + log(s)))
      g <- g - w[r] * nr[r] * mu
      H <- H - w[r] * nr[r] * (crossprod(Xr * P, Xr) - tcrossprod(mu))
    }
    list(ll = ll, g = g, H = H)
  }
  b <- rep(0, p); cur <- pieces(b); conv <- FALSE
  for (it in seq_len(maxit)) {
    step <- tryCatch(-solve(cur$H, cur$g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1; nw <- NULL
    for (h in 1:30) { nw <- pieces(b + s * step); if (is.finite(nw$ll) && nw$ll >= cur$ll) break; s <- s / 2 }
    if (!is.finite(nw$ll) || nw$ll < cur$ll) break
    d <- abs(nw$ll - cur$ll); b <- b + s * step; cur <- nw
    if (d < tol * (1 + abs(cur$ll)) && max(abs(s * step)) < 1e-6) { conv <- TRUE; break }
  }
  V <- tryCatch(solve(-cur$H), error = function(e) matrix(NA_real_, p, p))
  list(coef = b, se = sqrt(diag(V)), vcov = V, loglik = cur$ll, converged = conv,
       bad = !conv || any(abs(b) > 8) || !is.finite(cur$ll))
}

# reference-level bootstrap of the conditional estimator: resample references with replacement (weights = counts)
boot_cond_ref <- function(Xq, a, cells, Xs, B, progress = 50L) {
  R <- length(cells); out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    w <- as.numeric(stats::rmultinom(1, R, rep(1 / R, R)))
    f <- fit_cond_ref(Xq, a, cells, Xs, w)
    if (!f$bad) out[b] <- f$coef[length(f$coef)]
    if (b %% progress == 0L) message("    bootstrap ", b, " / ", B)
  }
  out
}

# unconditional profile with arbitrary cluster ids (one-step bootstrap)
profile_clusters <- function(label, sub, cl, ks, des, ex, Xq_list, a, B = 2000L) {
  Xs <- des$apply(sub)
  do.call(rbind, lapply(as.character(ks), function(k) {
    Xp <- cbind(int = 1, Xs, z = ex[[k]]$apply(sub))
    f <- fit_real(Xp, Xq_list[[k]], a)
    base <- data.frame(clusters = label, n_sites = nrow(sub), n_clusters = length(unique(cl)), k = as.integer(k))
    if (f$bad) return(cbind(base, data.frame(beta_iqr = NA, lo = NA, hi = NA, ratio_q25_vs_q75 = NA)))
    j <- length(f$coef); b <- boot_onestep_j(f, Xp, cl, j, B)
    ci <- ci_row(b, ex[[k]]$iqr_per_sd, ex[[k]]$sd_log, f$coef[j])
    cbind(base, ci[, c("beta_iqr", "lo", "hi", "ratio_q25_vs_q75")])
  }))
}