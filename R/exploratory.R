# R/exploratory.R -- functions for script 17 (EXPLORATORY items declared in the frozen design: distance shape, single-period models, period contrast)
# Needs: model_primary.R, diagnostics.R (boot_onestep_j, subset_profile), stepzero_b.R (cluster_ids)

# five equal-area distance classes from the DOMAIN distribution of the standardised log distance; class 3 is the reference
bin_design <- function(z_dom, z_site, d_dom_m, nb = 5L, ref = 3L) {
  br <- stats::quantile(z_dom, seq(0, 1, length.out = nb + 1L), names = FALSE); br[1] <- -Inf; br[nb + 1L] <- Inf
  bd <- cut(z_dom, br, labels = FALSE); bs <- cut(z_site, br, labels = FALSE)
  keep <- setdiff(seq_len(nb), ref)
  dum <- function(b) sapply(keep, function(j) as.numeric(b == j))
  list(Dq = dum(bd), Dp = dum(bs), keep = keep, bin_dom = bd, bin_site = bs,
       km_range = t(vapply(seq_len(nb), function(j) range(d_dom_m[bd == j]) / 1000, numeric(2))))
}

binned_fit <- function(label, k, sub, ex, des, dom, a, block_km, B = 2000L) {
  bd <- bin_design(ex$z_dom, ex$apply(sub), dom[[paste0("d_k", k)]])
  Xs <- des$apply(sub)
  Xq <- cbind(int = 1, des$dom, bd$Dq); Xp <- cbind(int = 1, Xs, bd$Dp)
  f <- fit_real(Xp, Xq, a)
  nb <- nrow(bd$km_range); p <- ncol(Xq); idx <- p - length(bd$keep) + seq_along(bd$keep)
  cl <- cluster_ids(as.matrix(sub[, c("x", "y")]), block_km * 1000)
  area <- vapply(seq_len(nb), function(j) sum(a[bd$bin_dom == j]), numeric(1))
  nsite <- vapply(seq_len(nb), function(j) sum(bd$bin_site == j), numeric(1))
  out <- data.frame(subset = label, k = k, bin = seq_len(nb), km_from = round(bd$km_range[, 1], 2), km_to = round(bd$km_range[, 2], 2),
                    area_km2 = round(area), n_sites = nsite, sites_per_100km2 = round(100 * nsite / area, 2),
                    rel_intensity = NA_real_, lo = NA_real_, hi = NA_real_, fit_ok = !f$bad)
  out$rel_intensity[out$bin == 3] <- 1
  if (!f$bad) for (q in seq_along(bd$keep)) {
    j <- idx[q]; b <- boot_onestep_j(f, Xp, cl, j, B); r <- out$bin == bd$keep[q]
    out$rel_intensity[r] <- round(exp(f$coef[j]), 2)
    if (sum(is.finite(b)) >= 50L) { qq <- stats::quantile(b, c(.025, .975), na.rm = TRUE, names = FALSE); out$lo[r] <- round(exp(qq[1]), 2); out$hi[r] <- round(exp(qq[2]), 2) }
  }
  out
}

# period contrast: resample 40 km blocks across all single-period sites, refit each period, record the differences in beta per IQR
boot_period_contrast <- function(subs, des, ex, Xq, a, block_km, B = 300L, progress = 50L) {
  allsub <- do.call(rbind, subs)
  blk_all <- cluster_ids(as.matrix(allsub[, c("x", "y")]), block_km * 1000); ub <- unique(blk_all)
  rows <- lapply(subs, function(s) split(seq_len(nrow(s)), factor(cluster_ids(as.matrix(s[, c("x", "y")]), block_km * 1000), levels = ub)))
  Xs <- lapply(subs, function(s) cbind(int = 1, des$apply(s), z = ex$apply(s)))
  nm <- names(subs); out <- matrix(NA_real_, B, length(nm), dimnames = list(NULL, nm))
  for (b in seq_len(B)) {
    pick <- sample(ub, length(ub), TRUE)
    for (p in nm) {
      i <- unlist(rows[[p]][pick], use.names = FALSE)
      if (length(i) < 8L) next
      f <- fit_real(Xs[[p]][i, , drop = FALSE], Xq, a)
      if (!f$bad) out[b, p] <- f$coef[length(f$coef)] * ex$iqr_per_sd
    }
    if (b %% progress == 0L) message("    bootstrap ", b, " / ", B)
  }
  out
}