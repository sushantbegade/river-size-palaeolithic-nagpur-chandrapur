# R/source_struct.R -- functions for script 21 (EXPLORATORY, post-freeze): source structure of the S1 result and of the S1/S3 labels.
# Uses the frozen helpers unchanged. Needs: config.R, quadrature.R (log_d), stepzero.R, stepzero_b.R (cluster_ids), stepzero_c.R,
# model_primary.R (fit_real, ci_row), domain_sens.R (ref_family).

# counts per first-named-author family (no coordinates): strings, S1, S3 and other sites inside the analysis domain
family_table <- function(sc_dom, fam) {
  d <- data.frame(family = fam, ref = sc_dom$ref_clean, s1 = sc_dom$S1_openair %in% TRUE, s3 = sc_dom$S3_riverbed %in% TRUE, stringsAsFactors = FALSE)
  t <- do.call(rbind, lapply(split(d, d$family), function(x)
    data.frame(family = x$family[1], n_strings = length(unique(x$ref)), n_sites = nrow(x), n_S1 = sum(x$s1), n_S3 = sum(x$s3), n_other = sum(!x$s1 & !x$s3))))
  t <- t[order(-t$n_S1, -t$n_S3), ]; rownames(t) <- NULL
  t$share_of_S1 <- round(t$n_S1 / sum(t$n_S1), 3)
  t
}

# leave-one-group-out refits at one k, with full-precision estimates (the frozen leave_one_out() rounds to 3 decimals); loglik is kept to detect true ties
loo_full <- function(sub, group, des, ex_k, Xq, a) {
  Xp <- cbind(int = 1, des$apply(sub), z = ex_k$apply(sub)); g <- sub[[group]]
  full <- fit_real(Xp, Xq, a); j <- length(full$coef)
  out <- do.call(rbind, lapply(sort(unique(g)), function(u) {
    keep <- g != u
    f <- fit_real(Xp[keep, , drop = FALSE], Xq, a)
    data.frame(left_out = u, n_left_out = sum(!keep), n_remaining = sum(keep),
               beta_iqr = if (f$bad) NA_real_ else f$coef[j] * ex_k$iqr_per_sd,
               loglik = if (f$bad) NA_real_ else f$loglik, fit_ok = !f$bad, stringsAsFactors = FALSE)
  }))
  out$change_vs_full <- out$beta_iqr - full$coef[j] * ex_k$iqr_per_sd
  out <- out[order(-abs(out$change_vs_full)), ]; rownames(out) <- NULL
  out$tied_beta_6dp <- duplicated(round(out$beta_iqr, 6)) | duplicated(round(out$beta_iqr, 6), fromLast = TRUE)
  out
}

# logistic regression of y on z (+ X). stratum constant: ordinary logistic regression with a common intercept.
# Several strata: CONDITIONAL logistic regression (survival::clogit), which removes the stratum effects exactly. An unconditional fixed-effect
# logit (conditional = FALSE) is kept only for comparison: with small strata it is biased away from 0. Returns the coefficient of z.
fe_logit <- function(y, z, stratum, X = NULL, conditional = TRUE) {
  n <- length(y); st <- as.integer(factor(stratum))
  if (max(st) == 1L || !conditional) {
    mm <- if (max(st) == 1L) matrix(1, n, 1) else 1 * outer(st, seq_len(max(st)), "==")
    colnames(mm) <- paste0("st", seq_len(ncol(mm)))
    Xall <- cbind(mm, if (!is.null(X)) X, z = z)
    f <- suppressWarnings(stats::glm.fit(Xall, y, family = stats::binomial(), control = stats::glm.control(maxit = 100)))
    b <- f$coefficients[ncol(Xall)]
    return(list(coef = unname(b), ok = isTRUE(f$converged) && all(is.finite(f$coefficients)) && abs(b) < 20))
  }
  strata <- survival::strata                                  # found by the formula below (environment of this function)
  d <- data.frame(y = y, z = z, st = st, one = 1)
  rhs <- "z"
  if (!is.null(X)) { colnames(X) <- paste0("v", seq_len(ncol(X))); d <- cbind(d, as.data.frame(X)); rhs <- paste(c("z", colnames(X)), collapse = " + ") }
  fo <- stats::as.formula(paste("survival::Surv(one, y) ~", rhs, "+ strata(st)"))
  # exact-ties Cox partial likelihood with one event time = conditional logistic likelihood (what survival::clogit fits internally)
  f <- tryCatch(suppressWarnings(survival::coxph(fo, data = d, ties = "exact")), error = function(e) NULL)
  if (is.null(f)) return(list(coef = NA_real_, ok = FALSE))
  b <- unname(stats::coef(f)["z"])
  list(coef = b, ok = is.finite(b) && abs(b) < 20 && all(is.finite(stats::coef(f))))
}

# cluster bootstrap over `cluster` (whole clusters resampled with replacement). fe = TRUE: every resampled cluster gets its own fixed effect
# (duplicates are separate strata); fe = FALSE: one common intercept (pooled model)
boot_logit <- function(y, z, cluster, X = NULL, B = 1000L, fe = TRUE) {
  uc <- unique(cluster); nc <- length(uc)
  by <- split(seq_along(y), factor(cluster, levels = uc))
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    pick <- sample.int(nc, nc, TRUE)
    idx <- unlist(by[pick], use.names = FALSE)
    st <- if (fe) rep(seq_along(pick), lengths(by[pick])) else rep(1L, length(idx))
    f <- fe_logit(y[idx], z[idx], st, if (!is.null(X)) X[idx, , drop = FALSE])
    if (f$ok) out[b] <- f$coef
  }
  out
}

# one table row: log-odds of riverbed per IQR of log distance (positive = riverbed label more likely farther from the river class)
label_row <- function(model, y, z, stratum, cluster, iqr_per_sd, X = NULL, fe = TRUE, B = 1000L) {
  f <- fe_logit(y, z, if (fe) stratum else rep(1L, length(y)), X)
  b <- boot_logit(y, z, cluster, X, B, fe)
  q <- stats::quantile(b[is.finite(b)], c(0.025, 0.975), names = FALSE)
  data.frame(model = model, n_sites = length(y), n_S1 = sum(y == 0), n_S3 = sum(y == 1), n_strata = length(unique(stratum)),
             log_odds_per_IQR = round(f$coef * iqr_per_sd, 3), lo = round(q[1] * iqr_per_sd, 3), hi = round(q[2] * iqr_per_sd, 3),
             odds_ratio_near_vs_far = round(exp(-f$coef * iqr_per_sd), 2), fit_ok = f$ok, n_boot_ok = sum(is.finite(b)), stringsAsFactors = FALSE)
}