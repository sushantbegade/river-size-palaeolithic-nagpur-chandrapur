# R/site_type.R -- functions for script 13 (EXPLORATORY: does the site-type label depend on distance to large rivers?)
# Needs: config.R, quadrature.R, stepzero.R, stepzero_b.R, stepzero_c.R (PRIOR_SD), model_primary.R, diagnostics.R

# ridge-penalised logistic regression (Newton / IRLS with step halving); column 1 = intercept (unpenalised)
fit_logit_ridge <- function(X, y, prior_sd = PRIOR_SD, maxit = 60L, tol = 1e-10) {
  p <- ncol(X); pen <- c(0, rep(1 / prior_sd^2, p - 1L))
  b <- rep(0, p)
  ll <- function(b) {
    eta <- drop(X %*% b)
    v <- sum(y * eta - log1p(exp(eta))) - 0.5 * sum(pen * b^2)
    if (is.finite(v)) v else -Inf
  }
  cur <- ll(b); conv <- FALSE
  for (it in seq_len(maxit)) {
    pr <- stats::plogis(drop(X %*% b)); w <- pr * (1 - pr)
    g <- crossprod(X, y - pr) - pen * b
    H <- -(crossprod(X * w, X) + diag(pen, p))
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1; new <- -Inf
    for (h in 1:30) { new <- ll(b + s * drop(step)); if (new >= cur) break; s <- s / 2 }
    if (!is.finite(new) || new < cur) break
    b <- b + s * drop(step)
    d <- abs(new - cur); cur <- new
    if (d < tol * (1 + abs(cur)) && max(abs(s * step)) < 1e-6) { conv <- TRUE; break }
  }
  pr <- stats::plogis(drop(X %*% b)); w <- pr * (1 - pr)
  V <- tryCatch(solve(crossprod(X * w, X) + diag(pen, p)), error = function(e) matrix(NA_real_, p, p))
  ll_plain <- sum(y * drop(X %*% b) - log1p(exp(drop(X %*% b))))
  list(coef = drop(b), se = sqrt(diag(V)), vcov = V, converged = conv, loglik = ll_plain, aic = -2 * ll_plain + 2 * p,
       bad = !conv || any(abs(b[-1]) > 8) || !is.finite(ll_plain))
}

# cluster bootstrap (full refits) of the last coefficient
boot_logit <- function(X, y, cl, B = 1000L) {
  ucl <- unique(cl); nc <- length(ucl)
  rows <- split(seq_along(y), factor(cl, levels = ucl))
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- unlist(rows[sample.int(nc, nc, TRUE)], use.names = FALSE)
    if (length(unique(y[idx])) < 2L) next
    f <- fit_logit_ridge(X[idx, , drop = FALSE], y[idx])
    if (!f$bad) out[b] <- f$coef[ncol(X)]
  }
  out
}