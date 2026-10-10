# tests/testthat/test-survey_bias.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "survey_bias.R"))

test_that("conditional likelihood with one all-cell footprint equals the point-process slopes", {
  set.seed(1); n <- 3000; a <- rep(0.25, n)
  X <- cbind(int = 1, x = stats::rnorm(n), z = stats::rnorm(n))
  cnt <- simulate_counts(0.3 * X[, 2] + 0.7 * X[, 3], a, 200); idx <- rep(seq_len(n), cnt)
  ppm <- fit_ppm(X[idx, ], X, a)
  cr <- fit_cond_ref(X[, -1], a, list(seq_len(n)), list(X[idx, -1]), prior_sd = 1e6)
  expect_true(cr$converged)
  expect_lt(max(abs(cr$coef - ppm$coef[-1])), 1e-4)
})

test_that("conditional estimator removes reference-specific effort levels", {
  set.seed(2); n <- 6000; a <- rep(0.25, n)
  X <- cbind(x = stats::rnorm(n), z = stats::rnorm(n))
  fp <- split(sample(n), rep(1:6, length.out = n)); fp <- lapply(fp, sort)         # six footprints
  effort <- c(1, 3, 0.3, 2, 0.5, 4)
  Xs <- lapply(1:6, function(r) { ix <- fp[[r]]; lam <- effort[r] * exp(0.8 * X[ix, 2]) * 0.25; k <- stats::rpois(length(ix), lam * 40 / sum(lam)); X[rep(ix, k), , drop = FALSE] })
  cr <- fit_cond_ref(X, a, fp, Xs, prior_sd = 1e6)
  expect_true(cr$converged)
  expect_lt(abs(cr$coef[2] - 0.8), 0.25)
})

test_that("reference bootstrap returns B draws with spread", {
  set.seed(3); n <- 3000; a <- rep(0.25, n); X <- cbind(x = stats::rnorm(n), z = stats::rnorm(n))
  fp <- lapply(split(sample(n), rep(1:8, length.out = n)), sort)
  Xs <- lapply(1:8, function(r) { ix <- fp[[r]]; lam <- exp(0.6 * X[ix, 2]); k <- stats::rpois(length(ix), lam * 25 / sum(lam)); X[rep(ix, k), , drop = FALSE] })
  b <- boot_cond_ref(X, a, fp, Xs, 30L, progress = 100L)
  expect_equal(length(b), 30L); expect_gt(stats::sd(b, na.rm = TRUE), 0.02)
})