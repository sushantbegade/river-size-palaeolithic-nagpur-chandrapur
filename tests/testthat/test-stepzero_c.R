# tests/testthat/test-stepzero_c.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))

test_that("ridge fit with a huge prior SD equals the ML fit", {
  set.seed(1); n <- 3000; a <- rep(0.25, n)
  X <- cbind(int = 1, x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  cnt <- simulate_counts(0.5 * X[, 2] - 0.3 * X[, 3], a, 250); idx <- rep(seq_len(n), cnt)
  m <- fit_ppm(X[idx, ], X, a); r <- fit_ppm_ridge(X[idx, ], X, a, prior_sd = 1e6)
  expect_lt(max(abs(m$coef - r$coef)), 1e-6)
})

test_that("ridge fit stays finite when a binary covariate has no events", {
  set.seed(2); n <- 3000; a <- rep(0.25, n)
  d <- as.numeric(seq_len(n) > 2700)
  X <- cbind(int = 1, d = d, x = stats::rnorm(n))
  cnt <- stats::rpois(n, a * ifelse(d == 1, 0, 0.02) * exp(0.3 * X[, 3])); idx <- rep(seq_len(n), cnt)
  r <- fit_ppm_ridge(X[idx, ], X, a)
  expect_true(r$converged)
  expect_true(all(abs(r$coef) < 8))
  expect_true(all(is.finite(r$se)))
})

test_that("one-step bootstrap works with the penalised fit", {
  set.seed(3); n <- 4000; a <- rep(0.25, n)
  X <- cbind(int = 1, x1 = stats::rnorm(n), z = stats::rnorm(n))
  cnt <- simulate_counts(0.4 * X[, 2] + 0.3 * X[, 3], a, 200); idx <- rep(seq_len(n), cnt)
  cl <- sample(paste0("b", 1:30), n, TRUE)
  r <- fit_ppm_ridge(X[idx, ], X, a)
  b <- boot_onestep(r, X[idx, ], cl[idx], 500L)
  expect_gt(stats::sd(b), 0.02); expect_lt(stats::sd(b), 0.5)
  expect_lt(abs(mean(b) - r$coef[3]), 0.05)
})