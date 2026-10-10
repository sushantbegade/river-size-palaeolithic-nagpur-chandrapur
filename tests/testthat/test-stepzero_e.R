# tests/testthat/test-stepzero_e.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "stepzero_e.R"))

test_that("K_CLASS bins low / mid / high", {
  expect_equal(K_CLASS(c(3, 5, 6, 7, 8, 9)), c("low", "low", "mid", "mid", "high", "high"))
})

test_that("profile picks the generating predictor when the effect is strong", {
  set.seed(5); n <- 5000; a <- rep(0.25, n)
  z1 <- stats::rnorm(n); z2 <- 0.4 * z1 + sqrt(1 - 0.16) * stats::rnorm(n)
  X <- matrix(stats::rnorm(n), n, 1, dimnames = list(NULL, "x"))
  Xk <- list(`1` = cbind(int = 1, X, z = z1), `2` = cbind(int = 1, X, z = z2))
  ik <- list(`1` = 1, `2` = 1)
  hit <- vapply(1:20, function(i) {
    cnt <- simulate_counts(1.0 * z1, a, 300)
    pr <- fit_profile(cnt, Xk, a, ik)
    pr$k[which.min(pr$aic)] == 1L
  }, logical(1))
  expect_gt(mean(hit), 0.9)
})