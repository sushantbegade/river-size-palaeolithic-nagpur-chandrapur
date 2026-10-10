# tests/testthat/test-site_type.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "site_type.R"))

test_that("ridge logistic with a huge prior SD equals glm", {
  set.seed(1); n <- 400; x <- stats::rnorm(n); z <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * x - 0.5 * z))
  X <- cbind(int = 1, x = x, z = z)
  f <- fit_logit_ridge(X, y, prior_sd = 1e6)
  g <- stats::glm(y ~ x + z, family = stats::binomial())
  expect_lt(max(abs(f$coef - stats::coef(g))), 1e-5)
  expect_lt(max(abs(f$se - sqrt(diag(stats::vcov(g))))), 1e-4)
})

test_that("ridge logistic stays finite under complete separation", {
  set.seed(2); n <- 200; x <- stats::rnorm(n); y <- as.numeric(x > 0)
  f <- fit_logit_ridge(cbind(int = 1, x = x), y)
  expect_true(all(is.finite(f$coef))); expect_true(all(abs(f$coef) < 40))
})

test_that("boot_logit returns B values with sensible spread", {
  set.seed(3); n <- 300; x <- stats::rnorm(n); y <- stats::rbinom(n, 1, stats::plogis(0.9 * x))
  cl <- sample(paste0("b", 1:20), n, TRUE)
  b <- boot_logit(cbind(int = 1, x = x), y, cl, 80L)
  expect_equal(length(b), 80L); expect_gt(stats::sd(b, na.rm = TRUE), 0.05); expect_lt(abs(mean(b, na.rm = TRUE) - 0.9), 0.4)
})