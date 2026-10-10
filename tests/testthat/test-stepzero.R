# tests/testthat/test-stepzero.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))

toy <- function(n = 4000, seed = 1) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n); a <- rep(0.25, n)
  X <- cbind(int = 1, x1 = x1, x2 = x2)
  list(X = X, a = a, eta = drop(X[, -1] %*% c(0.5, -0.3)))
}

test_that("fit_ppm equals lattice Poisson glm with offset (eps -> 0)", {
  d <- toy(); cnt <- simulate_counts(d$eta, d$a, 300)
  idx <- rep(seq_along(cnt), cnt)
  f <- fit_ppm(d$X[idx, , drop = FALSE], d$X, d$a)
  g <- stats::glm(cnt ~ d$X - 1 + offset(log(d$a)), family = stats::poisson())
  expect_true(f$converged)
  expect_lt(max(abs(f$coef - stats::coef(g))), 1e-5)
  expect_lt(max(abs(f$se - sqrt(diag(stats::vcov(g))))), 1e-4)
})

test_that("fit_ppm equals down-weighted Poisson glm.fit", {
  d <- toy(); cnt <- simulate_counts(d$eta, d$a, 300)
  idx <- rep(seq_along(cnt), cnt)
  f <- fit_ppm(d$X[idx, , drop = FALSE], d$X, d$a)
  y <- c(rep(1 / SIM_EPS, length(idx)), rep(0, nrow(d$X)))
  w <- c(rep(SIM_EPS, length(idx)), d$a)
  g <- stats::glm.fit(rbind(d$X[idx, , drop = FALSE], d$X), y, weights = w, family = stats::poisson(), start = f$coef)
  expect_lt(max(abs(f$coef - g$coefficients)), 1e-6)
})

test_that("simulate_counts hits the expected total", {
  d <- toy(); set.seed(2)
  tot <- replicate(300, sum(simulate_counts(d$eta, d$a, 90)))
  expect_lt(abs(mean(tot) - 90), 3)
})

test_that("rff_field has the requested marginal SD", {
  set.seed(4)
  xy <- cbind(stats::runif(500, 0, 5e4), stats::runif(500, 0, 5e4))
  v <- c(vapply(1:200, function(i) rff_field(0.8, 1e4)(xy)[1:20], numeric(20)))
  expect_lt(abs(stats::sd(v) - 0.8), 0.08)
})