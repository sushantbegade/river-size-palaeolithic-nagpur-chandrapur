# tests/testthat/test-stepzero_b.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))

test_that("one-step cluster bootstrap approximates the full-refit bootstrap SE", {
  set.seed(11); n <- 6000
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n); a <- rep(0.25, n)
  Xq <- cbind(int = 1, x1 = x1, z = x2)
  cl_cell <- sample(paste0("b", 1:40), n, TRUE)
  cnt <- simulate_counts(0.5 * x1 + 0.4 * x2, a, 250)
  idx <- rep(seq_along(cnt), cnt)
  Xp <- Xq[idx, , drop = FALSE]
  f <- fit_ppm(Xp, Xq, a)
  b1 <- boot_onestep(f, Xp, cl_cell[idx], 2000L)
  b2 <- boot_full(Xq, a, Xp, cl_cell[idx], 150L)
  r <- stats::sd(b1) / stats::sd(b2)
  expect_gt(r, 0.8); expect_lt(r, 1.25)
  expect_lt(abs(mean(b1) - f$coef[3]), 0.05)
})

test_that("nuisance_design_c drops alluvial and keeps 7 columns", {
  set.seed(3); n <- 400
  dom <- data.frame(elev_m = stats::rnorm(n, 200, 30), slope_deg = stats::rexp(n),
                    gm4 = sample(c("pediplain", "plateau_hills", "alluvial"), n, TRUE),
                    frac_crop = stats::runif(n), frac_built = stats::runif(n)^3,
                    d_road_major = stats::rexp(n, 1 / 1000), d_settle = stats::rexp(n, 1 / 1300))
  X <- nuisance_design_c(dom)
  expect_equal(ncol(X), 7L)
  expect_false("gm_alluvial" %in% colnames(X))
})

test_that("mde_table interpolates on log ratio", {
  ps <- data.frame(sigma = 0, n_target = 93, k = 6, ratio = c(1.25, 1.5, 2, 3), rej_ci = c(0.2, 0.5, 0.8, 0.95))
  expect_equal(mde_table(ps)$mde_ratio, 2)
})