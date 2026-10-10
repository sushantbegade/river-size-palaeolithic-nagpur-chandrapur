# tests/testthat/test-spatial_cv.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "spatial_cv.R"))

toy_cv <- function(seed = 1) {
  set.seed(seed); n <- 4000; a <- rep(0.25, n)
  x <- stats::rnorm(n); z <- stats::rnorm(n)
  blk <- sample(paste0("b", 1:8), n, TRUE)
  cnt <- simulate_counts(0.2 * x + 0.8 * z, a, 250); idx <- rep(seq_len(n), cnt)
  Xq <- list(M0 = cbind(int = rep(1, n)), M1_nuisance = cbind(int = 1, x = x), M2 = cbind(int = 1, x = x, z = z))
  Xp <- lapply(Xq, function(m) m[idx, , drop = FALSE])
  list(a = a, blk = blk, site_blk = blk[idx], Xq = Xq, Xp = Xp, n_sites = length(idx))
}

test_that("every site is held out exactly once", {
  d <- toy_cv()
  r <- cv_run(d$Xq, d$Xp, d$a, d$blk, d$site_blk)
  one <- r[r$model == "M1_nuisance", ]
  expect_equal(sum(one$n_test_sites), d$n_sites)
  expect_equal(sum(one$n_test_cells), length(d$a))
})

test_that("the true predictor improves held-out log-likelihood; the null model does not", {
  d <- toy_cv()
  r <- cv_run(d$Xq, d$Xp, d$a, d$blk, d$site_blk)
  s <- summarise_cv(r, d$n_sites)
  expect_gt(s$delta_vs_nuisance[s$model == "M2"], 5)
  expect_gt(s$folds_improved[s$model == "M2"], 5)
  expect_lt(s$delta_vs_nuisance[s$model == "M0"], 0)
})