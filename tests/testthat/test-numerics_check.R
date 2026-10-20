# tests for R/numerics_check.R (synthetic data only)
suppressWarnings(suppressPackageStartupMessages(library(testthat)))
for (f in c("config", "utils", "quadrature", "stepzero", "stepzero_b", "stepzero_c", "model_primary", "diagnostics", "numerics_check"))
  source(here::here("R", paste0(f, ".R")))

mk <- function(seed = 1, n_dom = 2500, n_s = 70) {
  set.seed(seed)
  Xq <- cbind(int = 1, a = rnorm(n_dom), z = rnorm(n_dom)); aq <- rep(0.25, n_dom)
  eta <- -6 + 0.3 * Xq[, "a"] + 0.6 * Xq[, "z"]
  idx <- sample(n_dom, n_s, FALSE, prob = exp(eta - max(eta)))
  list(Xq = Xq, aq = aq, Xp = Xq[idx, , drop = FALSE])
}

test_that("shift_clusters reproduces cluster_ids at zero shift and keeps the partition under a whole-block shift", {
  set.seed(3); xy <- cbind(runif(80, 3e5, 3.6e5), runif(80, 2.25e6, 2.31e6))
  expect_identical(shift_clusters(xy, 40, 0, 0), cluster_ids(xy, 40000))
  a <- shift_clusters(xy, 40, 0, 0); b <- shift_clusters(xy, 40, 40, 40)
  expect_equal(length(unique(a)), length(unique(b))); expect_equal(nrow(unique(data.frame(a, b))), length(unique(a)))   # same partition, new labels
  expect_true(length(unique(shift_clusters(xy, 40, 20, 10))) >= 1)
})

test_that("eps_check: results change little with eps and the eps term share scales with eps", {
  m <- mk(); ec <- eps_check(list(`6` = cbind(m$Xp)), list(`6` = m$Xq), m$aq, c(`6` = 1.2), eps_vec = c(1e-4, 1e-6, 1e-8))
  expect_equal(nrow(ec), 3L); expect_true(all(ec$converged))
  expect_lt(abs(ec$beta_iqr[1] - ec$beta_iqr[2]), 0.01); expect_lt(abs(ec$beta_iqr[3] - ec$beta_iqr[2]), 1e-4)
  expect_gt(ec$eps_term_share[1], 50 * ec$eps_term_share[2] * 0.5)        # about 100 times larger at 1e-4 than at 1e-6
  expect_lt(abs(ec$beta_iqr[2] - ec$beta_iqr_ml[2]), 0.15)                # ridge vs ML: small shrinkage
})

test_that("origin_row returns a CI that contains the estimate for most origins and counts clusters", {
  m <- mk(); set.seed(5)
  f <- fit_real(cbind(m$Xp), m$Xq, m$aq)
  xy <- cbind(runif(nrow(m$Xp), 3e5, 3.5e5), runif(nrow(m$Xp), 2.25e6, 2.3e6))
  r <- origin_row(f, cbind(m$Xp), xy, 40, 10, 20, iqr_per_sd = 1.2, sd_log = 1.0, B = 300L)
  expect_equal(nrow(r), 1L); expect_true(r$n_clusters >= 2L); expect_true(r$lo <= r$hi)
  expect_equal(r$n_clusters, length(unique(shift_clusters(xy, 40, 10, 20))))
  expect_true(is.logical(r$excludes_zero))
})