# tests/testthat/test-diagnostics.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))

test_that("boot_onestep_j equals boot_onestep for the last coefficient (same seed)", {
  set.seed(8); n <- 3000; a <- rep(0.25, n)
  X <- cbind(int = 1, x = stats::rnorm(n), z = stats::rnorm(n))
  cnt <- simulate_counts(0.4 * X[, 2] + 0.5 * X[, 3], a, 150); idx <- rep(seq_len(n), cnt)
  cl <- sample(paste0("b", 1:25), n, TRUE)
  f <- fit_ppm_ridge(X[idx, ], X, a)
  set.seed(1); b1 <- boot_onestep(f, X[idx, ], cl[idx], 300L)
  set.seed(1); b2 <- boot_onestep_j(f, X[idx, ], cl[idx], 3L, 300L)
  expect_equal(b1, b2)
})

test_that("leave_one_out drops exactly one group and reports change vs full", {
  set.seed(9); n <- 3000; a <- rep(0.25, n)
  dom <- data.frame(elev_m = stats::rnorm(n, 200, 30), slope_deg = stats::rexp(n), gm4 = sample(c("pediplain", "plateau_hills", "alluvial"), n, TRUE),
                    frac_crop = stats::runif(n), frac_built = stats::runif(n)^3, d_road_major = stats::rexp(n, 1/1000),
                    d_settle = stats::rexp(n, 1/1300), d_k6 = stats::rexp(n, 1/9000))
  des <- design_from_domain(dom); ex6 <- exposure_domain(dom, 6)
  Xq <- cbind(int = 1, des$dom, z = ex6$z_dom)
  idx <- sample(n, 60, prob = exp(-0.4 * ex6$z_dom))
  sub <- dom[idx, ]; sub$grp <- rep(paste0("g", 1:6), length.out = 60)
  r <- leave_one_out(sub, "grp", des, ex6, Xq, a)
  expect_equal(nrow(r), 6L); expect_true(all(r$n_remaining == 50L))
  expect_true(all(diff(abs(r$change_vs_full)) <= 1e-12))
})