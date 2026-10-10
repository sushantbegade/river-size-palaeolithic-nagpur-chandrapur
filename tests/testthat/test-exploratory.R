# tests/testthat/test-exploratory.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "exploratory.R"))

test_that("bin_design gives equal-area classes and drops the reference bin", {
  set.seed(1); z <- stats::rnorm(5000); d <- exp(z) * 1000
  bd <- bin_design(z, z[1:200], d)
  expect_equal(as.numeric(table(bd$bin_dom)), rep(1000, 5))
  expect_equal(ncol(bd$Dq), 4L); expect_false(3L %in% bd$keep)
  expect_true(all(diff(bd$km_range[, 1]) > 0))
})

test_that("period contrast bootstrap returns finite draws for each period", {
  set.seed(2); n <- 3000; a <- rep(0.25, n)
  dom <- data.frame(elev_m = stats::rnorm(n, 200, 30), slope_deg = stats::rexp(n), gm4 = sample(c("pediplain", "plateau_hills", "alluvial"), n, TRUE),
                    frac_crop = stats::runif(n), frac_built = stats::runif(n)^3, d_road_major = stats::rexp(n, 1/1000), d_settle = stats::rexp(n, 1/1300),
                    d_k6 = stats::rexp(n, 1/9000), x = stats::runif(n, 0, 1e5), y = stats::runif(n, 0, 1e5))
  des <- design_from_domain(dom); e6 <- exposure_domain(dom, 6); Xq <- cbind(int = 1, des$dom, z = e6$z_dom)
  mk <- function(m) { i <- sample(n, m, prob = exp(0.5 * e6$z_dom)); dom[i, ] }
  subs <- list(LP = mk(40), MP = mk(50), UP = mk(30))
  bm <- boot_period_contrast(subs, des, e6, Xq, a, 20, B = 20L, progress = 100L)
  expect_equal(dim(bm), c(20L, 3L)); expect_gt(sum(is.finite(bm)), 40)
})