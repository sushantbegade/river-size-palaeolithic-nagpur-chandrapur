# tests/testthat/test-model_primary.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))

mk_dom <- function(n, seed) {
  set.seed(seed)
  data.frame(elev_m = stats::rnorm(n, 200, 30), slope_deg = stats::rexp(n), gm4 = sample(c("pediplain", "plateau_hills", "alluvial"), n, TRUE),
             frac_crop = stats::runif(n), frac_built = ifelse(stats::runif(n) < .7, -1e-17, stats::runif(n)^2),
             d_road_major = stats::rexp(n, 1 / 1000), d_settle = stats::rexp(n, 1 / 1300), d_k6 = stats::rexp(n, 1 / 9000))
}

test_that("apply() on the domain reproduces the domain design exactly", {
  dom <- mk_dom(800, 1)
  des <- design_from_domain(dom)
  expect_equal(unname(des$apply(dom)), unname(des$dom), tolerance = 1e-10)
  expect_false("gm_alluvial" %in% colnames(des$dom))
})

test_that("sites are standardised with domain, not site, statistics", {
  dom <- mk_dom(800, 1); st <- dom[1:20, ]
  des <- design_from_domain(dom)
  expect_gt(max(abs(colMeans(des$apply(st)))), 1e-6)
})

test_that("exposure apply() matches domain z", {
  dom <- mk_dom(500, 2); e <- exposure_domain(dom, 6)
  expect_equal(e$apply(dom), e$z_dom, tolerance = 1e-12)
})

# append to tests/testthat/test-model_primary.R
test_that("fit_real does not flag a low intercept as separation", {
  set.seed(12); n <- 3000; a <- rep(25, n)                 # 75,000 km2: log intensity about -8.2
  X <- cbind(int = 1, x = stats::rnorm(n), z = stats::rnorm(n))
  cnt <- simulate_counts(0.3 * X[, 2] + 0.3 * X[, 3], a, 20); idx <- rep(seq_len(n), cnt)
  f <- fit_real(X[idx, ], X, a)
  expect_lt(f$coef[1], -8)
  expect_false(f$bad)
})