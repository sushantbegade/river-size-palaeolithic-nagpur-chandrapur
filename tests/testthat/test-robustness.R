# tests/testthat/test-robustness.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "robustness.R"))

test_that("thin_greedy enforces the minimum distance and keeps a maximal-ish subset", {
  set.seed(1); xy <- cbind(stats::runif(300, 0, 20000), stats::runif(300, 0, 20000))
  i <- thin_greedy(xy, 2000)
  d <- as.matrix(stats::dist(xy[i, ])); diag(d) <- Inf
  expect_gte(min(d), 2000)
  expect_lt(length(i), 300); expect_gt(length(i), 20)
})

test_that("thin_greedy keeps everything when points are far apart", {
  xy <- cbind(seq(0, 90000, by = 10000), 0)
  expect_equal(length(thin_greedy(xy, 2000)), 10L)
})