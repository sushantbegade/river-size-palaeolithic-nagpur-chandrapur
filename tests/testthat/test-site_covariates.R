# tests/testthat/test-site_covariates.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "detectability.R"))
source(here::here("R", "site_covariates.R"))

test_that("cluster_structure counts occupied blocks and singletons", {
  xy <- cbind(x = c(1000, 2000, 3000, 25000, 90000), y = c(1000, 1500, 500, 1000, 1000))
  cs <- cluster_structure(xy, c(10, 20), list(`10` = 20L, `20` = 9L))
  expect_equal(cs$n_clusters_occupied, c(3L, 3L))
  expect_equal(cs$singleton_clusters, c(2L, 2L))
  expect_equal(cs$max_cluster_share, c(0.6, 0.6))
})

test_that("qa_against_grid flags a real difference but not float32 noise", {
  nm <- c("elev_m", "slope_deg", paste0("d_k", SIZE_K), "d_lineament", "d_road_major", "d_road_any", "d_settle", "d_town",
          paste0("frac_", names(WC_GROUPS)))
  g <- as.data.frame(matrix(1, 10, length(nm), dimnames = list(NULL, nm))); g$gm4 <- "pediplain"
  n <- g
  n$elev_m <- n$elev_m + 7e-6; n$frac_crop <- n$frac_crop + 3e-8
  expect_equal(sum(qa_against_grid(n, g)$n_mismatch), 0)
  n$d_k6[3] <- n$d_k6[3] + 0.5
  expect_equal(qa_against_grid(n, g)$n_mismatch[qa_against_grid(n, g)$variable == "d_k6"], 1)
})