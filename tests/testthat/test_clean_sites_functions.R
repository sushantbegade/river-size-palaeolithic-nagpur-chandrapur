# Unit tests that need NO real data (synthetic inputs only)
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
source(here::here("R", "config.R")); source(here::here("R", "utils.R"))
source(here::here("R", "clean_sites.R"))

fake_raw <- function(chron = c("Lower Palaeolithic", "Middle and Upper Palaeolithic")) {
  n <- length(chron)
  tibble::tibble(
    site_id = 1:n, site_name = paste("S", 1:n), lat = rep(21, n), lon = rep(79, n),
    district = rep("Nagpur", n), state = "Maharashtra", river = c("-", "Kanhan")[seq_len(n) %% 2 + 1],
    basin = "Wainganga", precision = "GIS Map Based", assemblage = "x",
    site_type_raw = c("Open-Air", "Riverbed")[seq_len(n) %% 2 + 1], raw_material = "Not Reported",
    context = "Not Reported", culture = "x", chron = chron, pub_year = "2000",
    ref = "Ref, 2000", source_type = "Survey")
}
cw <- tibble::tibble(chron = c("Lower Palaeolithic", "Middle and Upper Palaeolithic"),
                     LP = c(1, 0), MP = c(0, 1), UP = c(0, 1))

test_that("period crosswalk works and flags multi-period", {
  d <- clean_sites(fake_raw(), cw)
  expect_equal(d$period_class, c("LP", "multi"))
  expect_equal(d$n_periods, c(1L, 2L))
})

test_that("unknown chronology term stops with a clear error", {
  expect_error(clean_sites(fake_raw(c("Lower Palaeolithic", "Neolithic")), cw),
               "missing from docs/crosswalk_period.csv")
})

test_that("site type spellings collapse", {
  r <- fake_raw(rep("Lower Palaeolithic", 3))
  r$site_type_raw <- c("Open-air", "Open-Air", "Riverbed")
  d <- clean_sites(r, cw)
  expect_equal(d$site_type, c("open_air", "open_air", "riverbed"))
  expect_equal(sum(d$S1_openair), 2L)
})

test_that("na_codes converts placeholders to NA", {
  expect_equal(na_codes(c("-", "Not Reported", " x ", "")), c(NA, NA, "x", NA))
})

test_that("greedy thinning respects the distance", {
  xy <- cbind(c(0, 500, 3000), c(0, 0, 0))
  expect_equal(greedy_thin_n(xy, 1000, 1), 2L)
  expect_equal(greedy_thin_n(xy[0, , drop = FALSE], 1000, 1), 0L)
})

test_that("coarsening snaps to 1 km cell centres and drops exact coordinates", {
  d <- clean_sites(fake_raw(), cw); d <- add_utm(d)
  d$x[1] <- 1234; d$y[1] <- 5678
  # recompute geographic centre from UTM for the first row only
  out <- coarsen_sites(d, 1000)
  expect_equal(out$x_c[1], 1500); expect_equal(out$y_c[1], 5500)
  expect_false(any(c("lat", "lon", "x", "y", "site_name") %in% names(out)))
})
