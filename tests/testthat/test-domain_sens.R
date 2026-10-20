# tests for R/domain_sens.R (synthetic data only; no real site or layer is read)
suppressPackageStartupMessages({ library(testthat); library(sf) })
for (f in c("config", "utils", "quadrature", "stepzero", "stepzero_b", "stepzero_c", "model_primary", "diagnostics", "domain_sens"))
  source(here::here("R", paste0(f, ".R")))

mk_dom <- function(n = 600, alluvial = TRUE, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 3e5, 3.4e5), y = runif(n, 2.25e6, 2.29e6))
  d$elev_m <- 300 + rnorm(n, 0, 20); d$slope_deg <- runif(n, 0.2, 6)
  d$gm4 <- sample(if (alluvial) c("pediplain", "plateau_hills", "alluvial") else c("pediplain", "plateau_hills"), n, TRUE, prob = if (alluvial) c(.6, .3, .1) else c(.7, .3))
  d$frac_crop <- runif(n); d$frac_built <- runif(n, 0, .1)
  d$d_road_major <- runif(n, 50, 9000); d$d_settle <- runif(n, 50, 6000)
  d$d_k6 <- runif(n, 20, 30000); d$d_k2 <- runif(n, 10, 9000); d$area_km2 <- 0.25
  d
}

test_that("ref_family takes the first-named author and keeps odd strings", {
  expect_equal(ref_family(c("Srinivasan, 1961; Joshi, 2018", "Badam et al., 2000", "IAR, 1960, pg. 72", "  Nikoshey, 2016")),
               c("Srinivasan", "Badam", "IAR", "Nikoshey"))
  expect_equal(ref_family("1960 survey"), "1960 survey")
})

test_that("design_for_domain equals design_from_domain when alluvial cells exist, and works when they do not", {
  d1 <- mk_dom(alluvial = TRUE); a <- design_for_domain(d1); b <- design_from_domain(d1)
  expect_equal(unname(as.matrix(a$dom)), unname(as.matrix(b$dom)), tolerance = 1e-12)
  expect_equal(colnames(a$dom), colnames(b$dom))
  d0 <- mk_dom(alluvial = FALSE)
  expect_error(design_from_domain(d0), "Constant or NA")
  g <- design_for_domain(d0)
  expect_equal(colnames(g$dom), colnames(b$dom))
  expect_equal(unname(g$apply(d0)), unname(as.matrix(g$dom)), tolerance = 1e-12)   # apply() on the domain reproduces the stored design
  expect_equal(unname(colMeans(g$dom)), rep(0, ncol(g$dom)), tolerance = 1e-10)
})

test_that("baseline-IQR rescaling equals ci_row with the baseline IQR (linearity)", {
  set.seed(1); b <- rnorm(400, 0.9, 0.2); est <- 0.85
  own <- ci_row(b, iqr_per_sd = 1.2, sd_log = 1.1, est_std = est)
  bas <- ci_row(b, iqr_per_sd = 1.5, sd_log = 1.1, est_std = est)
  prof <- data.frame(k = 6L, iqr_log = 1.2 * 1.1, beta_iqr = own$beta_iqr, lo = own$lo, hi = own$hi)
  r <- to_base_scale(prof, c(`6` = 1.5 * 1.1))
  expect_equal(r$beta_base, bas$beta_iqr, tolerance = 0.002); expect_equal(r$lo_base, bas$lo, tolerance = 0.002); expect_equal(r$hi_base, bas$hi, tolerance = 0.002)
})

test_that("a domain built from a subset of sites is inside the domain built from all sites", {
  set.seed(5)
  study <- sf::st_sfc(sf::st_polygon(list(rbind(c(3e5, 2.25e6), c(3.6e5, 2.25e6), c(3.6e5, 2.31e6), c(3e5, 2.31e6), c(3e5, 2.25e6)))), crs = CRS_ANALYSIS)
  n <- 80; ref <- sample(1:12, n, TRUE)
  cen <- cbind(runif(12, 3.05e5, 3.55e5), runif(12, 2.255e6, 2.305e6))
  xy <- cen[ref, ] + matrix(rnorm(2 * n, 0, 3000), ncol = 2)
  s <- sf::st_as_sf(data.frame(site_id = 1:n, ref_id = ref, x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
  full <- survey_domain(s, study)
  for (keep in list(s$ref_id %in% 1:6, seq_len(n) %% 3 == 0, s$ref_id == 1)) {
    sub <- subset_domain_poly(s, keep, study)
    out <- sum(as.numeric(sf::st_area(suppressWarnings(sf::st_difference(sub, full)))))
    expect_lt(out / as.numeric(sf::st_area(full)), 1e-4)   # only arc-approximation slivers lie outside
  }
  expect_error(subset_domain_poly(s, rep(FALSE, n), study))
})

test_that("write_s20 refuses coordinate-like columns", {
  tf <- tempfile(fileext = ".csv")
  expect_error(write_s20(data.frame(x = 1, v = 2), tf), "coordinate-like")
  expect_error(write_s20(data.frame(Latitude = 1), tf), "coordinate-like")
  expect_silent(suppressMessages(write_s20(data.frame(domain = "D0", beta = 1), tf)))
})

test_that("profile_domain and fullrefit_k return one row per k with consistent scales", {
  set.seed(9)
  dm <- mk_dom(1500, seed = 9); des <- design_for_domain(dm)
  dm$d_k3 <- runif(nrow(dm), 20, 20000)
  s <- dm[sample(nrow(dm), 60), ]; s$ref_id <- sample(1:15, 60, TRUE)
  p <- profile_domain("T", dm, s, c(3L, 6L), des, 40, 100L)
  expect_equal(p$k, c(3L, 6L)); expect_true(all(p$iqr_log > 0)); expect_true(all(p$fit_ok))
  base <- iqr_log_by_k(dm, c(3L, 6L))
  r <- to_base_scale(p, base)
  expect_equal(r$beta_base, r$beta_iqr, tolerance = 1e-3)         # own IQR = baseline IQR here
  fr <- fullrefit_k(dm, s, 6L, des, base[["6"]], 40, 20L)
  expect_equal(fr$n_boot_ok, 20L); expect_equal(fr$beta_iqr, fr$beta_base, tolerance = 1e-3)
  expect_true(fr$lo <= fr$beta_iqr + 1e-6 || fr$lo <= fr$hi)
})