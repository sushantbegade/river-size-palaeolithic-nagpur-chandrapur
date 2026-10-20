# tests for R/source_struct.R (synthetic data only)
suppressWarnings(suppressPackageStartupMessages(library(testthat)))
for (f in c("config", "utils", "quadrature", "stepzero", "stepzero_b", "stepzero_c", "model_primary", "diagnostics", "domain_sens", "source_struct"))
  source(here::here("R", paste0(f, ".R")))

test_that("family_table counts strings and sites per first-named author", {
  sc <- data.frame(ref_clean = c("A, 1960", "A, 1961; B, 2000", "B, 2000", "C, 1970", "C, 1970"),
                   S1_openair = c(TRUE, TRUE, FALSE, TRUE, FALSE), S3_riverbed = c(FALSE, FALSE, TRUE, FALSE, TRUE))
  ft <- family_table(sc, ref_family(sc$ref_clean))
  expect_equal(ft$family[ft$n_S1 == max(ft$n_S1)], "A")
  expect_equal(ft$n_strings[ft$family == "A"], 2L); expect_equal(ft$n_S1[ft$family == "A"], 2L)
  expect_equal(sum(ft$n_sites), 5L); expect_equal(sum(ft$share_of_S1), 1, tolerance = 0.01)
})

test_that("conditional logit recovers a within-stratum effect where the unconditional fixed-effect logit is biased", {
  set.seed(2)
  ns <- 150; sz <- 4; st <- rep(seq_len(ns), each = sz); n <- length(st)
  z <- rnorm(n); a <- rnorm(ns, 0, 1.5)[st]; y <- rbinom(n, 1, plogis(a + 1 * z))
  c1 <- fe_logit(y, z, st); u1 <- fe_logit(y, z, st, conditional = FALSE)
  expect_true(c1$ok); expect_lt(abs(c1$coef - 1), 0.25)
  expect_gt(u1$coef, c1$coef + 0.15)                      # unconditional estimate is pushed away from 0 with strata of 4
})

test_that("a single stratum gives ordinary logistic regression; nuisance columns are accepted", {
  set.seed(4); n <- 300; z <- rnorm(n); x <- cbind(a = rnorm(n), b = rnorm(n)); y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x[, 1]))
  f0 <- fe_logit(y, z, rep(1, n), x); g <- stats::glm(y ~ z + x[, 1] + x[, 2], family = binomial())
  expect_equal(f0$coef, unname(coef(g)["z"]), tolerance = 1e-6)
  st <- rep(1:30, each = 10); f1 <- fe_logit(y, z, st, x)
  expect_true(f1$ok); expect_true(is.finite(f1$coef))
})

test_that("boot_logit and label_row return finite intervals and respect cluster resampling", {
  set.seed(6)
  ns <- 40; st <- rep(seq_len(ns), each = 5); z <- rnorm(length(st)); y <- rbinom(length(st), 1, plogis(rnorm(ns, 0, 1)[st] - 1 * z))
  r <- label_row("t", y, z, st, st, iqr_per_sd = 1.3, fe = TRUE, B = 80L)
  expect_equal(r$n_strata, ns); expect_true(r$lo < r$log_odds_per_IQR + 1e-9 || r$lo < r$hi)
  expect_gt(r$n_boot_ok, 60L); expect_lt(r$log_odds_per_IQR, 0)
  expect_equal(r$n_S1 + r$n_S3, length(y))
})

test_that("loo_full returns one row per group, sorted by absolute change, with finite estimates", {
  set.seed(8)
  dm <- data.frame(x = runif(1500, 3e5, 3.4e5), y = runif(1500, 2.25e6, 2.29e6), elev_m = 300 + rnorm(1500, 0, 20), slope_deg = runif(1500, .2, 6),
                   gm4 = sample(c("pediplain", "plateau_hills", "alluvial"), 1500, TRUE, prob = c(.6, .3, .1)), frac_crop = runif(1500), frac_built = runif(1500, 0, .1),
                   d_road_major = runif(1500, 50, 9000), d_settle = runif(1500, 50, 6000), d_k6 = runif(1500, 20, 30000), area_km2 = 0.25)
  des <- design_from_domain(dm); ex <- exposure_domain(dm, 6); Xq <- cbind(int = 1, des$dom, z = ex$z_dom)
  s <- dm[sample(1500, 60), ]; s$grp <- sample(letters[1:6], 60, TRUE)
  r <- loo_full(s, "grp", des, ex, Xq, dm$area_km2)
  expect_equal(nrow(r), 6L); expect_true(all(is.finite(r$beta_iqr))); expect_false(is.unsorted(rev(abs(r$change_vs_full))))
  expect_equal(sum(r$n_left_out), 60L); expect_false(any(r$tied_beta_6dp))
})