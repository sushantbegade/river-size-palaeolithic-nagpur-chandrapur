# tests for R/figures.R (S18). Run from project root: testthat::test_file("tests/testthat/test-figures.R")
source(here::here("R", "config.R")); source(here::here("R", "utils.R")); source(here::here("R", "figures.R"))
need <- c("S11_profile", "S11_primary_k6", "S12_sample_sets_profile", "S12_tiers_profile", "S13_site_type_profile", "S13_decomposition_k6", "S14_cv_summary",
          "S15_cluster_definitions", "S15_conditional_profile", "S16_variants_profile", "S17_distance_shape")
tb <- setNames(lapply(need, function(n) rd_csv(paste0(n, ".csv"))), need)

test_that("near-zero limits keep their sign and precision", {
  expect_equal(fig_fmt(c(-0.002, 0.003, 0.014, 0.123, 0)), c("-0.002", "0.003", "0.01", "0.12", "0.00"))
  expect_equal(fig_est_ci(0.596, -0.002, 1.311), "0.60 (-0.002 to 1.31)")
  expect_equal(fig_est_ci(0.5, NA, NA), "0.50")
})
test_that("Fig2 values equal the saved tables and k8 lower limit is -0.002", {
  d <- fig2_data(tb); a8 <- d[d$panel == "a" & d$x == 8, ]
  expect_equal(a8$lo, -0.002); expect_equal(d$y[d$panel == "a" & d$x == 6], 1.236)
  expect_equal(sort(unique(d$series)), c("S1", "S2", "S3"))
})
test_that("Fig5 has no S2 rows for tiers or full-refit blocks (none exist in the tables)", {
  d <- fig5_data(tb); s2 <- d[d$panel == "b", ]
  expect_false(any(grepl("tier", s2$x_label))); expect_false(any(grepl("blocks", s2$x_label)))
  expect_equal(d$y[d$panel == "a" & d$x_label == "5 km survey domain"], 1.002)
})
test_that("Fig4 uses +-1 SE from se_delta_sum and keeps 40 and 20 km", {
  d <- fig4_data(tb); r <- d[d$series == "S1 40 km" & d$x_label == "k6", ]
  expect_equal(r$y, 28.92); expect_equal(r$hi - r$y, 10.61, tolerance = 1e-9); expect_true(all(c("S1 20 km", "S3 40 km") %in% d$series))
})
test_that("audit passes on written FigSource and catches a tampered value", {
  d <- rbind(fig2_data(tb), fig5_data(tb)); f <- tempfile(fileext = ".csv")
  w <- d; names(w)[names(w) == "x"] <- "x_pos"; names(w)[names(w) == "y"] <- "y_val"; utils::write.csv(w, f, row.names = FALSE)
  expect_equal(fig_audit(f, rd_csv)$n_bad, 0)
  w$y_val[5] <- w$y_val[5] + 0.01; utils::write.csv(w, f, row.names = FALSE)
  expect_equal(fig_audit(f, rd_csv)$n_bad, 1)
})
test_that("map refuses a public file that is not on the 1 km grid", {
  skip_if(Sys.getenv("PALEO_TEST_MAP") != "1", "rewrites data/public temporarily; set PALEO_TEST_MAP=1 to run on a project COPY")
  fp <- P("data", "public", "sites_public_1km.csv"); keep <- utils::read.csv(fp)
  bad <- keep; bad$x_c[1] <- bad$x_c[1] + 37; utils::write.csv(bad, fp, row.names = FALSE)
  on.exit(utils::write.csv(keep, fp, row.names = FALSE))
  expect_error(fig1_load(CRS_ANALYSIS), "1 km")
})
test_that("FigSource Fig1 has counts only and no coordinate columns", {
  m <- fig1_load(CRS_ANALYSIS); d <- fig1_data(m)
  expect_false(any(grepl("^(lon|lat|longitude|latitude|easting|northing|x_c$|y_c$)", names(d), ignore.case = TRUE))); expect_equal(d$y[d$series == "total points drawn"], nrow(m$pts))
})
# rd_csv(): reads one saved table from output/tables (searches sub-folders if not at top level); NULL if absent.
# Defined here only if your utils.R does not already define it (it was probably in your old R/figures.R).
if (!exists("rd_csv", mode = "function")) {
  rd_csv <- function(name) {
    f <- P("output", "tables", name)
    if (!file.exists(f)) { h <- list.files(P("output", "tables"), pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"), recursive = TRUE, full.names = TRUE)
    h <- h[!grepl("/paper/", h)]; if (!length(h)) return(NULL); f <- h[1] }
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  }
}