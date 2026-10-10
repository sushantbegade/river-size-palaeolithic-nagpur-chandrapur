# Regression test against the counts measured in the data audit (Oct 2026).
# Skipped automatically if S01 has not been run on this machine.
suppressPackageStartupMessages(library(testthat))
f <- here::here("data", "processed", "sites_clean.rds")

test_that("cleaned sites match audit counts", {
  skip_if_not(file.exists(f), "run scripts/01_clean_sites.R first")
  d <- readRDS(f)
  expect_equal(nrow(d), 197L)
  expect_equal(sum(d$S1_openair), 93L)
  expect_equal(sum(d$S3_riverbed), 88L)
  expect_equal(as.integer(table(d$coord_tier)[c("A_field_verified", "B_author_reported", "C_map_based")]),
               c(64L, 73L, 60L))
  expect_equal(as.integer(table(d$period_class)[c("LP", "MP", "UP", "multi")]),
               c(48L, 62L, 27L, 60L))
  expect_equal(c(sum(d$LP), sum(d$MP), sum(d$UP)), c(75L, 114L, 79L))
  expect_false(anyDuplicated(d$site_id) > 0)
})
