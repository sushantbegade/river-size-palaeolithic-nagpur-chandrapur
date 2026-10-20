# tests for the documented-change support in R/repro.R (added by scripts/patch_post_freeze.R); synthetic files only
suppressWarnings(suppressPackageStartupMessages(library(testthat)))
source(here::here("R", "repro.R"))

mk_proj <- function() {
  d <- tempfile("pf"); dir.create(file.path(d, "output", "logs"), recursive = TRUE); dir.create(file.path(d, "R")); dir.create(file.path(d, "docs"))
  for (n in c("a", "b", "c")) writeLines(n, file.path(d, "R", paste0(n, ".R")))
  fs <- file.path(d, "R", c("a.R", "b.R", "c.R"))
  utils::write.csv(data.frame(file = c("R/a.R", "R/b.R", "R/c.R"), md5 = unname(tools::md5sum(fs))), file.path(d, "output", "logs", "DESIGN_FREEZE_manifest.csv"), row.names = FALSE)
  d
}

test_that("a changed file is 'changed' without a change log and 'changed (documented)' when the log lists its current MD5", {
  d <- mk_proj()
  writeLines("a edited", file.path(d, "R", "a.R")); writeLines("b edited", file.path(d, "R", "b.R"))
  v0 <- repro_verify_freeze(d, expect_md5 = "x")
  expect_equal(v0$table$status, c("changed", "changed", "unchanged"))
  pc <- data.frame(file = c("R/a.R", "R/b.R"), md5_before = "", md5_after = c(unname(tools::md5sum(file.path(d, "R", "a.R"))), "not-the-current-md5"))
  utils::write.csv(pc, file.path(d, "docs", "post_freeze_changes.csv"), row.names = FALSE)
  v1 <- repro_verify_freeze(d, expect_md5 = "x")
  expect_equal(v1$table$status, c("changed (documented)", "changed", "unchanged"))
})

test_that("a further edit after the logged change makes the file 'changed' again; missing files stay missing", {
  d <- mk_proj(); writeLines("a edited", file.path(d, "R", "a.R"))
  utils::write.csv(data.frame(file = "R/a.R", md5_after = unname(tools::md5sum(file.path(d, "R", "a.R")))), file.path(d, "docs", "post_freeze_changes.csv"), row.names = FALSE)
  expect_equal(repro_verify_freeze(d, expect_md5 = "x")$table$status[1], "changed (documented)")
  writeLines("a edited again", file.path(d, "R", "a.R")); file.remove(file.path(d, "R", "c.R"))
  v <- repro_verify_freeze(d, expect_md5 = "x")
  expect_equal(v$table$status, c("changed", "unchanged", "missing"))
})

test_that("a malformed change log is ignored rather than breaking the check", {
  d <- mk_proj(); writeLines("a edited", file.path(d, "R", "a.R"))
  utils::write.csv(data.frame(foo = 1), file.path(d, "docs", "post_freeze_changes.csv"), row.names = FALSE)
  expect_equal(repro_verify_freeze(d, expect_md5 = "x")$table$status[1], "changed")
})