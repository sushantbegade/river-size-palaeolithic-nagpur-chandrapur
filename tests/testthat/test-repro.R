# tests/testthat/test-repro.R
source(here::here("R", "config.R")); source(here::here("R", "utils.R")); source(here::here("R", "repro.R"))

test_that("coordinate scan flags UTM and lon/lat pairs, reports line numbers only, and skips ordinary numbers", {
  f <- tempfile(fileext = ".csv")
  e <- paste0("3345", "00.7"); n <- paste0("23105", "00.2"); lo <- paste0("79.12", "3456"); la <- paste0("21.65", "4321")   # built at run time so this file holds no coordinate-like literals
  writeLines(c("a,b,c", "n1,1234,5.5", paste0("n2,", e, ",", n), paste0("n3,x,", lo, " and ", la), "md5,a1b2c3d4e5f6,100000", "ok,1500000,2500000"), f)
  h <- repro_scan_coords(f)
  expect_length(h, 1); expect_match(h, "line\\(s\\) 3,4"); expect_false(grepl("3345", h))
})
test_that("public 1 km grid check accepts cell centres and rejects exact positions", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(x_c = c(1500, 2500), y_c = c(3500, 4500)), f, row.names = FALSE); expect_true(repro_public_grid(f)$ok)
  utils::write.csv(data.frame(x_c = c(1512, 2500), y_c = c(3500, 4500)), f, row.names = FALSE); expect_false(repro_public_grid(f)$ok)
})
test_that("repro_check_repo flags forbidden types, neutral-named coordinate columns, UTM pairs and a bad public file", {
  skip_if(!nzchar(Sys.which("git")), "git not available")
  d <- tempfile("repo"); dir.create(d); dir.create(file.path(d, "data", "public"), recursive = TRUE); dir.create(file.path(d, "output"))
  old <- setwd(d); on.exit(setwd(old))
  system2("git", c("init", "-q")); system2("git", c("config", "user.email", "t@t")); system2("git", c("config", "user.name", "t"))
  writeLines("x,y\n1,2", "output/summary.csv")                                 # header x,y
  writeLines(paste0("e,n,v\n3345", "00,23105", "00,1\n"), "output/pair.csv")      # neutral names, UTM pair
  writeLines("a,b\n1,2", "output/fine.csv")
  saveRDS(1, "output/obj.rds")
  utils::write.csv(data.frame(x_c = 1512, y_c = 3500), file.path("data", "public", "sites_public_1km.csv"), row.names = FALSE)
  system2("git", c("add", "-A"))
  r <- repro_check_repo(d)
  txt <- paste(r$problems, collapse = "\n")
  expect_match(txt, "forbidden type: output/obj.rds"); expect_match(txt, "coordinate-like column\\(s\\) \\[x,y\\] in output/summary.csv")
  expect_match(txt, "possible exact coordinates.*output/pair.csv"); expect_match(txt, "not verified on the 1 km grid")
  expect_false(grepl("fine.csv", txt))
})
test_that("freeze verification classifies unchanged, changed and missing files", {
  d <- tempfile("fz"); dir.create(file.path(d, "output", "logs"), recursive = TRUE); dir.create(file.path(d, "R"), recursive = TRUE)
  writeLines("a", file.path(d, "R", "a.R")); writeLines("b", file.path(d, "R", "b.R")); writeLines("c", file.path(d, "R", "c.R"))
  fs <- file.path(d, "R", c("a.R", "b.R", "c.R"))
  utils::write.csv(data.frame(file = c("R/a.R", "R/b.R", "R/c.R", "R/gone.R"), md5 = c(unname(tools::md5sum(fs[1:2])), "zzz", "yyy")), file.path(d, "output", "logs", "DESIGN_FREEZE_manifest.csv"), row.names = FALSE)
  v <- repro_verify_freeze(d, expect_md5 = "x")
  expect_equal(v$table$status, c("unchanged", "unchanged", "changed", "missing")); expect_false(v$md5_ok)
})
test_that("package scan finds library(), require() and pkg:: calls", {
  d <- tempfile("pk"); dir.create(file.path(d, "R"), recursive = TRUE); writeLines(c("library(sf); library(dplyr)", "x <- terra::rast", "# readr::read_csv", "suppressPackageStartupMessages({ library(here) })"), file.path(d, "R", "a.R"))
  expect_equal(repro_pkgs(d), c("dplyr", "here", "sf", "terra"))
})