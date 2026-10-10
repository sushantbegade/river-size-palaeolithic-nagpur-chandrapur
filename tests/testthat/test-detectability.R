# tests/testthat/test-detectability.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "detectability.R"))

test_that("osm_tiles cover the bbox without gaps", {
  bb <- c(xmin = 78.1, ymin = 19.2, xmax = 80.4, ymax = 21.6)
  tl <- osm_tiles(bb, 0.5)
  expect_true(min(tl$w) <= bb[["xmin"]] && max(tl$e) >= bb[["xmax"]])
  expect_true(min(tl$s) <= bb[["ymin"]] && max(tl$n) >= bb[["ymax"]])
  expect_false(anyDuplicated(tl$tile) > 0)
})

test_that("wc_tile_ids follow lower-left 3-degree naming", {
  expect_equal(wc_tile_ids(c(xmin = 78.3, ymin = 19.2, xmax = 80.3, ymax = 21.6)), c("N18E078", "N21E078"))
  expect_equal(wc_tile_ids(c(xmin = 78.1, ymin = 19.0, xmax = 80.9, ymax = 20.9)), "N18E078")
})

test_that("highway and place classes", {
  expect_equal(classify_highway(c("primary", "residential", "footway", NA)), c("major", "minor", NA, NA))
})

test_that("osm_xml_ok rejects html, remark and truncated answers", {
  tf <- tempfile(fileext = ".osm")
  writeLines("<html><body>Too Many Requests</body></html>", tf); expect_false(osm_xml_ok(tf))
  writeLines(c("<?xml version=\"1.0\"?>", "<osm version=\"0.6\">", "<remark> runtime error </remark>", "</osm>"), tf)
  expect_false(osm_xml_ok(tf))
  writeLines(c("<?xml version=\"1.0\"?>", "<osm version=\"0.6\">", "<node id=\"1\" lat=\"1\" lon=\"1\"/>"), tf)
  expect_false(osm_xml_ok(tf))
  writeLines(c("<?xml version=\"1.0\"?>", "<osm version=\"0.6\">", "<node id=\"1\" lat=\"1\" lon=\"1\"/>", "</osm>"), tf)
  expect_true(osm_xml_ok(tf))
})

test_that("land-cover group fractions are exact on a toy raster", {
  r <- terra::rast(nrows = 110, ncols = 110, xmin = 79, xmax = 79 + 110 / 12000,
                   ymin = 21, ymax = 21 + 110 / 12000, crs = "EPSG:4326")
  terra::values(r) <- rep(c(40, 10), length.out = terra::ncell(r))
  tf <- tempfile(fileext = ".tif"); terra::writeRaster(r, tf, datatype = "INT1U")
  fr <- wc_group_fractions(tf, c(xmin = 79, ymin = 21, xmax = 79 + 110 / 12000, ymax = 21 + 110 / 12000), margin = 0)
  v <- terra::values(fr)
  expect_true(all(abs(rowSums(v, na.rm = TRUE) - 1) < 1e-9))
  expect_true(all(abs(v[, "frac_crop"] - 0.5) < 0.05))
})