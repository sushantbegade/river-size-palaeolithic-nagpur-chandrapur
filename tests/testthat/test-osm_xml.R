# tests/testthat/test-osm_xml.R
source(here::here("R", "config.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "detectability.R"))
source(here::here("R", "osm_xml.R"))

osm_fixture <- function() {
  tf <- tempfile(fileext = ".osm")
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<osm version="0.6" generator="Overpass API">',
    '  <node id="8000000001" lat="21.0000000" lon="79.0000000"/>',
    '  <node id="8000000002" lat="21.0000000" lon="79.0100000"/>',
    '  <node id="8000000003" lat="21.0100000" lon="79.0100000">',
    '    <tag k="highway" v="crossing"/>',
    '    <tag k="place" v="village"/>',
    '    <tag k="name" v="X"/>',
    '  </node>',
    '  <node id="8000000004" lat="21.0200000" lon="79.0200000">',
    '    <tag k="place" v="hamlet"/>',
    '  </node>',
    '  <node id="8000000005" lat="21.0300000" lon="79.0300000">',
    '    <tag k="place" v="suburb"/>',
    '  </node>',
    '  <way id="11">',
    '    <nd ref="8000000001"/>', '    <nd ref="8000000002"/>', '    <nd ref="8000000003"/>',
    '    <tag k="name" v="A"/>', '    <tag k="highway" v="primary"/>',
    '  </way>',
    '  <way id="12">',
    '    <nd ref="8000000001"/>', '    <nd ref="8000000003"/>',
    '    <tag k="highway" v="footway"/>',
    '  </way>',
    '  <way id="13">',
    '    <nd ref="8000000002"/>', '    <nd ref="9999999999"/>',
    '    <tag k="highway" v="tertiary"/>',
    '  </way>',
    '</osm>'), tf)
  tf
}

test_that("roads: keeps only wanted classes, ignores node tags, drops ways with <2 known nodes", {
  r <- parse_osm_roads_xml(osm_fixture())
  expect_equal(r$osm_id, "11")
  expect_equal(r$highway, "primary")
  expect_equal(nrow(sf::st_coordinates(r)), 3L)
})

test_that("places: keeps settlement classes only, one row per tagged node", {
  p <- parse_osm_places_xml(osm_fixture())
  expect_equal(sort(p$place), c("hamlet", "village"))
  expect_equal(sort(p$osm_id), c("8000000003", "8000000004"))
})

test_that("tile reader removes duplicate ways", {
  f <- osm_fixture()
  expect_equal(nrow(read_osm_xml_tiles(c(f, f), parse_osm_roads_xml)), 1L)
})