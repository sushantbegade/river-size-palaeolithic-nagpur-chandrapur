# scripts/08b_osm_diagnostic.R
suppressPackageStartupMessages({ library(here); library(readr); library(sf) })
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "external_data.R"))
source(here::here("R", "detectability.R"))

message("== S08b OSM road diagnostic: raw XML counts vs GDAL counts, per tile ==")
files <- list.files(EXT("osm"), pattern = "^roads_.*\\.osm$", full.names = TRUE)
if (!length(files)) stop("No road tiles in data/external/osm.", call. = FALSE)

count_xml <- function(f, chunk = 200000L) {
  con <- file(f, "r", encoding = "UTF-8"); on.exit(close(con))
  n <- c(node = 0, way = 0, hwy = 0, major = 0)
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE)
    if (!length(x)) break
    n["node"] <- n["node"] + sum(grepl("<node id=", x, fixed = TRUE))
    n["way"]  <- n["way"]  + sum(grepl("<way id=", x, fixed = TRUE))
    h <- x[grepl('k="highway"', x, fixed = TRUE)]
    n["hwy"] <- n["hwy"] + length(h)
    n["major"] <- n["major"] + sum(grepl(paste0('v="(', paste(HWY_MAJOR, collapse = "|"), ')"'), h))
  }
  n
}

out <- do.call(rbind, lapply(seq_along(files), function(i) {
  f <- files[i]
  x <- count_xml(f)
  g <- read_osm_layer(f, "lines")
  gh <- if (is.null(g)) NA else sum(g$highway %in% c(HWY_MAJOR, HWY_MINOR), na.rm = TRUE)
  gm <- if (is.null(g)) NA else sum(g$highway %in% HWY_MAJOR, na.rm = TRUE)
  gl <- if (is.null(g)) NA else nrow(g)
  message(sprintf("  [%2d/%2d] %s", i, length(files), basename(f)))
  data.frame(tile = sub("^roads_|\\.osm$", "", basename(f)), mb = round(file.size(f) / 1e6, 1),
             xml_nodes = x[["node"]], xml_ways = x[["way"]], xml_hwy_tags = x[["hwy"]], xml_major = x[["major"]],
             gdal_lines = gl, gdal_hwy = gh, gdal_major = gm)
}))
out$gdal_share_of_xml_hwy <- round(out$gdal_hwy / pmax(out$xml_hwy_tags, 1), 3)
write_csv_safe(out, P("output", "tables", "S08b_osm_road_diagnostic.csv"))

cat("\n========== S08b SUMMARY ==========\n")
print(out, row.names = FALSE)
cat("\nTOTAL xml ways:", sum(out$xml_ways), "| xml highway tags:", sum(out$xml_hwy_tags),
    "| GDAL highway lines:", sum(out$gdal_hwy, na.rm = TRUE), "\n")
cat("\nFirst 12 lines of the largest tile (for format check):\n")
big <- files[which.max(file.size(files))]
print(readLines(big, n = 12, warn = FALSE))
cat("\nDONE S08b. Paste this console output back.\n")