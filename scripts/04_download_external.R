suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "external_data.R"))

t0 <- Sys.time()
message("== S04 download external data ==")
gpkg_h <- P("data", "processed", "hydro_prep.gpkg")
if (!file.exists(gpkg_h)) stop("Run scripts/02_prepare_hydrography.R first.", call. = FALSE)
bnd <- sf::st_read(gpkg_h, layer = "boundary", quiet = TRUE)
bu  <- sf::st_union(sf::st_geometry(bnd))

url_hr <- "https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_as_shp.zip"
zip_hr <- EXT("HydroRIVERS_v10_as_shp.zip")
download_if_missing(url_hr, zip_hr)
unzip_once(zip_hr, EXT("HydroRIVERS_v10_as"))
hr <- read_region(find_shp(EXT("HydroRIVERS_v10_as"), "HydroRIVERS_v10_as"), "HydroRIVERS")
require_fields(hr, c("HYRIV_ID", "NEXT_DOWN", "MAIN_RIV", "LENGTH_KM", "UPLAND_SKM",
                     "DIS_AV_CMS", "ORD_STRA"), "HydroRIVERS")
hr$in_study <- lengths(sf::st_intersects(hr, bu)) > 0

zip_gdw <- EXT("GDW_v1_0_shp.zip")
gdw_md5 <- NA_character_; gdw_url <- "figshare article 25988293"
if (!is_valid_zip(zip_gdw)) {
  info <- figshare_file(25988293, "shp\\.zip$")
  message("figshare file: ", info$name, " (", round(info$size / 1e6, 1), " MB)")
  gdw_md5 <- info$md5; gdw_url <- info$url
  download_if_missing(info$url, zip_gdw, md5 = info$md5)
}
unzip_once(zip_gdw, EXT("GDW_v1_0"))
shp_all <- list.files(EXT("GDW_v1_0"), pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
message("GDW shapefiles found: ", paste(basename(shp_all), collapse = ", "))
barr <- read_region(find_shp(EXT("GDW_v1_0"), "barrier"), "GDW barriers")
resv <- read_region(find_shp(EXT("GDW_v1_0"), "reservoir"), "GDW reservoirs")

gpkg <- P("data", "processed", "external_region.gpkg")
unlink(gpkg)
sf::st_write(hr,   gpkg, layer = "hydrorivers",    quiet = TRUE)
sf::st_write(barr, gpkg, layer = "gdw_barriers",   quiet = TRUE)
sf::st_write(resv, gpkg, layer = "gdw_reservoirs", quiet = TRUE)
message("wrote data/processed/external_region.gpkg")

hs <- sf::st_drop_geometry(hr)
hr_sum <- dplyr::summarise(dplyr::group_by(hs, in_study, ORD_STRA), n_reaches = dplyr::n(),
                           length_km = round(sum(LENGTH_KM), 1),
                           median_upland_km2 = round(stats::median(UPLAND_SKM), 0), .groups = "drop")
write_csv_safe(hr_sum, P("output", "tables", "S04_hydrorivers_summary.csv"))

near <- lengths(sf::st_is_within_distance(barr, bu, dist = 100000)) > 0
barr_near <- drop_coord_cols(sf::st_drop_geometry(barr[near, ]))
write_csv_safe(barr_near, P("output", "tables", "S04_gdw_barriers_within100km.csv"))

src <- data.frame(
  item = c("HydroRIVERS v10 Asia shapefile", "Global Dam Watch v1.0 shapefile"),
  source = c(url_hr, gdw_url),
  file = c(basename(zip_hr), basename(zip_gdw)),
  bytes = c(file.size(zip_hr), file.size(zip_gdw)),
  md5 = c(unname(tools::md5sum(zip_hr)), unname(tools::md5sum(zip_gdw))),
  published_md5 = c(NA, gdw_md5),
  downloaded = format(Sys.Date()),
  licence_note = c("HydroSHEDS licence: attribution required, see hydrosheds.org",
                   "GDW v1.0: check licence on figshare/globaldamwatch.org"))
write_csv_safe(src, P("output", "logs", "S04_external_sources.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S04_sessionInfo.txt"))

cat("\n========== S04 SUMMARY ==========\n")
cat("HydroRIVERS reaches in region:", nrow(hr), "| intersecting districts:", sum(hr$in_study), "\n")
cat("HydroRIVERS fields:", paste(names(hr)[names(hr) != attr(hr, "sf_column")], collapse = ", "), "\n")
cat("\nReaches intersecting districts, by Strahler order (HydroRIVERS ORD_STRA):\n")
print(as.data.frame(hr_sum[hr_sum$in_study, ]), row.names = FALSE)
cat("\nGDW barriers in region:", nrow(barr), "| within 100 km of districts:", sum(near), "\n")
cat("GDW barrier fields:", paste(names(barr)[names(barr) != attr(barr, "sf_column")], collapse = ", "), "\n")
cat("GDW reservoir polygons in region:", nrow(resv), "\n")
cat("GDW reservoir fields:", paste(names(resv)[names(resv) != attr(resv, "sf_column")], collapse = ", "), "\n")
cat("\nBarriers within 100 km of the districts (first 60 rows, name/year/volume-like columns):\n")
cols <- grep("name|year|cap|vol|area|purp|type|height|river", names(barr_near), ignore.case = TRUE, value = TRUE)
if (!length(cols)) cols <- names(barr_near)[seq_len(min(10, ncol(barr_near)))]
print(utils::head(as.data.frame(barr_near[, cols, drop = FALSE]), 60), row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S04. Paste this console output back.\n")