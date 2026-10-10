suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf); library(terra)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "persistence.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S03 persistence index ==")

gpkg <- P("data", "processed", "hydro_prep.gpkg")
if (!file.exists(gpkg)) stop("Run scripts/02_prepare_hydrography.R first.", call. = FALSE)
s2dir <- DATA("Sentinel 2A Multispectral Data")
f_b03 <- file.path(s2dir, "B03_Nagpur_Chandrapur.tif")
f_b11 <- file.path(s2dir, "B11_Nagpur_Chandrapur.tif")
for (f in c(f_b03, f_b11)) if (!file.exists(f)) stop("Missing raster: ", f, call. = FALSE)

b03 <- terra::rast(f_b03); b11 <- terra::rast(f_b11)
code03 <- terra::crs(b03, describe = TRUE)$code
code11 <- terra::crs(b11, describe = TRUE)$code
if (!identical(as.character(code03), as.character(CRS_ANALYSIS)) ||
    !identical(as.character(code11), as.character(CRS_ANALYSIS))) {
  stop("Sentinel rasters are not EPSG:", CRS_ANALYSIS, " (found ", code03, " / ", code11, ").",
       call. = FALSE)
}
message("B03 res ", paste(terra::res(b03), collapse = "x"), " m | B11 res ",
        paste(terra::res(b11), collapse = "x"), " m")

rv   <- sf::st_read(gpkg, layer = "rivers", quiet = TRUE)
core <- sf::st_read(gpkg, layer = "removed_core", quiet = TRUE)
bnd  <- sf::st_read(gpkg, layer = "boundary", quiet = TRUE)
bu   <- sf::st_union(sf::st_geometry(bnd))

lg <- rv[!is.na(rv$ordsh) & rv$ordsh >= MIN_ORD_PI, ]
message("reaches sampled (ordsh >= ", MIN_ORD_PI, "): ", nrow(lg), " | ",
        round(sum(lg$length_km), 1), " km")

pts <- reach_points(lg)
pts$in_removed <- lengths(sf::st_intersects(
  sf::st_as_sf(pts, coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE), core)) > 0
message("points: ", nrow(pts), " | inside removed areas: ", sum(pts$in_removed))
message("extracting Sentinel values along reaches (a few minutes)...")
nb  <- neighbourhood_mndwi(b03, b11, as.matrix(pts[, c("x", "y")]))
pts <- add_metrics(pts, nb)
pts$valid <- !is.na(pts$mndwi_max40) & !pts$cloudy & !pts$in_removed
pts$og    <- order_group(pts$ordsh)

message("baseline land points...")
bxy <- baseline_points(bu, core, lg)
bnb <- neighbourhood_mndwi(b03, b11, as.matrix(bxy))
base <- add_metrics(bxy, bnb)
base$valid <- !is.na(base$mndwi_max40) & !base$cloudy
message("baseline points: ", nrow(base), " (valid ", sum(base$valid), ")")

sub <- subsegment_table(pts, thr = 0)

fp <- fp_table(pts, base)
write_csv_safe(fp, P("output", "tables", "S03_baseline_fp_by_threshold.csv"))

auc <- data.frame(
  group = c("order5", "order6to8", "order9plus"),
  auc_vs_baseline = vapply(c("order5", "order6to8", "order9plus"), function(g)
    round(auc_mw(pts$mndwi_max40[pts$valid & pts$og == g], base$mndwi_max40[base$valid]), 3), numeric(1)),
  n_points = vapply(c("order5", "order6to8", "order9plus"), function(g) sum(pts$valid & pts$og == g), numeric(1)))
write_csv_safe(auc, P("output", "tables", "S03_detectability_auc.csv"))

by_order <- pts |>
  dplyr::filter(valid) |>
  dplyr::group_by(ordsh) |>
  dplyr::summarise(n_points = dplyr::n(), length_km = round(dplyr::n() * SPACING_M / 1000, 1),
                   median_mndwi40 = round(stats::median(mndwi_max40), 3),
                   share_gt0 = round(mean(mndwi_max40 > 0), 3),
                   share_centre_gt0 = round(mean(mndwi_c > 0, na.rm = TRUE), 3), .groups = "drop")
write_csv_safe(by_order, P("output", "tables", "S03_pi_by_order.csv"))

cls_len <- sub |>
  dplyr::filter(usable) |>
  dplyr::mutate(og = order_group(ordsh)) |>
  dplyr::group_by(og, pi_class) |>
  dplyr::summarise(n_sub = dplyr::n(), length_km = round(sum(length_km_valid), 1), .groups = "drop")
write_csv_safe(cls_len, P("output", "tables", "S03_provisional_class_length.csv"))

by_river <- sub |>
  dplyr::filter(usable, !is.na(ordsh), ordsh >= 6L) |>
  dplyr::mutate(rivname = dplyr::coalesce(rivname, "<unnamed>")) |>
  dplyr::group_by(rivname) |>
  dplyr::summarise(max_ordsh = max(ordsh), n_sub = dplyr::n(),
                   length_km = round(sum(length_km_valid), 1),
                   median_pi = round(stats::median(pi), 3),
                   len_P_cand_km = round(sum(length_km_valid[pi_class == "P_candidate"]), 1),
                   len_I_km = round(sum(length_km_valid[pi_class == "I"]), 1),
                   len_S_km = round(sum(length_km_valid[pi_class == "S"]), 1), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(max_ordsh), dplyr::desc(length_km))
write_csv_safe(by_river, P("output", "tables", "S03_pi_by_river_order6plus.csv"))

saveRDS(list(points = pts, baseline = base), P("data", "processed", "pi_points.rds"))
write_csv_safe(sub, P("data", "processed", "pi_subsegments.csv"))

meta <- data.frame(item = c("B03_res_m", "B11_res_m", "spacing_m", "chunk_m", "neigh_radii_m",
                            "B03_cloud_dn", "PI_P", "PI_I", "n_points", "n_valid", "n_baseline_valid"),
                   value = c(paste(terra::res(b03), collapse = "x"), paste(terra::res(b11), collapse = "x"),
                             SPACING_M, CHUNK_M, paste(NEIGH_RADII, collapse = ","), B03_CLOUD,
                             PI_P, PI_I, nrow(pts), sum(pts$valid), sum(base$valid)))
write_csv_safe(meta, P("output", "logs", "S03_parameters.csv"))
message("computing checksums of B03 and B11 (may take a minute)...")
write_csv_safe(file_checksums(c(f_b03, f_b11)), P("output", "logs", "S03_input_checksums.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S03_sessionInfo.txt"))

cat("\n========== S03 SUMMARY ==========\n")
cat("Points:", nrow(pts), "| valid:", sum(pts$valid), "| inside removed:", sum(pts$in_removed),
    "| cloudy:", sum(pts$cloudy), "| NA:", sum(is.na(pts$mndwi_max40)), "\n")
cat("Baseline land points valid:", sum(base$valid), "\n")
cat("\nFalse-positive / water share by group and threshold (max within 40 m | centre pixel):\n")
print(as.data.frame(fp), row.names = FALSE)
cat("\nDetectability (AUC vs baseline land; 0.5 = no signal, 1 = perfect):\n")
print(auc, row.names = FALSE)
cat("\nBy stream order:\n"); print(as.data.frame(by_order), row.names = FALSE)
cat("\nProvisional class length (km) by order group (PROVISIONAL thresholds 0.5 / 0.15):\n")
print(as.data.frame(cls_len), row.names = FALSE)
cat("\nLarge rivers (order >= 6), provisional class length:\n")
print(as.data.frame(utils::head(by_river, 30)), row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S03. Paste this console output back.\n")