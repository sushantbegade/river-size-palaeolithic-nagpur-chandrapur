# scripts/08_detectability.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf); library(terra)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "external_data.R"))
source(here::here("R", "detectability.R"))
source(here::here("R", "osm_xml.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S08 detectability covariates (OSM roads/places, ESA WorldCover) ==")

f_h <- P("data", "processed", "hydro_prep.gpkg")
f_q <- P("data", "processed", "quadrature_cov.rds")
for (f in c(f_h, f_q)) if (!file.exists(f)) stop("Missing: ", f, "\n  Run S02, S06, S07 first.", call. = FALSE)

bnd   <- sf::st_read(f_h, layer = "boundary", quiet = TRUE)
study <- sf::st_union(sf::st_geometry(bnd))
grid  <- readRDS(f_q)
need  <- c("q_id", "x", "y", "in_survey", "in_core", "elev_m", "slope_deg", "gm4",
           paste0("d_k", SIZE_K), "d_lineament")
if (!all(need %in% names(grid))) stop("quadrature_cov.rds lacks: ", paste(setdiff(need, names(grid)), collapse = ", "),
                                      "\n  Re-run S06 (patched) and S07.", call. = FALSE)
message("grid points: ", nrow(grid))

# ---- 1. OpenStreetMap download (Overpass, tiled, resumable) ---------------------------
bb_osm <- study_bbox_ll(study, OSM_BUFFER_DEG)
message("OSM bbox (lon/lat, study + ", OSM_BUFFER_DEG, " deg): ", paste(round(bb_osm, 3), collapse = " "))
osm_dir <- EXT("osm")
jobs <- rbind(cbind(kind = "roads",  osm_tiles(bb_osm, OSM_TILE_DEG_ROADS)),
              cbind(kind = "places", osm_tiles(bb_osm, OSM_TILE_DEG_PLACES)))
jobs$file   <- file.path(osm_dir, sprintf("%s_%s.osm", jobs$kind, jobs$tile))
jobs$status <- NA_character_
message("Overpass requests needed: ", nrow(jobs), " (", sum(jobs$kind == "roads"), " road tiles, ",
        sum(jobs$kind == "places"), " place tiles). Finished tiles are skipped on re-run.")
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, ]
  st <- overpass_fetch(osm_query(j$kind, j$s, j$w, j$n, j$e), j$file)
  jobs$status[i] <- st
  message(sprintf("  [%2d/%2d] %-6s %s : %s", i, nrow(jobs), j$kind, j$tile, st))
  if (st == "downloaded") Sys.sleep(OSM_PAUSE_S)
}
jobs$bytes <- ifelse(file.exists(jobs$file), file.size(jobs$file), NA)
write_csv_safe(jobs[, c("kind", "tile", "status", "bytes")], P("output", "tables", "S08_osm_tile_status.csv"))
if (any(jobs$status == "failed")) {
  stop(sum(jobs$status == "failed"), " OSM tile(s) failed: ",
       paste(paste(jobs$kind, jobs$tile)[jobs$status == "failed"], collapse = ", "),
       "\n  Wait a few minutes and source this script again (finished tiles are kept).",
       "\n  If road tiles keep failing, set OSM_TILE_DEG_ROADS <- 0.25 in R/detectability.R and delete data/external/osm.",
       call. = FALSE)
}

# ---- 2. read and classify OSM ------------------------------------------------------------
message("reading OSM tiles with the direct XML parser (a few minutes)...")
rd <- read_osm_xml_tiles(jobs$file[jobs$kind == "roads"],  parse_osm_roads_xml)
pl <- read_osm_xml_tiles(jobs$file[jobs$kind == "places"], parse_osm_places_xml)
if (is.null(rd) || is.null(pl)) stop("OSM read returned no roads or no places; inspect data/external/osm.", call. = FALSE)
message("parsed ways (unique): ", nrow(rd), " | parsed place nodes (unique): ", nrow(pl))
rd <- rd[!is.na(classify_highway(rd$highway)), ]
pl <- pl[pl$place %in% PLACE_SETTLE, ]
rd$hclass <- classify_highway(rd$highway)
rd <- sf::st_transform(rd, CRS_ANALYSIS)
pl <- sf::st_transform(pl, CRS_ANALYSIS)
osm_counts <- osm_count_table(rd, pl)
write_csv_safe(osm_counts, P("output", "tables", "S08_osm_counts.csv"))
message("roads: ", nrow(rd), " (major ", sum(rd$hclass == "major"), ") | places: ", nrow(pl))

grid$d_road_major <- nearest_dist(grid, sf::st_geometry(rd[rd$hclass == "major", ]))
grid$d_road_any   <- nearest_dist(grid, sf::st_geometry(rd))
grid$d_settle     <- nearest_dist(grid, sf::st_geometry(pl))
grid$d_town       <- nearest_dist(grid, sf::st_geometry(pl[pl$place %in% PLACE_TOWN, ]))

# ---- 3. ESA WorldCover 2021 v200 ---------------------------------------------------------
bb_study <- study_bbox_ll(study, 0)
wc_ids   <- wc_tile_ids(bb_study)
wc_paths <- EXT("worldcover", wc_file(wc_ids))
message("WorldCover tiles needed: ", paste(wc_ids, collapse = ", "))
for (k in seq_along(wc_ids)) download_tif_if_missing(wc_url(wc_ids[k]), wc_paths[k])

f_fr <- P("data", "processed", "worldcover_frac_100m.tif")
fr_names <- paste0("frac_", names(WC_GROUPS))
if (file.exists(f_fr) && identical(names(terra::rast(f_fr)), fr_names)) {
  message("using cached land-cover fractions: ", basename(f_fr), " (delete it to recompute)")
  fr <- terra::rast(f_fr)
} else {
  message("class-group fractions (10 m -> 100 m -> 0.5 km window)... 10-25 min on real tiles")
  fr <- wc_group_fractions(wc_paths, bb_study)
  terra::writeRaster(fr, f_fr, overwrite = TRUE)
}
ll <- sf::st_coordinates(sf::st_transform(sf::st_geometry(grid), 4326))
ex <- terra::extract(fr, ll)
for (v in names(fr)) grid[[v]] <- ex[[v]]

# ---- 4. diagnostics (domain = survey domain minus removed core; as S06/S07) ------------------
dd  <- sf::st_drop_geometry(grid)
dom <- dd[dd$in_survey & !dd$in_core, ]
dvars <- c("d_road_major", "d_road_any", "d_settle", "d_town")
fvars <- paste0("frac_", names(WC_GROUPS))
fsum  <- rowSums(dom[, fvars])
message("domain points: ", nrow(dom), " | NA by covariate: ",
        paste(paste0(c(dvars, fvars), "=", vapply(dom[, c(dvars, fvars)], function(z) sum(is.na(z)), numeric(1))),
              collapse = ", "))
message("sum of land-cover fractions: min ", round(min(fsum, na.rm = TRUE), 3),
        " | median ", round(stats::median(fsum, na.rm = TRUE), 3), " (1 - sum = other classes / no data)")

summ <- cov_summary(dom, c(dvars, fvars), km_vars = dvars)
write_csv_safe(summ, P("output", "tables", "S08_detect_summary.csv"))

size_vars <- paste0("d_k", SIZE_K)
X <- data.frame(elev = dom$elev_m, slope = dom$slope_deg,
                as.data.frame(lapply(dom[, c(size_vars, "d_lineament", dvars)], log_d)),
                dom[, fvars])
cmat <- stats::cor(X, use = "pairwise.complete.obs")
cor_out <- cbind(variable = rownames(cmat), as.data.frame(round(cmat, 3)))
write_csv_safe(cor_out, P("output", "tables", "S08_detect_correlations.csv"))
hp <- high_cor_pairs(cmat, 0.6)
write_csv_safe(hp, P("output", "tables", "S08_high_correlation_pairs.csv"))

vif <- vif_by_k(dom, size_vars, base = c("elev_m", "slope_deg", "gm4"),
                extra = c("d_road_any", "d_settle", "frac_crop", "frac_built", "frac_tree"))
write_csv_safe(vif, P("output", "tables", "S08_vif_by_k.csv"))

by_gm <- dom |>
  dplyr::group_by(gm4) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   med_d_road_major_km = round(stats::median(d_road_major) / 1000, 2),
                   med_d_road_any_km   = round(stats::median(d_road_any) / 1000, 2),
                   med_d_settle_km     = round(stats::median(d_settle) / 1000, 2),
                   med_frac_crop  = round(stats::median(frac_crop, na.rm = TRUE), 3),
                   med_frac_built = round(stats::median(frac_built, na.rm = TRUE), 3),
                   med_frac_tree  = round(stats::median(frac_tree, na.rm = TRUE), 3), .groups = "drop")
write_csv_safe(by_gm, P("output", "tables", "S08_detect_by_geomorph.csv"))

saveRDS(grid, P("data", "processed", "quadrature_det.rds"))
write_csv_safe(sf::st_drop_geometry(grid), P("data", "processed", "quadrature_det.csv"))

src <- data.frame(
  item = c("OSM Overpass endpoint", "OSM roads query", "OSM places query", "OSM download date",
           paste0("WorldCover ", wc_ids)),
  value = c(OSM_ENDPOINT, osm_query("roads", 0, 0, 0, 0), osm_query("places", 0, 0, 0, 0), format(Sys.time(), "%Y-%m-%d %H:%M"),
            wc_url(wc_ids)),
  note = c("OSM data (c) OpenStreetMap contributors, ODbL", "bbox replaced per tile", "bbox replaced per tile",
           "live database: re-download gives a different snapshot",
           rep("ESA WorldCover 10 m 2021 v200, CC BY 4.0 (confirm on ESA page)", length(wc_ids))))
write_csv_safe(src, P("output", "logs", "S08_external_sources.csv"))
write_csv_safe(file_checksums(c(jobs$file, wc_paths)), P("output", "logs", "S08_input_checksums.csv"))
write_csv_safe(data.frame(item = c("OSM_TILE_DEG_ROADS", "OSM_TILE_DEG_PLACES", "OSM_BUFFER_DEG", "HWY_MAJOR", "HWY_MINOR",
                                   "PLACE_SETTLE", "PLACE_TOWN", "WC_AGG", "WC_FOCAL_W", "WC_GROUPS"),
                          value = c(OSM_TILE_DEG_ROADS, OSM_TILE_DEG_PLACES, OSM_BUFFER_DEG,
                                    paste(HWY_MAJOR, collapse = ","), paste(HWY_MINOR, collapse = ","),
                                    paste(PLACE_SETTLE, collapse = ","), paste(PLACE_TOWN, collapse = ","),
                                    paste(WC_AGG, collapse = "x"), WC_FOCAL_W,
                                    paste(vapply(names(WC_GROUPS), function(g) paste0(g, "=", paste(WC_GROUPS[[g]], collapse = "+")),
                                                 character(1)), collapse = ";"))),
               P("output", "logs", "S08_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S08_sessionInfo.txt"))

cat("\n========== S08 SUMMARY ==========\n")
cat("OSM requests:", nrow(jobs), "| downloaded now:", sum(jobs$status == "downloaded"),
    "| cached:", sum(jobs$status == "cached"), "\n")
cat("\nOSM feature counts:\n"); print(as.data.frame(osm_counts), row.names = FALSE)
cat("\nDetectability covariates in analysis domain (distances in km):\n"); print(summ, row.names = FALSE)
cat("\nPairs with |r| >= 0.6 (log distances; fractions raw):\n"); print(hp, row.names = FALSE)
cat("\nVIF of each size predictor given terrain + geomorph (+ planned detectability set):\n"); print(vif, row.names = FALSE)
cat("\nBy geomorphology class:\n"); print(as.data.frame(by_gm), row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S08. Paste this console output back.\n")