suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(stringr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))

set.seed(SEED)
message("== S02 prepare hydrography ==")
t0 <- Sys.time()

f_bnd <- DATA("District Boundary", "Nagpur_and_Chandrapur_Boundry.shp")
f_riv <- DATA("Rivers", "Rivers.shp")
f_wb  <- DATA("Waterbody", "Waterbody.shp")
f_gm  <- DATA("Geomorphology", "Geomorphology.shp")

bnd  <- prep_boundary(read_vector_utm(f_bnd, "District boundary"))
bu   <- sf::st_union(sf::st_geometry(bnd))
study_km2 <- sum(bnd$area_km2)
message("study area (districts): ", round(study_km2, 1), " km2")

rv <- prep_rivers(read_vector_utm(f_riv, "Rivers"), bu)
message("river segments: ", nrow(rv), " | in study: ", sum(rv$in_study),
        " | in 25 km buffer: ", sum(rv$in_buf25))

wb <- prep_waterbodies(read_vector_utm(f_wb, "Waterbody"), bu)
message("waterbody polygons: ", nrow(wb), " | in study: ", sum(wb$in_study))

gm <- prep_geomorph(read_vector_utm(f_gm, "Geomorphology"))
message("geomorphology polygons: ", nrow(gm))

masks <- build_removed_masks(gm, wb, bu)
core_km2 <- sum(as.numeric(sf::st_area(masks$core))) / 1e6
chan_km2 <- if (nrow(masks$channel)) sum(as.numeric(sf::st_area(masks$channel))) / 1e6 else 0

gpkg <- P("data", "processed", "hydro_prep.gpkg")
dir.create(dirname(gpkg), recursive = TRUE, showWarnings = FALSE)
unlink(gpkg)
sf::st_write(bnd, gpkg, layer = "boundary", quiet = TRUE)
sf::st_write(rv,  gpkg, layer = "rivers", quiet = TRUE)
sf::st_write(wb,  gpkg, layer = "waterbodies", quiet = TRUE)
sf::st_write(gm,  gpkg, layer = "geomorph", quiet = TRUE)
sf::st_write(masks$core, gpkg, layer = "removed_core", quiet = TRUE)
if (nrow(masks$channel)) sf::st_write(masks$channel, gpkg, layer = "removed_channel", quiet = TRUE)
message("wrote data/processed/hydro_prep.gpkg")

write_csv_safe(network_by_order(rv),      P("output", "tables", "S02_network_by_order.csv"))
write_csv_safe(large_river_table(rv, 6L), P("output", "tables", "S02_large_rivers_order6plus.csv"))
write_csv_safe(waterbody_by_code(wb),     P("output", "tables", "S02_waterbody_by_code.csv"))
write_csv_safe(reservoir_table(wb),       P("output", "tables", "S02_reservoirs_code02.csv"))
write_csv_safe(geomorph_table(gm, study_km2), P("output", "tables", "S02_geomorph_classes.csv"))
removed <- data.frame(item = c("study_area_km2", "removed_core_km2", "removed_core_pct",
                               "river_channel_km2", "river_channel_pct"),
                      value = round(c(study_km2, core_km2, 100 * core_km2 / study_km2,
                                      chan_km2, 100 * chan_km2 / study_km2), 3))
write_csv_safe(removed, P("output", "tables", "S02_removed_area.csv"))

shp_parts <- unlist(lapply(c(f_bnd, f_riv, f_wb, f_gm), function(f)
  paste0(tools::file_path_sans_ext(f), c(".shp", ".dbf", ".shx", ".prj"))))
shp_parts <- shp_parts[file.exists(shp_parts)]
write_csv_safe(file_checksums(shp_parts), P("output", "logs", "S02_input_checksums.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S02_sessionInfo.txt"))

cat("\n========== S02 SUMMARY ==========\n")
cat("Study area km2:", round(study_km2, 1), "\n")
cat("\nRivers: CRS of output =", sf::st_crs(rv)$epsg, "| segments", nrow(rv),
    "| NA ordsh:", sum(is.na(rv$ordsh)), "\n")
cat("Network length in study (km): all", round(sum(rv$length_km[rv$in_study]), 1),
    "| large5", round(sum(rv$length_km[rv$in_study & rv$large5]), 1),
    "| large6", round(sum(rv$length_km[rv$in_study & rv$large6]), 1),
    "| large7", round(sum(rv$length_km[rv$in_study & rv$large7]), 1),
    "| trunk9", round(sum(rv$length_km[rv$in_study & rv$trunk9]), 1), "\n")
cat("\nLength (km) by stream order, inside study:\n")
print(as.data.frame(stats::aggregate(length_km ~ ordsh, rv[rv$in_study, ] |> sf::st_drop_geometry(), sum)),
      row.names = FALSE)
cat("\nLarge rivers (order >= 6) in study:\n")
print(as.data.frame(utils::head(large_river_table(rv, 6L), 30)), row.names = FALSE)
cat("\nWaterbodies by wbcode:\n"); print(as.data.frame(waterbody_by_code(wb)), row.names = FALSE)
cat("\nReservoirs (wbcode 02), largest 15:\n")
print(as.data.frame(utils::head(reservoir_table(wb), 15)), row.names = FALSE)
cat("\nWaterbody area check, layer area_ha vs computed (ratio median):",
    round(stats::median(wb$area_ha / (wb$area_km2 * 100), na.rm = TRUE), 3), "\n")
cat("\nGeomorphology classes (km2):\n")
print(as.data.frame(stats::aggregate(area_km2 ~ geom_class, sf::st_drop_geometry(gm), sum)), row.names = FALSE)
cat("\nGeomorphology descriptions mapped to 'other':\n")
print(unique(gm$descriptio[gm$geom_class == "other"]))
cat("\nRemoved areas: core", round(core_km2, 1), "km2 (", round(100 * core_km2 / study_km2, 2),
    "% ) | river channel", round(chan_km2, 1), "km2 (", round(100 * chan_km2 / study_km2, 2), "% )\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S02. Paste this console output back.\n")