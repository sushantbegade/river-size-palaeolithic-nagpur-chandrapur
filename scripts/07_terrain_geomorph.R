suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf); library(terra)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "covariates.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S07 terrain and geomorphology ==")

f_h <- P("data", "processed", "hydro_prep.gpkg")
f_q <- P("data", "processed", "quadrature.rds")
f_d <- DATA("DEM", "DEM.tif")
for (f in c(f_h, f_q, f_d)) if (!file.exists(f)) stop("Missing: ", f, "\n  Run S02 and S06 first.", call. = FALSE)

bnd   <- sf::st_read(f_h, layer = "boundary", quiet = TRUE)
study <- sf::st_union(sf::st_geometry(bnd))
gm    <- sf::st_read(f_h, layer = "geomorph", quiet = TRUE)
grid  <- readRDS(f_q)
grid$q_id <- seq_len(nrow(grid))

message("preparing DEM (crop, zero to NA, project to UTM 30 m, fill small holes)... a few minutes")
dem_utm <- prep_dem_utm(f_d, study)
message("terrain (slope) and aggregation to ~510 m...")
tr <- terrain_stack(dem_utm)
terra::writeRaster(tr, P("data", "processed", "terrain_510m.tif"), overwrite = TRUE)

xy  <- cbind(grid$x, grid$y)
ex  <- terra::extract(tr, xy)
grid$elev_m    <- ex[["elev_m"]]
grid$slope_deg <- ex[["slope_deg"]]
grid$geom_class <- join_geomorph(grid, gm)
grid$gm4 <- geomorph4(grid$geom_class)

dd  <- sf::st_drop_geometry(grid)
dom <- dd[dd$in_survey & !dd$in_core, ]
message("domain points: ", nrow(dom), " | NA elevation: ", sum(is.na(dom$elev_m)),
        " | NA slope: ", sum(is.na(dom$slope_deg)), " | NA geomorph class: ", sum(is.na(dom$gm4)))

terr <- data.frame(variable = c("elev_m", "slope_deg"),
                   p01 = c(stats::quantile(dom$elev_m, .01, na.rm = TRUE), stats::quantile(dom$slope_deg, .01, na.rm = TRUE)),
                   p25 = c(stats::quantile(dom$elev_m, .25, na.rm = TRUE), stats::quantile(dom$slope_deg, .25, na.rm = TRUE)),
                   median = c(stats::median(dom$elev_m, na.rm = TRUE), stats::median(dom$slope_deg, na.rm = TRUE)),
                   p75 = c(stats::quantile(dom$elev_m, .75, na.rm = TRUE), stats::quantile(dom$slope_deg, .75, na.rm = TRUE)),
                   p99 = c(stats::quantile(dom$elev_m, .99, na.rm = TRUE), stats::quantile(dom$slope_deg, .99, na.rm = TRUE)),
                   n_na = c(sum(is.na(dom$elev_m)), sum(is.na(dom$slope_deg))))
terr[, 2:6] <- round(terr[, 2:6], 2)
write_csv_safe(terr, P("output", "tables", "S07_terrain_summary.csv"))

gshare <- dom |>
  dplyr::count(gm4, geom_class, name = "n_cells") |>
  dplyr::mutate(pct = round(100 * n_cells / sum(n_cells), 2))
write_csv_safe(gshare, P("output", "tables", "S07_geomorph_share_domain.csv"))

vars <- c(paste0("d_k", SIZE_K), "d_lineament")
ld <- as.data.frame(lapply(dom[, vars], log_d))
cmat <- stats::cor(cbind(elev = dom$elev_m, slope = dom$slope_deg, ld), use = "pairwise.complete.obs")
cor_cov <- data.frame(variable = c("elev", "slope"), round(cmat[c("elev", "slope"), vars], 3))
write_csv_safe(cor_cov, P("output", "tables", "S07_terrain_vs_distance_correlations.csv"))

by_gm <- dom |>
  dplyr::group_by(gm4) |>
  dplyr::summarise(n_cells = dplyr::n(), median_elev = round(stats::median(elev_m, na.rm = TRUE)),
                   median_slope = round(stats::median(slope_deg, na.rm = TRUE), 2),
                   median_d_k5_km = round(stats::median(d_k5) / 1000, 2),
                   median_d_k6_km = round(stats::median(d_k6) / 1000, 2),
                   median_d_k7_km = round(stats::median(d_k7) / 1000, 2), .groups = "drop")
write_csv_safe(by_gm, P("output", "tables", "S07_covariates_by_geomorph.csv"))

saveRDS(grid, P("data", "processed", "quadrature_cov.rds"))
write_csv_safe(sf::st_drop_geometry(grid), P("data", "processed", "quadrature_cov.csv"))
write_csv_safe(data.frame(item = c("AGG_FACTOR", "DEM_MARGIN_DEG", "dem_resolution_m", "hole_fill"),
                          value = c(AGG_FACTOR, DEM_MARGIN_DEG, 30, "focal mean w=5, NA cells only")),
               P("output", "logs", "S07_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S07_sessionInfo.txt"))

cat("\n========== S07 SUMMARY ==========\n")
cat("Domain grid points:", nrow(dom), "\n")
cat("\nTerrain in domain:\n"); print(terr, row.names = FALSE)
cat("\nGeomorphology share of domain:\n"); print(as.data.frame(gshare), row.names = FALSE)
cat("\nCovariates by geomorphology class:\n"); print(as.data.frame(by_gm), row.names = FALSE)
cat("\nCorrelation of terrain with log distance predictors:\n"); print(cor_cov, row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S07. Paste this console output back.\n")