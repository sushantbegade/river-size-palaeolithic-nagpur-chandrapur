suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "quadrature.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S06 quadrature and predictors ==")

f_sites <- P("data", "processed", "sites_clean.rds")
f_h     <- P("data", "processed", "hydro_prep.gpkg")
f_e     <- P("data", "processed", "external_region.gpkg")
f_cal   <- P("output", "tables", "S05_wris_order_vs_hr_upland.csv")
for (f in c(f_sites, f_h, f_e, f_cal)) if (!file.exists(f)) stop("Missing: ", f, "\n  Run S01-S05 first.", call. = FALSE)

sites <- readRDS(f_sites)
sites_sf <- sf::st_as_sf(sites[, c("site_id", "ref_id", "x", "y")], coords = c("x", "y"),
                         crs = CRS_ANALYSIS, remove = FALSE)
bnd   <- sf::st_read(f_h, layer = "boundary", quiet = TRUE)
study <- sf::st_union(sf::st_geometry(bnd))
rv    <- sf::st_read(f_h, layer = "rivers", quiet = TRUE)
core  <- sf::st_union(sf::st_geometry(sf::st_read(f_h, layer = "removed_core", quiet = TRUE)))
chan  <- if ("removed_channel" %in% sf::st_layers(f_h)$name)
  sf::st_union(sf::st_geometry(sf::st_read(f_h, layer = "removed_channel", quiet = TRUE))) else
    sf::st_sfc(sf::st_polygon(), crs = CRS_ANALYSIS)
hr    <- sf::st_read(f_e, layer = "hydrorivers", quiet = TRUE)

cal <- readr::read_csv(f_cal, show_col_types = FALSE)
upl <- stats::setNames(cal$upland_median, as.character(cal$ordsh))
need <- as.character(SIZE_K[SIZE_K >= 5L])
if (!all(need %in% names(upl))) stop("Calibration table lacks WRIS orders: ", paste(setdiff(need, names(upl)), collapse = ","), call. = FALSE)
message("HydroRIVERS catchment thresholds used for the outside extension (km2): ",
        paste(paste0("k", need, "=", upl[need]), collapse = ", "))

dom_survey <- survey_domain(sites_sf, study)
grid <- make_grid(study)
grid$in_survey  <- lengths(sf::st_intersects(grid, dom_survey)) > 0
dom_survey5 <- survey_domain(sites_sf, study, buf_m = 5000)
grid$in_survey5 <- lengths(sf::st_intersects(grid, dom_survey5)) > 0
grid$in_core    <- lengths(sf::st_intersects(grid, core)) > 0
grid$in_channel <- if (length(chan) && !all(sf::st_is_empty(chan))) lengths(sf::st_intersects(grid, chan)) > 0 else FALSE
grid$area_km2 <- (CELL_M / 1000)^2
message("grid points in study: ", nrow(grid), " | in survey domain: ", sum(grid$in_survey))

hro  <- hr_outside(hr, study)
message("HydroRIVERS reaches outside the study within ", EXT_BUFFER_M / 1000, " km: ", nrow(hro))
sets <- size_line_sets(rv, hro, upl)
for (k in SIZE_K) {
  kk <- as.character(k)
  grid[[paste0("d_k", k)]]       <- nearest_dist(grid, sets$with_ext[[kk]])
  grid[[paste0("d_k", k, "_wris")]] <- nearest_dist(grid, sets$wris_only[[kk]])
  message("  distances done for k = ", k)
}
lin_path <- DATA("Linement", "Linement.shp")
if (file.exists(lin_path)) {
  lin <- sf::st_geometry(read_vector_utm(lin_path, "Lineaments"))
  grid$d_lineament <- nearest_dist(grid, lin)
} else { grid$d_lineament <- NA_real_; message("lineament file not found: negative control skipped") }

dd  <- sf::st_drop_geometry(grid)
dom <- dd[dd$in_survey & !dd$in_core, ]
dvars <- c(paste0("d_k", SIZE_K), "d_lineament")
ld <- as.data.frame(lapply(dom[, dvars], log_d)); names(ld) <- dvars

cor_tab <- as.data.frame(round(stats::cor(ld, use = "pairwise.complete.obs"), 3))
cor_tab <- cbind(variable = rownames(cor_tab), cor_tab)
write_csv_safe(cor_tab, P("output", "tables", "S06_distance_correlations_logkm.csv"))

qs <- t(sapply(dom[, dvars], function(v) round(stats::quantile(v / 1000, c(.1, .25, .5, .75, .9), na.rm = TRUE), 2)))
qs <- cbind(variable = rownames(qs), as.data.frame(qs))
write_csv_safe(qs, P("output", "tables", "S06_distance_quantiles_km.csv"))

pairs <- t(utils::combn(paste0("d_k", SIZE_K), 2))
disc <- data.frame(a = pairs[, 1], b = pairs[, 2])
disc$share_differ_gt1km <- vapply(seq_len(nrow(disc)), function(i)
  round(mean(abs(dom[[disc$a[i]]] - dom[[disc$b[i]]]) > 1000, na.rm = TRUE), 3), numeric(1))
disc$rank_cor <- vapply(seq_len(nrow(disc)), function(i)
  round(stats::cor(dom[[disc$a[i]]], dom[[disc$b[i]]], method = "spearman", use = "pairwise.complete.obs"), 3), numeric(1))
write_csv_safe(disc, P("output", "tables", "S06_discriminating_area.csv"))

border <- data.frame(k = SIZE_K)
border$share_ext_closer_gt500m <- vapply(SIZE_K, function(k)
  round(mean(dom[[paste0("d_k", k, "_wris")]] - dom[[paste0("d_k", k)]] > 500, na.rm = TRUE), 3), numeric(1))
border$median_gain_m_where_closer <- vapply(SIZE_K, function(k) {
  g <- dom[[paste0("d_k", k, "_wris")]] - dom[[paste0("d_k", k)]]
  round(stats::median(g[g > 500], na.rm = TRUE))
}, numeric(1))
write_csv_safe(border, P("output", "tables", "S06_border_effect.csv"))

area_tab <- data.frame(
  item = c("study_km2", "survey_domain_10km_km2", "survey_domain_5km_km2", "grid_points_study",
           "grid_points_survey_10km", "grid_points_survey_5km", "grid_points_survey_minus_core",
           "grid_points_survey_minus_core_channel", "grid_cell_km2"),
  value = c(round(as.numeric(sf::st_area(study)) / 1e6, 1), round(as.numeric(sf::st_area(dom_survey)) / 1e6, 1),
            round(as.numeric(sf::st_area(dom_survey5)) / 1e6, 1), nrow(grid), sum(grid$in_survey),
            sum(grid$in_survey5), sum(grid$in_survey & !grid$in_core),
            sum(grid$in_survey & !grid$in_core & !grid$in_channel), CELL_M^2 / 1e6))
write_csv_safe(area_tab, P("output", "tables", "S06_domain_summary.csv"))

saveRDS(grid, P("data", "processed", "quadrature.rds"))
write_csv_safe(dd, P("data", "processed", "quadrature.csv"))
gp <- P("data", "processed", "domains.gpkg"); unlink(gp)
sf::st_write(sf::st_sf(domain = "survey", geometry = dom_survey), gp, layer = "survey_domain", quiet = TRUE)
write_csv_safe(data.frame(item = c("SIZE_K", "CELL_M", "SURVEY_BUF_M", "SURVEY_MIN_N", "EXT_BUFFER_M", "upland_thresholds_km2"),
                          value = c(paste(SIZE_K, collapse = ","), CELL_M, SURVEY_BUF_M, SURVEY_MIN_N, EXT_BUFFER_M,
                                    paste(paste0("k", need, "=", upl[need]), collapse = ";"))),
               P("output", "logs", "S06_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S06_sessionInfo.txt"))

cat("\n========== S06 SUMMARY ==========\n")
print(as.data.frame(area_tab), row.names = FALSE)
cat("\nDistance quantiles (km) inside analysis domain (survey domain minus removed core):\n")
print(as.data.frame(qs), row.names = FALSE)
cat("\nCorrelations among log distances (km + 0.1):\n"); print(cor_tab, row.names = FALSE)
cat("\nWhere predictors differ (share of domain with |difference| > 1 km), and rank correlation:\n")
print(disc, row.names = FALSE)
cat("\nBorder effect: share of domain where the HydroRIVERS extension makes the nearest reach >500 m closer:\n")
print(border, row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S06. Paste this console output back.\n")