# R/site_covariates.R -- functions for script 10 (covariates at arbitrary points; same code path as the grid)
# Needs: config.R, quadrature.R (nearest_dist, SIZE_K), covariates.R (join_geomorph, geomorph4),
#        detectability.R (HWY_MAJOR, classify_highway, PLACE_TOWN), WC_GROUPS

# Covariates at points `pts` (sf, CRS_ANALYSIS). ctx holds every surface; nothing else is read.
covariates_at <- function(pts, ctx) {
  n <- nrow(pts)
  xy <- sf::st_coordinates(pts)
  out <- data.frame(row = seq_len(n))
  ex <- terra::extract(ctx$terrain, xy)
  out$elev_m <- ex[["elev_m"]]; out$slope_deg <- ex[["slope_deg"]]
  g <- sf::st_sf(q_id = seq_len(n), geometry = sf::st_geometry(pts))
  out$geom_class <- join_geomorph(g, ctx$geomorph)
  out$gm4 <- geomorph4(out$geom_class)
  for (k in SIZE_K) out[[paste0("d_k", k)]] <- nearest_dist(pts, ctx$sets$with_ext[[as.character(k)]])
  out$d_lineament   <- nearest_dist(pts, ctx$lineament)
  out$d_road_major  <- nearest_dist(pts, ctx$roads_major)
  out$d_road_any    <- nearest_dist(pts, ctx$roads_any)
  out$d_settle      <- nearest_dist(pts, ctx$places_all)
  out$d_town        <- nearest_dist(pts, ctx$places_town)
  ll <- sf::st_coordinates(sf::st_transform(sf::st_geometry(pts), 4326))
  fx <- terra::extract(ctx$landcover, ll)
  for (v in names(ctx$landcover)) out[[v]] <- fx[[v]]
  out$row <- NULL
  out
}

# compare covariates recomputed at grid cell centres with the values stored by S06-S08
# rasters are stored as 32-bit floats, so values re-read from file differ from the in-memory values used when
# the grid was built by up to ~1e-5; tolerances below are set above that and far below any real error
qa_against_grid <- function(cov_new, grid_rows, tol_m = 1e-6, tol_frac = 1e-6, tol_terrain = 1e-3) {
  num <- c("elev_m", "slope_deg", paste0("d_k", SIZE_K), "d_lineament", "d_road_major", "d_road_any",
           "d_settle", "d_town", paste0("frac_", names(WC_GROUPS)))
  res <- do.call(rbind, lapply(num, function(v) {
    d <- abs(cov_new[[v]] - grid_rows[[v]])
    tol <- if (startsWith(v, "frac_")) tol_frac else if (v %in% c("elev_m", "slope_deg")) tol_terrain else tol_m
    data.frame(variable = v, n_compared = sum(!is.na(d)), max_abs_diff = signif(max(d, na.rm = TRUE), 3),
               n_mismatch = sum(d > tol | xor(is.na(cov_new[[v]]), is.na(grid_rows[[v]])), na.rm = TRUE))
  }))
  gm <- data.frame(variable = "gm4", n_compared = nrow(cov_new), max_abs_diff = NA_real_,
                   n_mismatch = sum(cov_new$gm4 != grid_rows$gm4 | xor(is.na(cov_new$gm4), is.na(grid_rows$gm4)), na.rm = TRUE))
  rbind(res, gm)
}

# cluster structure of a set of points for several block sizes (counts only; no covariates)
cluster_structure <- function(xy, blocks_km, domain_blocks = NULL) {
  do.call(rbind, lapply(blocks_km, function(b) {
    id <- cluster_ids(xy, b * 1000)
    tb <- table(id)
    data.frame(block_km = b, n_sites = nrow(xy), n_clusters_occupied = length(tb),
               n_blocks_in_domain = if (is.null(domain_blocks)) NA_integer_ else domain_blocks[[as.character(b)]],
               max_cluster_share = round(max(tb) / sum(tb), 3), median_sites_per_cluster = stats::median(as.numeric(tb)),
               singleton_clusters = sum(tb == 1L))
  }))
}