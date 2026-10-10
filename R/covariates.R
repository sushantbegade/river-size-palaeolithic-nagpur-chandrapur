AGG_FACTOR <- 17L
DEM_MARGIN_DEG <- 0.1

prep_dem_utm <- function(dem_path, study) {
  dem <- terra::rast(dem_path)
  bb  <- sf::st_bbox(sf::st_transform(sf::st_sf(geometry = study), 4326))
  e   <- terra::ext(bb["xmin"] - DEM_MARGIN_DEG, bb["xmax"] + DEM_MARGIN_DEG,
                    bb["ymin"] - DEM_MARGIN_DEG, bb["ymax"] + DEM_MARGIN_DEG)
  dem <- terra::crop(dem, e)
  dem <- terra::ifel(dem == 0, NA, dem)
  dem <- terra::project(dem, paste0("EPSG:", CRS_ANALYSIS), res = 30, method = "bilinear")
  terra::focal(dem, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
}

terrain_stack <- function(dem_utm, fact = AGG_FACTOR) {
  slope <- terra::terrain(dem_utm, v = "slope", unit = "degrees")
  st <- c(dem_utm, slope)
  names(st) <- c("elev_m", "slope_deg")
  terra::aggregate(st, fact = fact, fun = "mean", na.rm = TRUE)
}

geomorph4 <- function(geom_class) {
  dplyr::case_when(
    geom_class == "pediplain"                                 ~ "pediplain",
    geom_class == "plateau_hills"                             ~ "plateau_hills",
    geom_class %in% c("alluvial_floodplain", "river_channel") ~ "alluvial",
    TRUE                                                      ~ NA_character_
  )
}

join_geomorph <- function(grid, gm) {
  g <- grid[, "q_id"]
  j <- suppressWarnings(sf::st_join(g, gm[, "geom_class"], join = sf::st_intersects, left = TRUE))
  j <- j[!duplicated(j$q_id), ]
  j$geom_class[match(grid$q_id, j$q_id)]
}