SIZE_K       <- c(1L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)   # k = 1 means "any drainage"; else WRIS ordsh >= k
CELL_M       <- 500
SURVEY_BUF_M <- 10000
SURVEY_MIN_N <- 5L
EXT_BUFFER_M <- 30000

survey_domain <- function(sites_sf, study, min_n = SURVEY_MIN_N, buf_m = SURVEY_BUF_M) {
  by_ref <- split(seq_len(nrow(sites_sf)), sites_sf$ref_id)
  parts <- lapply(by_ref, function(ix) {
    u <- sf::st_union(sf::st_geometry(sites_sf)[ix])
    h <- if (length(ix) >= min_n) sf::st_convex_hull(u) else u
    sf::st_buffer(h, buf_m)
  })
  dom <- sf::st_union(do.call(c, parts))
  sf::st_intersection(dom, study)
}

make_grid <- function(study, cell_m = CELL_M) {
  bb <- sf::st_bbox(study)
  xs <- seq(bb["xmin"] + cell_m / 2, bb["xmax"], by = cell_m)
  ys <- seq(bb["ymin"] + cell_m / 2, bb["ymax"], by = cell_m)
  g  <- expand.grid(x = xs, y = ys)
  p  <- sf::st_as_sf(g, coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
  p[lengths(sf::st_intersects(p, study)) > 0, ]
}

nearest_dist <- function(pts, lines_geom) {
  if (length(lines_geom) == 0) return(rep(NA_real_, nrow(pts)))
  i <- sf::st_nearest_feature(pts, lines_geom)
  as.numeric(sf::st_distance(pts, lines_geom[i], by_element = TRUE))
}

hr_outside <- function(hr, study, buf_m = EXT_BUFFER_M) {
  near <- lengths(sf::st_is_within_distance(hr, study, dist = buf_m)) > 0
  h <- hr[near, c("HYRIV_ID", "UPLAND_SKM")]
  h <- suppressWarnings(sf::st_difference(h, study))
  if (!nrow(h)) return(h)
  h <- h[!sf::st_is_empty(h), ]
  h[sf::st_geometry_type(h) %in% c("LINESTRING", "MULTILINESTRING"), ]
}

size_line_sets <- function(rv, hr_out, upland_by_k, ks = SIZE_K) {
  sets <- list(); sets_wris <- list()
  for (k in ks) {
    w <- if (k == 1L) rv else rv[!is.na(rv$ordsh) & rv$ordsh >= k, ]
    h <- if (k < 5L) hr_out[0, ] else hr_out[hr_out$UPLAND_SKM >= upland_by_k[as.character(k)], ]
    sets_wris[[as.character(k)]] <- sf::st_geometry(w)
    sets[[as.character(k)]] <- c(sf::st_geometry(w), sf::st_geometry(h))
  }
  list(with_ext = sets, wris_only = sets_wris)
}

log_d <- function(d_m) log(d_m / 1000 + 0.1)