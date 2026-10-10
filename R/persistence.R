SPACING_M   <- 20
CHUNK_M     <- 1000
NEIGH_RADII <- c(20, 40)
B03_CLOUD   <- 5000
MIN_ORD_PI  <- 5L

reach_points <- function(rv, spacing_m = SPACING_M, chunk_m = CHUNK_M) {
  len_m <- rv$length_km * 1000
  n <- pmax(1L, as.integer(round(len_m / spacing_m)))
  mp <- sf::st_line_sample(sf::st_geometry(rv), n = n, type = "regular")
  cc <- sf::st_coordinates(mp)
  feat <- as.integer(cc[, "L1"])
  cnt  <- tabulate(feat, nbins = nrow(rv))[feat]
  k    <- stats::ave(seq_along(feat), feat, FUN = seq_along)
  dist_along <- ((k - 0.5) / cnt) * len_m[feat]
  data.frame(seg_id = rv$seg_id[feat], k = k,
             x = unname(cc[, "X"]), y = unname(cc[, "Y"]),
             dist_along_m = dist_along, chunk = floor(dist_along / chunk_m),
             rivname = rv$rivname[feat], ordsh = rv$ordsh[feat], rilcode = rv$rilcode[feat],
             stringsAsFactors = FALSE)
}

neighbourhood_offsets <- function(radii = NEIGH_RADII, n_dir = 8L) {
  ang <- seq(0, 2 * pi, length.out = n_dir + 1L)[-(n_dir + 1L)]
  off <- do.call(rbind, lapply(radii, function(r) cbind(r * cos(ang), r * sin(ang))))
  rbind(c(0, 0), off)
}

extract_vals <- function(r, xy) {
  v <- terra::extract(r, xy)
  as.numeric(v[[ncol(v)]])
}

neighbourhood_mndwi <- function(b03, b11, xy, offsets = neighbourhood_offsets(), chunk_n = 100000L) {
  n_all <- nrow(xy); m <- nrow(offsets)
  mn  <- matrix(NA_real_, n_all, m)
  b3c <- rep(NA_real_, n_all)
  starts <- seq(1L, n_all, by = chunk_n)
  for (s in starts) {
    idx <- s:min(n_all, s + chunk_n - 1L)
    n <- length(idx)
    xa <- rep(xy[idx, 1], times = m) + rep(offsets[, 1], each = n)
    ya <- rep(xy[idx, 2], times = m) + rep(offsets[, 2], each = n)
    xyall <- cbind(xa, ya)
    g <- extract_vals(b03, xyall)
    h <- extract_vals(b11, xyall)
    den <- g + h
    v <- ifelse(is.na(den) | den <= 0, NA_real_, (g - h) / den)
    mn[idx, ] <- matrix(v, nrow = n, ncol = m)
    b3c[idx]  <- g[seq_len(n)]
    message("  extracted ", min(n_all, max(idx)), " / ", n_all, " points")
  }
  list(mndwi = mn, b03_centre = b3c)
}

row_max <- function(m) {
  out <- suppressWarnings(apply(m, 1, max, na.rm = TRUE))
  out[is.infinite(out)] <- NA_real_
  out
}

add_metrics <- function(pts, nb) {
  pts$mndwi_c     <- nb$mndwi[, 1]
  pts$mndwi_max20 <- row_max(nb$mndwi[, 1:9])
  pts$mndwi_max40 <- row_max(nb$mndwi)
  pts$b03_c       <- nb$b03_centre
  pts$cloudy      <- !is.na(pts$b03_c) & pts$b03_c > B03_CLOUD
  pts
}

baseline_points <- function(bnd_union, core, large_lines, n = 40000L, min_dist_m = 1000, keep = 20000L) {
  p <- sf::st_sf(geometry = sf::st_sample(bnd_union, n, type = "random"))
  in_core <- lengths(sf::st_intersects(p, core)) > 0
  near    <- lengths(sf::st_is_within_distance(p, large_lines, dist = min_dist_m)) > 0
  p <- p[!in_core & !near, ]
  p <- p[seq_len(min(nrow(p), keep)), ]
  xy <- sf::st_coordinates(p)
  data.frame(x = unname(xy[, 1]), y = unname(xy[, 2]))
}

PI_P <- 0.5; PI_I <- 0.15

subsegment_table <- function(pts, thr = 0) {
  pts$sub_id <- paste(pts$seg_id, pts$chunk, sep = "_")
  d <- pts[pts$valid, ]
  tot <- dplyr::summarise(dplyr::group_by(pts, sub_id, seg_id, chunk, rivname, ordsh, rilcode),
                          n_total = dplyr::n(), n_removed = sum(in_removed), .groups = "drop")
  s <- dplyr::summarise(dplyr::group_by(d, sub_id), n_valid = dplyr::n(),
                        med_mndwi40 = stats::median(mndwi_max40),
                        pi = mean(mndwi_max40 > thr),
                        pi_centre = mean(mndwi_c > thr, na.rm = TRUE),
                        x_mid = mean(x), y_mid = mean(y), .groups = "drop")
  out <- dplyr::left_join(tot, s, by = "sub_id")
  out$usable <- !is.na(out$n_valid) & out$n_valid >= 10L
  out$pi_class <- dplyr::case_when(!out$usable ~ NA_character_,
                                   out$pi >= PI_P ~ "P_candidate",
                                   out$pi >= PI_I ~ "I",
                                   TRUE ~ "S")
  out$length_km_valid <- ifelse(is.na(out$n_valid), 0, out$n_valid * SPACING_M / 1000)
  out
}

auc_mw <- function(pos, neg) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

order_group <- function(o) {
  dplyr::case_when(o >= 9 ~ "order9plus", o >= 6 ~ "order6to8", o >= 5 ~ "order5", TRUE ~ "below5")
}

fp_table <- function(pts, base, thresholds = c(0, -0.05, -0.1, -0.2)) {
  g <- list(baseline_land = base$mndwi_max40[base$valid],
            order5 = pts$mndwi_max40[pts$valid & pts$og == "order5"],
            order6to8 = pts$mndwi_max40[pts$valid & pts$og == "order6to8"],
            order9plus = pts$mndwi_max40[pts$valid & pts$og == "order9plus"])
  gc <- list(baseline_land = base$mndwi_c[base$valid],
             order5 = pts$mndwi_c[pts$valid & pts$og == "order5"],
             order6to8 = pts$mndwi_c[pts$valid & pts$og == "order6to8"],
             order9plus = pts$mndwi_c[pts$valid & pts$og == "order9plus"])
  out <- expand.grid(group = names(g), threshold = thresholds, stringsAsFactors = FALSE)
  out$n <- vapply(out$group, function(k) length(g[[k]]), numeric(1))
  out$share_max40_gt <- mapply(function(k, t) round(mean(g[[k]] > t, na.rm = TRUE), 4), out$group, out$threshold)
  out$share_centre_gt <- mapply(function(k, t) round(mean(gc[[k]] > t, na.rm = TRUE), 4), out$group, out$threshold)
  out
}