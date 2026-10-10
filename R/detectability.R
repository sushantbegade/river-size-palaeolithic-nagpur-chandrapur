# R/detectability.R -- functions for script 08 (no data are read at source time)
# Needs: quadrature.R (nearest_dist, log_d, SIZE_K), config.R (CRS_ANALYSIS)

# ---- OpenStreetMap via Overpass API ----------------------------------------------
OSM_ENDPOINT        <- "https://overpass-api.de/api/interpreter"
OSM_TILE_DEG_ROADS  <- 0.5     # lower to 0.25 if road tiles time out
OSM_TILE_DEG_PLACES <- 1.0
OSM_BUFFER_DEG      <- 0.3     # about 33 km beyond the districts (border effect)
OSM_PAUSE_S         <- 10      # pause after each successful download
OSM_RETRY_WAIT_S    <- 30      # Overpass policy: wait 30 s after 429 / 504
OSM_MAX_TRY         <- 4L
OSM_USER_AGENT      <- "palaeolithic-waterbody-study/1.0 (academic research; R download.file)"

HWY_MAJOR <- c("motorway", "trunk", "primary", "secondary", "tertiary",
               "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link")
HWY_MINOR <- c("unclassified", "residential")          # tracks excluded on purpose
PLACE_SETTLE <- c("city", "town", "village", "hamlet")
PLACE_TOWN   <- c("city", "town")

# ---- ESA WorldCover 2021 v200 -------------------------------------------------------
WC_BASE       <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map"
WC_MARGIN_DEG <- 0.02
WC_AGG        <- c(11L, 11L)   # 10 m pixels -> about 100 m cells (cols, rows)
WC_FOCAL_W    <- 5L            # 5 x 5 of those cells -> about 0.5 km window
WC_GROUPS <- list(tree = 10, shrubgrass = c(20, 30), crop = 40,
                  built = 50, bare = 60, wet = c(80, 90))

study_bbox_ll <- function(study, buf_deg = 0) {
  bb <- sf::st_bbox(sf::st_transform(sf::st_sf(geometry = study), 4326))
  c(xmin = unname(bb["xmin"]) - buf_deg, ymin = unname(bb["ymin"]) - buf_deg,
    xmax = unname(bb["xmax"]) + buf_deg, ymax = unname(bb["ymax"]) + buf_deg)
}

osm_tiles <- function(bb, tile_deg) {
  ix <- seq.int(floor(bb[["xmin"]] / tile_deg), ceiling(bb[["xmax"]] / tile_deg) - 1L)
  iy <- seq.int(floor(bb[["ymin"]] / tile_deg), ceiling(bb[["ymax"]] / tile_deg) - 1L)
  g <- expand.grid(ix = ix, iy = iy)
  w <- g$ix * tile_deg
  s <- g$iy * tile_deg
  data.frame(tile = sprintf("E%05d_N%05d", as.integer(round(w * 100)), as.integer(round(s * 100))),
             w = w, s = s, e = w + tile_deg, n = s + tile_deg, stringsAsFactors = FALSE)
}

osm_query <- function(kind, s, w, n, e) {
  box <- sprintf("(%.5f,%.5f,%.5f,%.5f)", s, w, n, e)
  if (kind == "roads") {
    rx <- paste(c(HWY_MAJOR, HWY_MINOR), collapse = "|")
    return(paste0('[out:xml][timeout:180];(way["highway"~"^(', rx, ')$"]', box, ';);(._;>;);out body qt;'))
  }
  if (kind == "places") {
    rx <- paste(PLACE_SETTLE, collapse = "|")
    return(paste0('[out:xml][timeout:120];(node["place"~"^(', rx, ')$"]', box, ';);out body;'))
  }
  stop("Unknown OSM kind: ", kind, call. = FALSE)
}

.raw_txt <- function(r) rawToChar(r[r != as.raw(0)])

# TRUE only for a complete Overpass XML answer (no <remark> = no timeout / memory error)
osm_xml_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  sz <- file.size(path)
  if (is.na(sz) || sz < 50) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  head <- .raw_txt(readBin(con, "raw", min(sz, 2000)))
  seek(con, where = max(0, sz - 1000), origin = "start")
  tail <- .raw_txt(readBin(con, "raw", 1000))
  grepl("<osm", head, fixed = TRUE, useBytes = TRUE) &&
    grepl("</osm>", tail, fixed = TRUE, useBytes = TRUE) &&
    !grepl("<remark", tail, fixed = TRUE, useBytes = TRUE)
}

overpass_fetch <- function(query, dest, tries = OSM_MAX_TRY, wait = OSM_RETRY_WAIT_S) {
  if (osm_xml_ok(dest)) return("cached")
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  url <- paste0(OSM_ENDPOINT, "?data=", utils::URLencode(query, reserved = TRUE))
  old <- options(timeout = 900, HTTPUserAgent = OSM_USER_AGENT)
  on.exit(options(old), add = TRUE)
  tmp <- paste0(dest, ".part")
  for (i in seq_len(tries)) {
    res <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl"),
                    error = function(e) e, warning = function(w) w)
    if (!inherits(res, "condition") && osm_xml_ok(tmp)) {
      file.rename(tmp, dest)
      return("downloaded")
    }
    if (file.exists(tmp)) unlink(tmp)
    if (i < tries) {
      message("    attempt ", i, " failed; waiting ", wait, " s")
      Sys.sleep(wait)
    }
  }
  "failed"
}

read_osm_layer <- function(path, layer) {
  x <- tryCatch(suppressWarnings(sf::st_read(path, layer = layer, quiet = TRUE)),
                error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x
}

read_osm_tiles <- function(files, layer, keep) {
  parts <- lapply(files, function(f) {
    x <- read_osm_layer(f, layer)
    if (is.null(x)) return(NULL)
    for (k in setdiff(keep, names(x))) x[[k]] <- NA_character_
    x[, keep]
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  x <- do.call(rbind, parts)
  x[!duplicated(x$osm_id), ]          # ways crossing tile borders appear in several tiles
}

classify_highway <- function(h) {
  dplyr::case_when(h %in% HWY_MAJOR ~ "major", h %in% HWY_MINOR ~ "minor", TRUE ~ NA_character_)
}

osm_count_table <- function(rd, pl) {
  tr <- table(rd$highway)
  r <- data.frame(layer = "roads", class = names(tr), n = as.integer(tr), stringsAsFactors = FALSE)
  rl <- tapply(as.numeric(sf::st_length(rd)) / 1000, rd$highway, sum)
  r$length_km <- round(as.numeric(rl[r$class]), 1)
  tp <- table(pl$place)
  p <- data.frame(layer = "places", class = names(tp), n = as.integer(tp), length_km = NA_real_,
                  stringsAsFactors = FALSE)
  rbind(r, p)
}

# ---- WorldCover -------------------------------------------------------------------
wc_tile_ids <- function(bb) {
  ix <- seq.int(floor(bb[["xmin"]] / 3), floor((bb[["xmax"]] - 1e-9) / 3))
  iy <- seq.int(floor(bb[["ymin"]] / 3), floor((bb[["ymax"]] - 1e-9) / 3))
  g <- expand.grid(ix = ix, iy = iy)
  lon <- g$ix * 3L; lat <- g$iy * 3L
  sprintf("%s%02d%s%03d", ifelse(lat >= 0, "N", "S"), abs(lat), ifelse(lon >= 0, "E", "W"), abs(lon))
}

wc_file <- function(id) sprintf("ESA_WorldCover_10m_2021_v200_%s_Map.tif", id)
wc_url  <- function(id) paste0(WC_BASE, "/", wc_file(id))

is_valid_tif <- function(path) {
  if (!file.exists(path) || is.na(file.size(path)) || file.size(path) < 1e5) return(FALSE)
  ok <- tryCatch({ r <- terra::rast(path); terra::ncol(r) > 0 }, error = function(e) FALSE)
  isTRUE(ok)
}

download_tif_if_missing <- function(url, dest) {
  if (is_valid_tif(dest)) { message("already present: ", basename(dest)); return(invisible(dest)) }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(7200, getOption("timeout"))); on.exit(options(old), add = TRUE)
  message("downloading ", basename(dest), " (about 100-200 MB) ...")
  res <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl"),
                  error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition") || !is_valid_tif(dest)) {
    if (file.exists(dest)) unlink(dest)
    stop("Download failed: ", url,
         "\n  MANUAL OPTION: open the link in a browser, save the file exactly as",
         "\n  ", dest, "\n  then run this script again.", call. = FALSE)
  }
  invisible(dest)
}

# class-group fractions: 10 m pixels -> 100 m cells (mean of 0/1) -> 5 x 5 focal mean
wc_group_fractions <- function(paths, bb, margin = WC_MARGIN_DEG, agg = WC_AGG, w = WC_FOCAL_W) {
  r <- if (length(paths) == 1L) terra::rast(paths) else terra::vrt(paths)
  e <- terra::ext(bb[["xmin"]] - margin, bb[["xmax"]] + margin, bb[["ymin"]] - margin, bb[["ymax"]] + margin)
  r <- terra::crop(r, e)
  out <- list()
  for (g in names(WC_GROUPS)) {
    message("  land-cover group: ", g)
    cls <- WC_GROUPS[[g]]
    b <- terra::classify(r, cbind(cls - 0.5, cls + 0.5, 1), others = 0)
    a <- terra::aggregate(b, fact = agg, fun = "mean", na.rm = TRUE)
    out[[g]] <- terra::focal(a, w = w, fun = "mean", na.rm = TRUE)
  }
  fr <- terra::rast(out)
  names(fr) <- paste0("frac_", names(WC_GROUPS))
  fr
}

# ---- summaries ----------------------------------------------------------------------
cov_summary <- function(dom, vars, km_vars) {
  do.call(rbind, lapply(vars, function(v) {
    x <- dom[[v]]
    if (v %in% km_vars) x <- x / 1000
    q <- stats::quantile(x, c(.01, .25, .5, .75, .99), na.rm = TRUE, names = FALSE)
    data.frame(variable = v, unit = ifelse(v %in% km_vars, "km", "fraction"),
               p01 = round(q[1], 3), p25 = round(q[2], 3), median = round(q[3], 3),
               p75 = round(q[4], 3), p99 = round(q[5], 3), n_na = sum(is.na(x)))
  }))
}

high_cor_pairs <- function(cmat, thr = 0.6) {
  ij <- which(upper.tri(cmat) & abs(cmat) >= thr, arr.ind = TRUE)
  if (!nrow(ij)) return(data.frame(a = character(0), b = character(0), r = numeric(0)))
  out <- data.frame(a = rownames(cmat)[ij[, 1]], b = colnames(cmat)[ij[, 2]], r = round(cmat[ij], 3))
  out[order(-abs(out$r)), ]
}

# R2 / VIF of each size predictor (log distance) regressed on the other covariates
vif_by_k <- function(dom, size_vars, base, extra) {
  d <- dom
  d$gm4 <- factor(d$gm4)
  for (v in c(size_vars, extra[grepl("^d_", extra)])) d[[v]] <- log_d(d[[v]])
  r2 <- function(y, rhs) {
    f <- stats::reformulate(rhs, response = y)
    summary(stats::lm(f, data = d))$r.squared
  }
  out <- do.call(rbind, lapply(size_vars, function(v) {
    a <- r2(v, base); b <- r2(v, c(base, extra))
    data.frame(size_var = v, r2_terrain_gm = round(a, 3), vif_terrain_gm = round(1 / (1 - a), 2),
               r2_plus_detect = round(b, 3), vif_plus_detect = round(1 / (1 - b), 2))
  }))
  out
}