# R/clean_sites.R -- functions for script 01 (no data are read at source time)

HEADER_PREFIX <- c("site id", "site name", "latitude", "longitude", "district", "state",
                   "associated river", "river basin", "location precision",
                   "lithic assemblage", "site type", "raw material",
                   "reported surface context", "cultural attribution",
                   "relative chronology", "publication year", "primary reference",
                   "data source type")

COL_NAMES <- c("site_id", "site_name", "lat", "lon", "district", "state", "river",
               "basin", "precision", "assemblage", "site_type_raw", "raw_material",
               "context", "culture", "chron", "pub_year", "ref", "source_type")

read_sites_raw <- function(path) {
  raw <- readxl::read_excel(path, .name_repair = "minimal")
  got <- tolower(stringr::str_squish(names(raw)))
  if (length(got) != length(HEADER_PREFIX)) {
    stop("Expected ", length(HEADER_PREFIX), " columns, found ", length(got), ".", call. = FALSE)
  }
  ok <- startsWith(got, HEADER_PREFIX)
  if (!all(ok)) {
    stop("Header mismatch at column(s): ", paste(which(!ok), collapse = ", "),
         "\n  found: ", paste(got[!ok], collapse = " | "), call. = FALSE)
  }
  names(raw) <- COL_NAMES
  raw
}

na_codes <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x %in% c("", "-", "Not Reported")] <- NA_character_
  x
}

clean_sites <- function(raw, crosswalk) {
  stopifnot(all(c("chron", "LP", "MP", "UP") %in% names(crosswalk)))
  stopifnot(!anyDuplicated(crosswalk$chron))
  d <- as.data.frame(raw, stringsAsFactors = FALSE)

  chr <- names(d)[vapply(d, is.character, logical(1))]
  d[chr] <- lapply(d[chr], stringr::str_squish)
  d$site_id  <- as.integer(d$site_id)
  d$lat      <- as.numeric(d$lat)
  d$lon      <- as.numeric(d$lon)
  d$pub_year <- as.character(d$pub_year)
  for (v in c("river", "raw_material", "context", "culture", "assemblage", "source_type")) {
    d[[v]] <- na_codes(d[[v]])
  }

  # --- site type ---------------------------------------------------------------
  lt <- stringr::str_to_lower(d$site_type_raw)
  d$site_type <- dplyr::case_when(
    lt == "open-air"                     ~ "open_air",
    lt == "riverbed"                     ~ "riverbed",
    lt == "lake shore"                   ~ "lake_shore",
    stringr::str_detect(lt, "trench")    ~ "insitu_trench",
    lt == "factory"                      ~ "factory",
    lt == "not reported"                 ~ "not_reported",
    TRUE                                 ~ NA_character_
  )
  if (anyNA(d$site_type)) {
    stop("Unmapped site types: ",
         paste(unique(d$site_type_raw[is.na(d$site_type)]), collapse = "; "), call. = FALSE)
  }

  # --- coordinate tier ---------------------------------------------------------
  d$coord_tier <- dplyr::case_when(
    stringr::str_detect(d$precision, "Field Verified")                 ~ "A_field_verified",
    stringr::str_detect(d$precision, "Reported By Referenced Author")  ~ "B_author_reported",
    stringr::str_detect(d$precision, "GIS Map Based")                  ~ "C_map_based",
    TRUE                                                               ~ NA_character_
  )
  if (anyNA(d$coord_tier)) {
    stop("Unmapped precision terms: ",
         paste(unique(d$precision[is.na(d$coord_tier)]), collapse = "; "), call. = FALSE)
  }

  # --- period via crosswalk (docs/crosswalk_period.csv) ------------------------
  cw <- as.data.frame(crosswalk[, c("chron", "LP", "MP", "UP")])
  cw$chron <- stringr::str_squish(cw$chron)
  for (p in c("LP", "MP", "UP")) cw[[p]] <- as.integer(cw[[p]])
  bad <- setdiff(unique(d$chron), cw$chron)
  if (length(bad)) {
    stop("Chronology terms missing from docs/crosswalk_period.csv:\n  ",
         paste(bad, collapse = "\n  "), call. = FALSE)
  }
  d <- dplyr::left_join(d, cw, by = "chron")
  d$n_periods     <- d$LP + d$MP + d$UP
  d$single_period <- d$n_periods == 1L
  d$period_class  <- dplyr::case_when(
    d$n_periods > 1L ~ "multi",
    d$LP == 1L       ~ "LP",
    d$MP == 1L       ~ "MP",
    d$UP == 1L       ~ "UP",
    TRUE             ~ NA_character_
  )

  # --- pre-registered sample sets (design v3.1, Part 4) ------------------------
  d$S1_openair  <- d$site_type == "open_air"
  d$S2_all      <- TRUE
  d$S3_riverbed <- d$site_type == "riverbed"
  d$S4_other    <- !(d$S1_openair | d$S3_riverbed)

  # --- flags and cluster id ----------------------------------------------------
  d$flag_megalith <- stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(d$context, "")),
                                         "megalith|cairn")
  d$ref_clean <- dplyr::coalesce(d$ref, "unknown")
  d$ref_id    <- as.integer(factor(d$ref_clean))
  d
}

add_utm <- function(d) {
  pts <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = CRS_GEOGRAPHIC, remove = FALSE)
  pts <- sf::st_transform(pts, CRS_ANALYSIS)
  xy  <- sf::st_coordinates(pts)
  d$x <- xy[, 1]
  d$y <- xy[, 2]
  d
}

validate_sites <- function(d) {
  stopifnot(!anyNA(d$site_id), !anyDuplicated(d$site_id))
  stopifnot(!anyNA(d$lat), !anyNA(d$lon))
  inside <- d$lon >= STUDY_BBOX["xmin"] & d$lon <= STUDY_BBOX["xmax"] &
            d$lat >= STUDY_BBOX["ymin"] & d$lat <= STUDY_BBOX["ymax"]
  if (!all(inside)) {
    stop("Sites outside study bounding box: site_id ",
         paste(d$site_id[!inside], collapse = ", "), call. = FALSE)
  }
  stopifnot(!anyNA(d$coord_tier), !anyNA(d$period_class), !anyNA(d$site_type))
  stopifnot(all(d$district %in% c("Nagpur", "Chandrapur")))
  dup <- sum(duplicated(d[c("lat", "lon")]))
  if (dup > 0) message("NOTE: ", dup, " site(s) share identical coordinates with another site.")
  invisible(TRUE)
}

# ---- thinning counts (numbers only; no thinned dataset is saved) ---------------
greedy_thin_n <- function(xy, d_m, seed) {
  n <- nrow(xy)
  if (n == 0L) return(0L)
  set.seed(seed)
  keep <- integer(0)
  for (i in sample.int(n)) {
    if (!length(keep) ||
        all(sqrt((xy[keep, 1] - xy[i, 1])^2 + (xy[keep, 2] - xy[i, 2])^2) >= d_m)) {
      keep <- c(keep, i)
    }
  }
  length(keep)
}

thinning_table <- function(d, dists_m = c(2000, 5000), n_rep = 50L, seed = 1L) {
  groups <- list(
    pooled     = rep(TRUE, nrow(d)),
    S1_openair = d$S1_openair,
    LP_any     = d$LP == 1L,
    MP_any     = d$MP == 1L,
    UP_any     = d$UP == 1L,
    LP_only    = d$period_class == "LP",
    MP_only    = d$period_class == "MP",
    UP_only    = d$period_class == "UP"
  )
  out <- expand.grid(group = names(groups), d_km = dists_m / 1000, stringsAsFactors = FALSE)
  out$n_raw <- vapply(out$group, function(g) sum(groups[[g]]), numeric(1))
  out$n_thinned_median <- mapply(function(g, dk) {
    xy <- as.matrix(d[groups[[g]], c("x", "y")])
    stats::median(vapply(seq_len(n_rep), function(r) greedy_thin_n(xy, dk * 1000, seed + r),
                         numeric(1)))
  }, out$group, out$d_km)
  out
}

# ---- coarsened public version ---------------------------------------------------
coarsen_sites <- function(d, cell_m = 1000) {
  d$x_c <- (floor(d$x / cell_m) + 0.5) * cell_m
  d$y_c <- (floor(d$y / cell_m) + 0.5) * cell_m
  ll <- sf::st_as_sf(d[c("site_id", "x_c", "y_c")], coords = c("x_c", "y_c"),
                     crs = CRS_ANALYSIS, remove = FALSE)
  ll <- sf::st_transform(ll, CRS_GEOGRAPHIC)
  xy <- sf::st_coordinates(ll)
  d$lon_c <- round(xy[, 1], 5)
  d$lat_c <- round(xy[, 2], 5)
  keep <- c("site_id", "district", "basin", "coord_tier", "site_type", "LP", "MP", "UP",
            "period_class", "S1_openair", "S3_riverbed", "S4_other", "river", "ref_id",
            "x_c", "y_c", "lon_c", "lat_c")
  d[keep]
}

# ---- summary tables (contain NO coordinates; safe to commit) ---------------------
summarise_sites <- function(d) {
  list(
    S01_type_by_district = tidyr::pivot_wider(dplyr::count(d, site_type, district),
                             names_from = district, values_from = n, values_fill = 0),
    S01_tier_by_period = dplyr::summarise(dplyr::group_by(d, coord_tier),
                           n = dplyr::n(), LP = sum(LP), MP = sum(MP), UP = sum(UP),
                           single_period = sum(single_period)),
    S01_period_class_by_set = rbind(
      cbind(set = "S1_openair",  dplyr::count(d[d$S1_openair, ],  period_class)),
      cbind(set = "S2_all",      dplyr::count(d,                  period_class)),
      cbind(set = "S3_riverbed", dplyr::count(d[d$S3_riverbed, ], period_class)),
      cbind(set = "S4_other",    dplyr::count(d[d$S4_other, ],    period_class))),
    S01_tier_by_type = tidyr::pivot_wider(dplyr::count(d, coord_tier, site_type),
                         names_from = site_type, values_from = n, values_fill = 0),
    S01_reference_counts = dplyr::summarise(dplyr::group_by(d, ref_id, ref_clean),
                             n = dplyr::n(), n_open_air = sum(S1_openair),
                             n_riverbed = sum(S3_riverbed), .groups = "drop") |>
                           dplyr::arrange(dplyr::desc(n)),
    S01_culture_vs_chron_QA = dplyr::arrange(dplyr::count(d, culture, chron), culture, chron),
    S01_flag_megalith = d[d$flag_megalith, c("site_id", "site_type", "period_class")]
  )
}
