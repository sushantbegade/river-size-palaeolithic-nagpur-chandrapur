DATA <- function(...) file.path(paleo_data_root(), ...)

require_fields <- function(x, fields, layer) {
  miss <- setdiff(fields, names(x))
  if (length(miss)) {
    stop(layer, ": missing field(s): ", paste(miss, collapse = ", "),
         "\n  fields found: ", paste(names(x), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

read_vector_utm <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
  x <- sf::st_read(path, quiet = TRUE)
  if (is.na(sf::st_crs(x))) stop(label, ": layer has no CRS", call. = FALSE)
  x <- sf::st_transform(x, CRS_ANALYSIS)
  if (any(!sf::st_is_valid(x))) x <- sf::st_make_valid(x)
  x[!sf::st_is_empty(x), ]
}

prep_boundary <- function(bnd) {
  chr <- names(bnd)[vapply(sf::st_drop_geometry(bnd), is.character, logical(1))]
  bnd$unit <- if (length(chr)) as.character(bnd[[chr[1]]]) else as.character(seq_len(nrow(bnd)))
  bnd$area_km2 <- as.numeric(sf::st_area(bnd)) / 1e6
  bnd[, c("unit", "area_km2")]
}

prep_rivers <- function(rv, bnd_union, buf_m = 25000) {
  require_fields(rv, c("rivname", "rilcode", "ordsh"), "Rivers")
  rv <- rv[sf::st_geometry_type(rv) %in% c("LINESTRING", "MULTILINESTRING"),
           c("rivname", "rilcode", "ordsh")]
  rv <- suppressWarnings(sf::st_cast(rv, "LINESTRING"))
  rv$rivname <- stringr::str_squish(as.character(rv$rivname))
  rv$rivname[rv$rivname %in% c("", "NA")] <- NA_character_
  rv$rilcode <- ifelse(is.na(rv$rilcode), "NA", as.character(rv$rilcode))
  rv$ordsh   <- as.integer(rv$ordsh)
  rv$seg_id    <- seq_len(nrow(rv))
  rv$length_km <- as.numeric(sf::st_length(rv)) / 1000
  rv$in_study  <- lengths(sf::st_intersects(rv, bnd_union)) > 0
  rv$in_buf25  <- lengths(sf::st_intersects(rv, sf::st_buffer(bnd_union, buf_m))) > 0
  ok <- !is.na(rv$ordsh)
  rv$large5  <- ok & rv$ordsh >= 5L
  rv$large6  <- ok & rv$ordsh >= 6L
  rv$large7  <- ok & rv$ordsh >= 7L
  rv$trunk9  <- ok & rv$ordsh >= 9L
  rv
}

prep_waterbodies <- function(wb, bnd_union) {
  require_fields(wb, c("wbcode", "wbname", "area_ha"), "Waterbody")
  wb <- wb[, c("wbcode", "wbname", "area_ha")]
  wb$wbcode   <- as.character(wb$wbcode)
  wb$wbname   <- stringr::str_squish(as.character(wb$wbname))
  wb$area_km2 <- as.numeric(sf::st_area(wb)) / 1e6
  wb$in_study <- lengths(sf::st_intersects(wb, bnd_union)) > 0
  wb$is_reservoir_code02 <- wb$wbcode == "02"   # inference: confirm against WRIS legend
  wb
}

classify_geomorph <- function(desc) {
  d <- stringr::str_to_lower(stringr::str_squish(as.character(desc)))
  dplyr::case_when(
    stringr::str_detect(d, "^river$")                                    ~ "river_channel",
    stringr::str_detect(d, "dam|reservoir|pond|watbod|waterbod|lake")    ~ "still_water",
    stringr::str_detect(d, "quarry|mine dump|embankment")                ~ "disturbance",
    stringr::str_detect(d, "alluvial|flood plain")                       ~ "alluvial_floodplain",
    stringr::str_detect(d, "pediment|pediplain")                         ~ "pediplain",
    stringr::str_detect(d, "plateau|hills|valley")                       ~ "plateau_hills",
    TRUE                                                                 ~ "other"
  )
}

prep_geomorph <- function(gm) {
  require_fields(gm, "descriptio", "Geomorphology")
  gm <- gm[, "descriptio"]
  gm$descriptio <- stringr::str_squish(as.character(gm$descriptio))
  gm$geom_class <- classify_geomorph(gm$descriptio)
  gm$area_km2   <- as.numeric(sf::st_area(gm)) / 1e6
  gm
}

build_removed_masks <- function(gm, wb, bnd_union, quarry_buffer_m = 100) {
  parts <- list(
    sf::st_geometry(gm[gm$geom_class == "still_water", ]),
    sf::st_geometry(wb[wb$in_study, ]),
    sf::st_buffer(sf::st_geometry(gm[gm$geom_class == "disturbance", ]), quarry_buffer_m)
  )
  parts <- parts[vapply(parts, length, integer(1)) > 0]
  core_u <- sf::st_make_valid(sf::st_union(do.call(c, parts)))
  core <- sf::st_intersection(core_u, bnd_union)
  ch   <- gm[gm$geom_class == "river_channel", ]
  chan <- if (nrow(ch)) sf::st_intersection(sf::st_make_valid(sf::st_union(sf::st_geometry(ch))), bnd_union) else
    sf::st_sfc(crs = CRS_ANALYSIS)
  list(core = sf::st_sf(mask = "core", geometry = core),
       channel = sf::st_sf(mask = "channel", geometry = chan))
}

network_by_order <- function(rv) {
  x <- sf::st_drop_geometry(rv)
  dplyr::summarise(dplyr::group_by(x, ordsh, rilcode, in_study),
                   n = dplyr::n(), length_km = round(sum(length_km), 2), .groups = "drop")
}

large_river_table <- function(rv, min_ord = 6L) {
  x <- sf::st_drop_geometry(rv[rv$in_study & !is.na(rv$ordsh) & rv$ordsh >= min_ord, ])
  x$rivname <- dplyr::coalesce(x$rivname, "<unnamed>")
  x <- dplyr::summarise(dplyr::group_by(x, rivname), max_ordsh = max(ordsh),
                        n_seg = dplyr::n(), length_km = round(sum(length_km), 2), .groups = "drop")
  dplyr::arrange(x, dplyr::desc(max_ordsh), dplyr::desc(length_km))
}

waterbody_by_code <- function(wb) {
  x <- sf::st_drop_geometry(wb)
  dplyr::summarise(dplyr::group_by(x, wbcode, in_study), n = dplyr::n(),
                   area_km2 = round(sum(area_km2), 2), median_ha = round(stats::median(area_ha), 2),
                   max_ha = round(max(area_ha), 1), .groups = "drop")
}

reservoir_table <- function(wb) {
  x <- sf::st_drop_geometry(wb[wb$is_reservoir_code02, ])
  x <- x[order(-x$area_ha), c("wbname", "area_ha", "in_study")]
  x$area_ha <- round(x$area_ha, 1)
  x
}

geomorph_table <- function(gm, study_km2) {
  x <- sf::st_drop_geometry(gm)
  out <- dplyr::summarise(dplyr::group_by(x, geom_class, descriptio), n = dplyr::n(),
                          area_km2 = round(sum(area_km2), 2), .groups = "drop")
  out$pct_of_study_polygons <- round(100 * out$area_km2 / study_km2, 2)
  dplyr::arrange(out, geom_class, dplyr::desc(area_km2))
}