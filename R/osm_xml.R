# R/osm_xml.R -- direct reader for Overpass XML (replaces GDAL OSM driver, which dropped >99% of ways)
# Needs: detectability.R (HWY_MAJOR, HWY_MINOR, PLACE_SETTLE)

.xml_attr <- function(x, attr) sub(paste0('.* ', attr, '="([^"]*)".*'), "\\1", x)

# Ways: lines between <way ...> and </way>; nodes carry lat/lon. Tags inside nodes are ignored.
parse_osm_roads_xml <- function(path, keep = c(HWY_MAJOR, HWY_MINOR)) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  is_node  <- grepl("<node ", x, fixed = TRUE)
  is_wopen <- grepl("<way ", x, fixed = TRUE)
  is_wclose <- grepl("</way>", x, fixed = TRUE)
  way_idx <- cumsum(is_wopen)
  in_way  <- (cumsum(is_wopen) - cumsum(c(FALSE, is_wclose[-length(is_wclose)]))) > 0
  
  nl  <- x[is_node]
  nid <- as.numeric(.xml_attr(nl, "id"))
  nxy <- cbind(lon = as.numeric(.xml_attr(nl, "lon")), lat = as.numeric(.xml_attr(nl, "lat")))
  n_way <- sum(is_wopen)
  if (n_way == 0L) return(NULL)
  way_id <- .xml_attr(x[is_wopen], "id")
  
  hwy <- rep(NA_character_, n_way)
  ti  <- which(in_way & grepl(' k="highway" ', x, fixed = TRUE))
  hwy[way_idx[ti]] <- .xml_attr(x[ti], "v")
  
  ni  <- which(in_way & grepl("<nd ", x, fixed = TRUE))
  ref <- as.numeric(.xml_attr(x[ni], "ref"))
  w   <- way_idx[ni]
  m   <- match(ref, nid)
  ok  <- !is.na(m)
  sp  <- split(m[ok], w[ok])
  widx <- as.integer(names(sp))
  keep_w <- !is.na(hwy[widx]) & hwy[widx] %in% keep & lengths(sp) >= 2L
  sp <- sp[keep_w]; widx <- widx[keep_w]
  if (!length(sp)) return(NULL)
  geom <- lapply(sp, function(i) sf::st_linestring(nxy[i, , drop = FALSE]))
  sf::st_sf(osm_id = way_id[widx], highway = hwy[widx],
            geometry = sf::st_sfc(geom, crs = 4326))
}

# Place nodes: <node ...> followed by <tag k="place" v="..."/> lines
parse_osm_places_xml <- function(path, keep = PLACE_SETTLE) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  is_node <- grepl("<node ", x, fixed = TRUE)
  node_idx <- cumsum(is_node)
  ti <- which(grepl(' k="place" ', x, fixed = TRUE) & node_idx > 0L)
  if (!length(ti)) return(NULL)
  ni <- node_idx[ti]
  nl <- x[is_node]
  place <- .xml_attr(x[ti], "v")
  sel <- place %in% keep
  ni <- ni[sel]; place <- place[sel]
  if (!length(ni)) return(NULL)
  d <- data.frame(osm_id = .xml_attr(nl[ni], "id"), place = place,
                  lon = as.numeric(.xml_attr(nl[ni], "lon")), lat = as.numeric(.xml_attr(nl[ni], "lat")),
                  stringsAsFactors = FALSE)
  sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
}

read_osm_xml_tiles <- function(files, parser) {
  parts <- Filter(Negate(is.null), lapply(files, parser))
  if (!length(parts)) return(NULL)
  x <- do.call(rbind, parts)
  x[!duplicated(x$osm_id), ]
}