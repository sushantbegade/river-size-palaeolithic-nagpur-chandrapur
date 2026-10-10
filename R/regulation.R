CAP_FALLBACK_MCM <- 50
MATCH_MAX_M      <- 1500
EXCL_TIMELINE    <- c("Planned", "Under construction", "Removed", "Destroyed", "Subsumed", "Replaced")

REG_SCENARIOS <- data.frame(
  scen       = 1:8,
  label      = c("primary_dor10_orderstop_50km", "dor10_orderstop_10km", "dor10_orderstop_100km",
                 "dor10_noorderstop_50km", "dor10_noorderstop_100km", "dor2_orderstop_50km",
                 "dor50_orderstop_50km", "allinstream_orderstop_50km"),
  dor_min    = c(10, 10, 10, 10, 10, 2, 50, NA),
  dmax_km    = c(50, 10, 100, 50, 100, 50, 50, 50),
  order_stop = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

select_barriers <- function(bd, dor_min) {
  tl  <- as.character(bd$TIMELINE)
  ok  <- (bd$INSTREAM %in% "Instream") & !(tl %in% EXCL_TIMELINE)
  dor <- ifelse(bd$DOR_PC < 0, NA_real_, bd$DOR_PC)
  cap <- ifelse(bd$CAP_MCM < 0, NA_real_, bd$CAP_MCM)
  if (is.na(dor_min)) return(which(ok))
  strong <- (!is.na(dor) & dor >= dor_min) | (is.na(dor) & !is.na(cap) & cap >= CAP_FALLBACK_MCM)
  which(ok & strong)
}

trace_regulated <- function(hr, bd, rows, dmax_km, order_stop) {
  ids <- hr$HYRIV_ID; nxt <- match(hr$NEXT_DOWN, ids); len <- hr$LENGTH_KM; ord <- hr$ORD_STRA
  start <- match(bd$HYRIV_ID[rows], ids)
  oi <- integer(0); od <- numeric(0); og <- numeric(0)
  for (k in seq_along(rows)) {
    s <- start[k]
    if (is.na(s)) next
    cur <- s; cum <- 0.5 * len[s]; steps <- 0L
    oi <- c(oi, s); od <- c(od, 0); og <- c(og, bd$GDW_ID[rows[k]])
    repeat {
      nx <- nxt[cur]; steps <- steps + 1L
      if (is.na(nx) || cum > dmax_km || steps > 5000L) break
      if (order_stop && ord[nx] > ord[s]) break
      oi <- c(oi, nx); od <- c(od, cum); og <- c(og, bd$GDW_ID[rows[k]])
      cum <- cum + len[nx]; cur <- nx
    }
  }
  all <- data.frame(idx = oi, dist_km = od, GDW_ID = og)
  by_dam <- if (nrow(all)) {
    stats::aggregate(list(n_reaches = all$idx), list(GDW_ID = all$GDW_ID), length)
  } else data.frame(GDW_ID = numeric(0), n_reaches = integer(0))
  if (nrow(all)) {
    by_dam$km_flagged <- round(vapply(by_dam$GDW_ID, function(g) sum(len[all$idx[all$GDW_ID == g]]), numeric(1)), 1)
  } else by_dam$km_flagged <- numeric(0)
  all <- all[order(all$idx, all$dist_km), ]
  list(by_reach = all[!duplicated(all$idx), ], by_dam = by_dam,
       start_found = sum(!is.na(start)), start_missing = sum(is.na(start)))
}

match_to_hr <- function(pts_sf, hr_sf, max_m = MATCH_MAX_M) {
  i <- sf::st_nearest_feature(pts_sf, hr_sf)
  d <- as.numeric(sf::st_distance(pts_sf, hr_sf[i, ], by_element = TRUE))
  i[d > max_m] <- NA_integer_
  list(idx = i, dist_m = d)
}

mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

top_names <- function(names_vec, weights, n = 3) {
  ok <- !is.na(names_vec)
  if (!any(ok)) return(NA_character_)
  w <- tapply(weights[ok], names_vec[ok], sum)
  paste(names(sort(w, decreasing = TRUE))[seq_len(min(n, length(w)))], collapse = "; ")
}