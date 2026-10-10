suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "persistence.R"))
source(here::here("R", "regulation.R"))

t0 <- Sys.time()
message("== S05 regulation flags ==")
f_h <- P("data", "processed", "hydro_prep.gpkg")
f_e <- P("data", "processed", "external_region.gpkg")
f_p <- P("data", "processed", "pi_points.rds")
f_s <- P("data", "processed", "pi_subsegments.csv")
for (f in c(f_h, f_e, f_p, f_s)) if (!file.exists(f)) stop("Missing: ", f, "\n  Run S02, S03, S04 first.", call. = FALSE)

bnd  <- sf::st_read(f_h, layer = "boundary", quiet = TRUE)
bu   <- sf::st_union(sf::st_geometry(bnd))
hr   <- sf::st_read(f_e, layer = "hydrorivers", quiet = TRUE)
barr <- sf::st_read(f_e, layer = "gdw_barriers", quiet = TRUE)
resv <- sf::st_read(f_e, layer = "gdw_reservoirs", quiet = TRUE)
hr$HYRIV_ID <- as.numeric(hr$HYRIV_ID); hr$NEXT_DOWN <- as.numeric(hr$NEXT_DOWN)
bd <- sf::st_drop_geometry(barr)
bd$HYRIV_ID <- as.numeric(bd$HYRIV_ID); bd$GDW_ID <- as.numeric(bd$GDW_ID)
dam_name <- dplyr::coalesce(as.character(bd$DAM_NAME), as.character(bd$RES_NAME), paste0("GDW_", bd$GDW_ID))
message("HydroRIVERS reaches: ", nrow(hr), " | GDW barriers: ", nrow(barr),
        " | instream: ", sum(bd$INSTREAM %in% "Instream"))

P3   <- readRDS(f_p)$points
P3$sub_id <- paste(P3$seg_id, P3$chunk, sep = "_")
sub  <- readr::read_csv(f_s, show_col_types = FALSE)
crd  <- stats::aggregate(cbind(x, y) ~ sub_id, P3, mean)
bpts <- sf::st_as_sf(P3[, c("sub_id", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS)
near_b <- lengths(sf::st_is_within_distance(bpts, sf::st_boundary(bu), dist = 1000)) > 0
nb <- stats::aggregate(list(near_bnd_share = near_b), list(sub_id = P3$sub_id), mean)
sub <- sub |>
  dplyr::left_join(crd, by = "sub_id") |>
  dplyr::left_join(nb, by = "sub_id")
sp <- sf::st_as_sf(sub[, c("sub_id", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS)

m <- match_to_hr(sp, hr)
sub$hr_id      <- ifelse(is.na(m$idx), NA_real_, hr$HYRIV_ID[m$idx])
sub$hr_dist_m  <- round(m$dist_m)
sub$hr_upland  <- ifelse(is.na(m$idx), NA_real_, hr$UPLAND_SKM[m$idx])
sub$hr_order   <- ifelse(is.na(m$idx), NA_integer_, hr$ORD_STRA[m$idx])
message("sub-segments: ", nrow(sub), " | matched to HydroRIVERS within ", MATCH_MAX_M, " m: ",
        sum(!is.na(sub$hr_id)))

reg_ids <- list(); dams_primary <- NULL
for (s in REG_SCENARIOS$scen) {
  sc   <- REG_SCENARIOS[s, ]
  rows <- select_barriers(bd, sc$dor_min)
  tr   <- trace_regulated(hr, bd, rows, sc$dmax_km, sc$order_stop)
  reg_ids[[s]] <- hr$HYRIV_ID[tr$by_reach$idx]
  sub[[paste0("reg_", s)]] <- !is.na(sub$hr_id) & sub$hr_id %in% reg_ids[[s]]
  message(sprintf("scenario %d %-30s barriers %4d | start reach found %4d | flagged HR reaches %5d",
                  s, sc$label, length(rows), tr$start_found, nrow(tr$by_reach)))
  if (s == 1) {
    prim <- tr
    flag_tab <- tr$by_reach
    sub$dam_id_primary <- flag_tab$GDW_ID[match(sub$hr_id, hr$HYRIV_ID[flag_tab$idx])]
    sub$dist_below_km_primary <- round(flag_tab$dist_km[match(sub$hr_id, hr$HYRIV_ID[flag_tab$idx])], 1)
    sub$dam_name_primary <- dam_name[match(sub$dam_id_primary, bd$GDW_ID)]
    rows1 <- rows
  }
}

sub$len_total_km <- sub$n_total * SPACING_M / 1000
sub$rivname_f <- dplyr::coalesce(sub$rivname, "<unnamed>")
big <- sub[!is.na(sub$ordsh) & sub$ordsh >= 6L, ]

mi <- match(bd$HYRIV_ID[rows1], hr$HYRIV_ID)
snap_d <- rep(NA_real_, length(rows1))
okm <- !is.na(mi)
if (any(okm)) snap_d[okm] <- as.numeric(sf::st_distance(barr[rows1[okm], ], hr[mi[okm], ], by_element = TRUE))
dist_study <- as.numeric(sf::st_distance(barr[rows1, ], bu)) / 1000
dams <- data.frame(GDW_ID = bd$GDW_ID[rows1], dam = dam_name[rows1], river = bd$RIVER[rows1],
                   year = ifelse(bd$YEAR_DAM[rows1] < 0, NA, bd$YEAR_DAM[rows1]),
                   cap_mcm = ifelse(bd$CAP_MCM[rows1] < 0, NA, bd$CAP_MCM[rows1]),
                   dor_pc = ifelse(bd$DOR_PC[rows1] < 0, NA, round(bd$DOR_PC[rows1], 1)),
                   quality = bd$QUALITY[rows1], snap_dist_m = round(snap_d),
                   km_from_study = round(dist_study, 1), stringsAsFactors = FALSE)
dams <- dplyr::left_join(dams, prim$by_dam, by = "GDW_ID")
dams$wris_subsegs_flagged <- vapply(dams$GDW_ID, function(g) sum(big$dam_id_primary == g & big$reg_1, na.rm = TRUE), numeric(1))
dams <- dams[order(-dams$wris_subsegs_flagged, dams$km_from_study), ]
write_csv_safe(dams[dams$km_from_study <= 150 | dams$wris_subsegs_flagged > 0, ],
               P("output", "tables", "S05_dams_primary_set.csv"))

by_river <- big |>
  dplyr::group_by(rivname_f) |>
  dplyr::summarise(max_ordsh = max(ordsh), n_sub = dplyr::n(), km = round(sum(len_total_km), 1),
                   reg_km_primary = round(sum(len_total_km[reg_1]), 1),
                   reg_share_primary = round(sum(len_total_km[reg_1]) / sum(len_total_km), 2),
                   reg_share_10km = round(sum(len_total_km[reg_2]) / sum(len_total_km), 2),
                   reg_share_100km = round(sum(len_total_km[reg_3]) / sum(len_total_km), 2),
                   reg_share_noorder100 = round(sum(len_total_km[reg_5]) / sum(len_total_km), 2),
                   reg_share_dor2 = round(sum(len_total_km[reg_6]) / sum(len_total_km), 2),
                   top_dams = top_names(dam_name_primary, len_total_km),
                   median_hr_dist_m = stats::median(hr_dist_m, na.rm = TRUE),
                   near_boundary_share = round(sum(len_total_km * near_bnd_share, na.rm = TRUE) / sum(len_total_km), 2),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(max_ordsh), dplyr::desc(km))
write_csv_safe(by_river, P("output", "tables", "S05_regulation_by_river_order6plus.csv"))

use <- big[big$usable %in% TRUE, ]
gate <- do.call(rbind, lapply(REG_SCENARIOS$scen, function(s) {
  reg <- use[[paste0("reg_", s)]]
  cls <- ifelse(reg, "I_regulated", use$pi_class)
  km  <- use$length_km_valid
  by_r <- function(k) { x <- tapply(km[cls == k], use$rivname_f[cls == k], sum); sum(x >= 10, na.rm = TRUE) }
  P_km <- sum(km[cls == "P_candidate"]); S_km <- sum(km[cls == "S"])
  data.frame(scen = s, label = REG_SCENARIOS$label[s], P_km = round(P_km, 1), S_km = round(S_km, 1),
             I_km = round(sum(km[cls == "I"]), 1), I_regulated_km = round(sum(km[cls == "I_regulated"]), 1),
             named_rivers_P_ge10km = by_r("P_candidate"), named_rivers_S_ge10km = by_r("S"),
             gate0b_pass = P_km >= 150 & S_km >= 150 & by_r("P_candidate") >= 3 & by_r("S") >= 3)
}))
write_csv_safe(gate, P("output", "tables", "S05_gate0b_preview_by_scenario.csv"))

cls_tab <- use |>
  dplyr::mutate(og = order_group(ordsh), class = ifelse(reg_1, "I_regulated", pi_class)) |>
  dplyr::group_by(og, class) |>
  dplyr::summarise(n_sub = dplyr::n(), length_km = round(sum(length_km_valid), 1), .groups = "drop")
write_csv_safe(cls_tab, P("output", "tables", "S05_class_after_regulation.csv"))

cal <- sub[!is.na(sub$hr_id) & !is.na(sub$ordsh), ] |>
  dplyr::group_by(ordsh) |>
  dplyr::summarise(n_sub = dplyr::n(), upland_q10 = round(stats::quantile(hr_upland, .1)),
                   upland_median = round(stats::median(hr_upland)),
                   upland_q90 = round(stats::quantile(hr_upland, .9)),
                   hr_order_mode = mode_chr(as.character(hr_order)), .groups = "drop")
write_csv_safe(cal, P("output", "tables", "S05_wris_order_vs_hr_upland.csv"))

named <- sub[!is.na(sub$rivname) & !is.na(sub$ordsh) & sub$ordsh >= 8L, ]
un    <- sub[is.na(sub$rivname) & !is.na(sub$ordsh) & sub$ordsh >= 8L, ]
if (nrow(un) && nrow(named)) {
  ns <- sf::st_as_sf(named[, c("x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS)
  us <- sf::st_as_sf(un[, c("x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS)
  ni <- sf::st_nearest_feature(us, ns)
  un$nearest_named <- named$rivname[ni]
  un$dist_named_m  <- round(as.numeric(sf::st_distance(us, ns[ni, ], by_element = TRUE)))
  un$near_resv     <- lengths(sf::st_is_within_distance(us, resv, dist = 500)) > 0
  unt <- un |>
    dplyr::group_by(nearest_named) |>
    dplyr::summarise(n_sub = dplyr::n(), km = round(sum(len_total_km), 1),
                     median_pi = round(stats::median(pi, na.rm = TRUE), 2),
                     share_near_reservoir = round(mean(near_resv), 2),
                     median_dist_to_named_m = stats::median(dist_named_m),
                     share_regulated_primary = round(mean(reg_1), 2), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(km))
} else unt <- data.frame()
write_csv_safe(unt, P("output", "tables", "S05_unnamed_order8plus_reaches.csv"))

keep <- c("sub_id", "hr_id", "hr_dist_m", paste0("reg_", REG_SCENARIOS$scen),
          "dam_id_primary", "dist_below_km_primary", "near_bnd_share")
write_csv_safe(sub[, keep], P("data", "processed", "regulation_flags.csv"))
write_csv_safe(cbind(REG_SCENARIOS, CAP_FALLBACK_MCM = CAP_FALLBACK_MCM, MATCH_MAX_M = MATCH_MAX_M),
               P("output", "logs", "S05_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S05_sessionInfo.txt"))

cat("\n========== S05 SUMMARY ==========\n")
cat("Barriers (primary set, DOR >= 10 % or capacity fallback):", length(rows1),
    "| start reach found:", prim$start_found, "| missing:", prim$start_missing, "\n")
cat("Barrier snap distance to its HydroRIVERS reach (m): median",
    stats::median(dams$snap_dist_m, na.rm = TRUE), "| 90th pct",
    round(stats::quantile(dams$snap_dist_m, .9, na.rm = TRUE)), "\n")
cat("Sub-segment match distance to HydroRIVERS (m): median", stats::median(sub$hr_dist_m, na.rm = TRUE),
    "| 90th pct", round(stats::quantile(sub$hr_dist_m, .9, na.rm = TRUE)), "| unmatched", sum(is.na(sub$hr_id)), "\n")
cat("\nDams that flag WRIS large sub-segments (primary scenario):\n")
print(as.data.frame(utils::head(dams[dams$wris_subsegs_flagged > 0, c("dam", "river", "year", "cap_mcm",
                                                                      "dor_pc", "km_from_study", "km_flagged", "wris_subsegs_flagged")], 40)), row.names = FALSE)
cat("\nBARRIERS NEAR THE DISTRICTS with capacity >= 100 MCM (check names against what you know):\n")
print(as.data.frame(dams[!is.na(dams$cap_mcm) & dams$cap_mcm >= 100 & dams$km_from_study <= 100,
                         c("dam", "river", "year", "cap_mcm", "dor_pc", "km_from_study", "wris_subsegs_flagged")]), row.names = FALSE)
cat("\nRegulation by river (order >= 6):\n")
print(as.data.frame(utils::head(by_river, 35)), row.names = FALSE)
cat("\nGate 0b preview after regulation (PROVISIONAL thresholds):\n")
print(as.data.frame(gate), row.names = FALSE)
cat("\nWRIS order vs HydroRIVERS upland area (km2):\n"); print(as.data.frame(cal), row.names = FALSE)
cat("\nUnnamed order >= 8 reaches (nearest named river):\n"); print(as.data.frame(utils::head(unt, 15)), row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S05. Paste this console output back.\n")