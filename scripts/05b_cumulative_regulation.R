suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "persistence.R"))
source(here::here("R", "regulation.R"))
source(here::here("R", "dor.R"))

t0 <- Sys.time()
message("== S05b cumulative regulation (DOR) ==")
f_e <- P("data", "processed", "external_region.gpkg")
f_s <- P("data", "processed", "pi_subsegments.csv")
f_r <- P("data", "processed", "regulation_flags.csv")
for (f in c(f_e, f_s, f_r)) if (!file.exists(f)) stop("Missing: ", f, "\n  Run S03, S04, S05 first.", call. = FALSE)

hr <- sf::st_drop_geometry(sf::st_read(f_e, layer = "hydrorivers", quiet = TRUE))
bd <- sf::st_drop_geometry(sf::st_read(f_e, layer = "gdw_barriers", quiet = TRUE))
hr$HYRIV_ID <- as.numeric(hr$HYRIV_ID); hr$NEXT_DOWN <- as.numeric(hr$NEXT_DOWN)
bd$HYRIV_ID <- as.numeric(bd$HYRIV_ID); bd$GDW_ID <- as.numeric(bd$GDW_ID)
dam_name <- dplyr::coalesce(as.character(bd$DAM_NAME), as.character(bd$RES_NAME), paste0("GDW_", bd$GDW_ID))

nxt <- match(hr$NEXT_DOWN, hr$HYRIV_ID)
viol <- sum(!is.na(nxt) & hr$UPLAND_SKM[nxt] < hr$UPLAND_SKM)
message("reaches: ", nrow(hr), " | with a downstream reach in region: ", sum(!is.na(nxt)),
        " | upland-area order violations: ", viol)
if (viol > 0) warning(viol, " reach pairs have downstream upland < upstream upland; check network.", call. = FALSE)

rows <- dam_rows_for_dor(bd)
est  <- is_estimated_capacity(bd, rows)
message("instream dams with usable capacity: ", length(rows), " of ", sum(bd$INSTREAM %in% "Instream"),
        " instream | estimated-capacity dams: ", sum(est))
st_all <- cumulative_storage(hr, bd, rows)
st_rep <- cumulative_storage(hr, bd, rows[!est])
hr$stor_all_mcm <- st_all; hr$stor_rep_mcm <- st_rep
hr$dor_all <- dor_index(hr, st_all)
hr$dor_rep <- dor_index(hr, st_rep)

subs <- readr::read_csv(f_s, show_col_types = FALSE)
flg  <- readr::read_csv(f_r, show_col_types = FALSE)
subs <- dplyr::left_join(subs, flg[, c("sub_id", "hr_id", "reg_1")], by = "sub_id")
k <- match(subs$hr_id, hr$HYRIV_ID)
subs$dor_all   <- hr$dor_all[k]
subs$dor_rep   <- hr$dor_rep[k]
subs$stor_mcm  <- hr$stor_all_mcm[k]
subs$dis_cms   <- hr$DIS_AV_CMS[k]
subs$len_total_km <- subs$n_total * SPACING_M / 1000
subs$rivname_f <- dplyr::coalesce(subs$rivname, "<unnamed>")
big <- subs[!is.na(subs$ordsh) & subs$ordsh >= 6L, ]

sc <- rbind(
  data.frame(label = "S05_rule_distance_order_stop", var = "reg_1", thr = NA_real_),
  data.frame(label = paste0("dor_all_ge_", DOR_THRESHOLDS), var = "dor_all", thr = DOR_THRESHOLDS),
  data.frame(label = "dor_reported_only_ge_5", var = "dor_rep", thr = 5))
sc$primary <- sc$label == paste0("dor_all_ge_", DOR_PRIMARY)

use <- big[big$usable %in% TRUE, ]
gate <- do.call(rbind, lapply(seq_len(nrow(sc)), function(i) {
  reg <- if (sc$var[i] == "reg_1") use$reg_1 else !is.na(use[[sc$var[i]]]) & use[[sc$var[i]]] >= sc$thr[i]
  reg[is.na(reg)] <- FALSE
  cls <- ifelse(reg, "I_regulated", use$pi_class)
  km  <- use$length_km_valid
  by_r <- function(kk) { x <- tapply(km[cls == kk], use$rivname_f[cls == kk], sum); sum(x >= 10, na.rm = TRUE) }
  P_km <- sum(km[cls == "P_candidate"]); S_km <- sum(km[cls == "S"])
  data.frame(label = sc$label[i], primary = sc$primary[i], P_km = round(P_km, 1), S_km = round(S_km, 1),
             I_km = round(sum(km[cls == "I"]), 1), I_regulated_km = round(sum(km[cls == "I_regulated"]), 1),
             named_rivers_P_ge10km = by_r("P_candidate"), named_rivers_S_ge10km = by_r("S"),
             gate0b_pass = P_km >= 150 & S_km >= 150 & by_r("P_candidate") >= 3 & by_r("S") >= 3)
}))
write_csv_safe(gate, P("output", "tables", "S05b_gate0b_dor_scenarios.csv"))

share <- function(flag, w) round(sum(w[flag %in% TRUE]) / sum(w), 2)
by_river <- big |>
  dplyr::group_by(rivname_f) |>
  dplyr::summarise(max_ordsh = max(ordsh), km = round(sum(len_total_km), 1),
                   median_dis_cms = round(stats::median(dis_cms, na.rm = TRUE), 1),
                   dor_p10 = round(stats::quantile(dor_all, .1, na.rm = TRUE), 1),
                   dor_median = round(stats::median(dor_all, na.rm = TRUE), 1),
                   dor_p90 = round(stats::quantile(dor_all, .9, na.rm = TRUE), 1),
                   max_storage_mcm = round(max(stor_mcm, na.rm = TRUE)),
                   share_ge2 = share(dor_all >= 2, len_total_km), share_ge5 = share(dor_all >= 5, len_total_km),
                   share_ge10 = share(dor_all >= 10, len_total_km), share_ge25 = share(dor_all >= 25, len_total_km),
                   share_ge50 = share(dor_all >= 50, len_total_km),
                   share_ge5_reported_only = share(dor_rep >= 5, len_total_km), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(max_ordsh), dplyr::desc(km))
write_csv_safe(by_river, P("output", "tables", "S05b_dor_by_river_order6plus.csv"))

prim_reg <- !is.na(use$dor_all) & use$dor_all >= DOR_PRIMARY
use$class_nat <- ifelse(prim_reg, "I_regulated", use$pi_class)
cls_river <- use |>
  dplyr::group_by(rivname_f, class_nat) |>
  dplyr::summarise(km = round(sum(length_km_valid), 1), .groups = "drop") |>
  tidyr::pivot_wider(names_from = class_nat, values_from = km, values_fill = 0)
write_csv_safe(cls_river, P("output", "tables", "S05b_class_by_river_primary.csv"))

write_csv_safe(subs[, c("sub_id", "hr_id", "dor_all", "dor_rep", "stor_mcm", "dis_cms")],
               P("data", "processed", "regulation_dor.csv"))
write_csv_safe(data.frame(item = c("DOR_PRIMARY_pct", "DOR_THRESHOLDS_pct", "MCM_PER_CMS_YR",
                                   "n_dams_used", "n_estimated_capacity"),
                          value = c(DOR_PRIMARY, paste(DOR_THRESHOLDS, collapse = ","), MCM_PER_CMS_YR,
                                    length(rows), sum(est))),
               P("output", "logs", "S05b_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S05b_sessionInfo.txt"))

cat("\n========== S05b SUMMARY ==========\n")
cat("Dams used:", length(rows), "| estimated-capacity dams:", sum(est), "\n")
cat("Sub-segments with a DOR value:", sum(!is.na(subs$dor_all)), "of", nrow(subs), "\n")
cat("\nDOR (%) distribution on HydroRIVERS reaches intersecting the districts:\n")
inst <- hr$in_study %in% TRUE
print(round(stats::quantile(hr$dor_all[inst], c(.1, .25, .5, .75, .9, .99)), 1))
cat("\nGate 0b preview under cumulative-DOR rules (PROVISIONAL thresholds):\n")
print(as.data.frame(gate), row.names = FALSE)
cat("\nDOR by river (order >= 6): percentiles of cumulative DOR along the river, natural discharge, shares regulated:\n")
print(as.data.frame(utils::head(by_river, 35)), row.names = FALSE)
cat("\nClass length (km) by river, PRIMARY rule (DOR >=", DOR_PRIMARY, "% => I_regulated):\n")
ord_r <- if ("P_candidate" %in% names(cls_river)) order(-cls_river$P_candidate) else seq_len(nrow(cls_river))
print(as.data.frame(utils::head(cls_river[ord_r, ], 25)), row.names = FALSE)
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S05b. Paste this console output back.\n")