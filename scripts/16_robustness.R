# scripts/16_robustness.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "robustness.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S16 robustness matrix (S1 and S2; baseline = S11/S12 design) ==")

KS       <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)
BLOCK_KM <- as.numeric(Sys.getenv("S16_BLOCK", "40"))
THIN_R   <- as.integer(Sys.getenv("S16_THIN_R", "100"))

# reading rule, fixed BEFORE the result: a variant is "consistent" for S1 k=6 if beta stays positive with a CI that excludes 0, and for
# S2 k=6 if the CI still includes 0. Any variant that changes sign or CI exclusion is reported as sensitive.

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f <- list(q = P("data", "processed", "quadrature_det.rds"), s = P("data", "processed", "sites_cov.rds"), h = P("data", "processed", "hydro_prep.gpkg"),
          e = P("data", "processed", "external_region.gpkg"), cal = P("output", "tables", "S05_wris_order_vs_hr_upland.csv"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S02, S05, S08 and S10 first.", call. = FALSE)
grid <- sf::st_drop_geometry(readRDS(f$q))
sc_all <- readRDS(f$s)
dom  <- grid[grid$in_survey & !grid$in_core, ]
sc   <- sc_all[sc_all$in_domain, ]
mk_sets <- function(d) list(S1_openair = d[d$S1_openair %in% TRUE, ], S2_all = d)
sets <- mk_sets(sc)
ex_d <- function(dm) function(k) exposure_domain(dm, k)
message("domain cells: ", nrow(dom), " | S1 ", nrow(sets$S1_openair), " | S2 ", nrow(sets$S2_all))

res <- list()
add <- function(...) { r <- run_variant(...); res[[length(res) + 1L]] <<- r; message("  done: ", r$variant[1]); invisible(r) }
add("baseline", dom, sets, KS, ex_d(dom), block_km = BLOCK_KM)

# --- pre-registered: without elevation; road variant; prior SD; --------------------------------------------------------------
add("no elevation", dom, sets, KS, ex_d(dom), drop_cols = c("gm_alluvial", "elev"), block_km = BLOCK_KM)
dom_r <- dom; dom_r$d_road_major <- dom_r$d_road_any
sets_r <- lapply(sets, function(d) { d$d_road_major <- d$d_road_any; d })
add("d_road_any instead of d_road_major", dom_r, sets_r, KS, ex_d(dom_r), block_km = BLOCK_KM)
for (psd in c(1, 5, 1e6)) {
  old <- get("PRIOR_SD", envir = globalenv()); assign("PRIOR_SD", psd, envir = globalenv())
  add(paste0("prior SD ", ifelse(psd >= 1e5, "huge (about ML)", psd)), dom, sets, KS, ex_d(dom), block_km = BLOCK_KM)
  assign("PRIOR_SD", old, envir = globalenv())
}

# --- pre-registered: HydroRIVERS extension off (distances to the WRIS network only) -----------------------------------------------
study <- sf::st_union(sf::st_geometry(sf::st_read(f$h, layer = "boundary", quiet = TRUE)))
rv <- sf::st_read(f$h, layer = "rivers", quiet = TRUE); hr <- sf::st_read(f$e, layer = "hydrorivers", quiet = TRUE)
cal <- readr::read_csv(f$cal, show_col_types = FALSE); upl <- stats::setNames(cal$upland_median, as.character(cal$ordsh))
lines <- size_line_sets(rv, hr_outside(hr, study), upl)
spts <- sf::st_as_sf(sc[, c("site_id", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
dom_w <- dom; sets_w <- sets
sc_w <- sc
for (k in KS) { sc_w[[paste0("d_k", k)]] <- nearest_dist(spts, lines$wris_only[[as.character(k)]]); dom_w[[paste0("d_k", k)]] <- dom[[paste0("d_k", k, "_wris")]] }
sets_w <- mk_sets(sc_w)
add("HydroRIVERS extension off (WRIS network only)", dom_w, sets_w, KS, ex_d(dom_w), block_km = BLOCK_KM)

# --- pre-registered: removed core retained in the domain ------------------------------------------------------------------------------
dom_c <- grid[grid$in_survey, ]; dom_c$gm4[is.na(dom_c$gm4)] <- "pediplain"
sc_c <- sc_all[sc_all$in_survey, ]; sc_c$gm4[is.na(sc_c$gm4)] <- "pediplain"
add("removed core retained", dom_c, mk_sets(sc_c), KS, ex_d(dom_c), block_km = BLOCK_KM)

# --- pre-registered: 5 km survey domain (site membership by nearest grid cell) ------------------------------------------------------------
if ("in_survey5" %in% names(grid)) {
  gsf <- sf::st_as_sf(grid[, c("x", "y", "in_survey5")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
  nn <- sf::st_nearest_feature(spts, gsf)
  keep5 <- grid$in_survey5[nn]
  dom5 <- grid[grid$in_survey5 & !grid$in_core, ]
  add("5 km survey domain", dom5, mk_sets(sc[keep5, ]), KS, ex_d(dom5), block_km = BLOCK_KM)
}

# --- additions NOT in the pre-registration (labelled): fewer nuisance terms -----------------------------------------------------------------
add("ADDITIONAL: without settle", dom, sets, KS, ex_d(dom), drop_cols = c("gm_alluvial", "settle"), block_km = BLOCK_KM)
add("ADDITIONAL: without settle and road", dom, sets, KS, ex_d(dom), drop_cols = c("gm_alluvial", "settle", "road"), block_km = BLOCK_KM)

vp <- do.call(rbind, res)
write_csv_safe(vp, P("output", "tables", "S16_variants_profile.csv"))

# --- pre-registered: thinning to 2 km and 5 km (random greedy, repeated) ----------------------------------------------------------------------
des0 <- design_from_domain(dom); ex6 <- exposure_domain(dom, 6); Xq6 <- cbind(int = 1, des0$dom, z = ex6$z_dom)
thin <- do.call(rbind, lapply(names(sets), function(nm) do.call(rbind, lapply(c(2000, 5000), function(dm) {
  message("  thinning ", nm, " at ", dm / 1000, " km (", THIN_R, " replicates)")
  thin_profile(nm, sets[[nm]], dm, 6, des0, ex6, Xq6, dom$area_km2, THIN_R) }))))
write_csv_safe(thin, P("output", "tables", "S16_thinning_k6.csv"))
write_csv_safe(data.frame(item = c("BLOCK_KM", "THIN_R", "variants", "prior_sd_baseline"),
                          value = c(BLOCK_KM, THIN_R, paste(unique(vp$variant), collapse = " | "), PRIOR_SD)), P("output", "logs", "S16_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S16_sessionInfo.txt"))

cat("\n========== S16 SUMMARY ==========\n")
cat("beta = log intensity ratio per IQR of log distance (positive = intensity RISES with distance); one-step cluster CI, ", BLOCK_KM, " km blocks.\n", sep = "")
show <- c("variant", "subset", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75")
cat("\nk = 6, S1 open-air:\n"); print(vp[vp$k == 6 & vp$subset == "S1_openair", show], row.names = FALSE)
cat("\nk = 6, S2 all sites:\n"); print(vp[vp$k == 6 & vp$subset == "S2_all", show], row.names = FALSE)
cat("\nS1, profile peak and trunk end (k = 5, 6, 7, 8, 9) per variant:\n")
print(reshape(vp[vp$subset == "S1_openair" & vp$k >= 5, c("variant", "k", "beta_iqr")], idvar = "variant", timevar = "k", direction = "wide"), row.names = FALSE)
cat("\nTHINNING (k = 6; random greedy thinning, replicates of beta per IQR):\n"); print(thin, row.names = FALSE)
cs1 <- vp[vp$k == 6 & vp$subset == "S1_openair", ]; cs2 <- vp[vp$k == 6 & vp$subset == "S2_all", ]
s1_ok <- cs1$beta_iqr > 0 & cs1$lo > 0; s2_ok <- cs2$lo <= 0 & cs2$hi >= 0
cat("\nREADING RULE (fixed before the result): S1 k6 positive with CI excluding 0 in", sum(s1_ok, na.rm = TRUE), "of", nrow(cs1), "variants;",
    "S2 k6 CI including 0 in", sum(s2_ok, na.rm = TRUE), "of", nrow(cs2), "variants.\n")
if (any(!s1_ok | is.na(s1_ok))) cat("  S1 sensitive variants:", paste(cs1$variant[!s1_ok | is.na(s1_ok)], collapse = "; "), "\n")
if (any(!s2_ok | is.na(s2_ok))) cat("  S2 sensitive variants:", paste(cs2$variant[!s2_ok | is.na(s2_ok)], collapse = "; "), "\n")
cat("\nNOTE: robustness of the S1 pattern does not address the label dependence found in S13.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S16. Paste this console output back.\n")