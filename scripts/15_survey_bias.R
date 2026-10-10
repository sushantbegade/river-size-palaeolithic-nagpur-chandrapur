# scripts/15_survey_bias.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "survey_bias.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S15 survey-effort sensitivity: reference-aware clusters and a reference-conditional likelihood ==")

B_COND   <- as.integer(Sys.getenv("S15_B", "300"))
KS_BOOT  <- as.integer(strsplit(Sys.getenv("S15_KBOOT", "4,6,8"), ",")[[1]])
KS       <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)
BLOCK_KM <- as.numeric(Sys.getenv("S15_BLOCK", "40"))
BLOCK2   <- as.numeric(Sys.getenv("S15_BLOCK2", "20"))

# reading rule for the conditional analysis, fixed BEFORE the result:
#   same sign and >= half the unconditional S1 k=6 estimate with a reference-bootstrap CI that excludes 0 -> NOT explained by between-reference effort;
#   same sign but CI includes 0 -> inconclusive (low information within footprints); sign flips -> depends on between-reference differences.

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f <- list(q = P("data", "processed", "quadrature_det.rds"), s = P("data", "processed", "sites_cov.rds"), h = P("data", "processed", "hydro_prep.gpkg"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S02, S08 and S10 first.", call. = FALSE)
grid  <- sf::st_drop_geometry(readRDS(f$q))
dom   <- grid[grid$in_survey & !grid$in_core, ]
a     <- dom$area_km2
sc_all <- readRDS(f$s)                       # all 197 sites: footprints use every site of a reference, whatever its type
study <- sf::st_union(sf::st_geometry(sf::st_read(f$h, layer = "boundary", quiet = TRUE)))
sc    <- sc_all[sc_all$in_domain, ]
des   <- design_from_domain(dom)
ex    <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)
Xq_list <- lapply(KS, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)); names(Xq_list) <- KS
sets  <- list(S1_openair = sc[sc$S1_openair %in% TRUE, ], S2_all = sc)
message("domain cells: ", nrow(dom), " | sites: all ", nrow(sc_all), ", in domain ", nrow(sc), " | S1 ", nrow(sets$S1_openair))

# ---- A. unconditional model with reference-aware cluster definitions (one-step bootstrap) ----------------------------
cl_defs <- function(sub) {
  blk  <- cluster_ids(as.matrix(sub[, c("x", "y")]), BLOCK_KM * 1000)
  blk2 <- cluster_ids(as.matrix(sub[, c("x", "y")]), BLOCK2 * 1000)
  list(`block 40 km (primary)` = blk, `reference` = sub$ref_clean,
       `reference x block 40 km` = paste(sub$ref_clean, blk), `reference x block 20 km` = paste(sub$ref_clean, blk2))
}
clu_prof <- do.call(rbind, lapply(names(sets), function(nm) {
  cd <- cl_defs(sets[[nm]])
  do.call(rbind, lapply(names(cd), function(cn) cbind(sample_set = nm, profile_clusters(cn, sets[[nm]], cd[[cn]], KS, des, ex, Xq_list, a))))
}))
write_csv_safe(clu_prof, P("output", "tables", "S15_cluster_definitions.csv"))

# ---- B. reference-conditional likelihood ------------------------------------------------------------------------------------
sites_all_sf <- sf::st_as_sf(sc_all[, c("site_id", "ref_clean", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
fps   <- ref_footprints(sites_all_sf, study)
cells <- footprint_cells(fps, as.matrix(dom[, c("x", "y")]), CRS_ANALYSIS)
fp_sum <- data.frame(reference = names(cells), n_sites_all = as.integer(table(sc_all$ref_clean)[names(cells)]),
                     n_cells = vapply(cells, length, integer(1)))
fp_sum$share_of_domain <- round(fp_sum$n_cells / nrow(dom), 3)
write_csv_safe(fp_sum, P("output", "tables", "S15_footprint_summary.csv"))

cond_one <- function(nm, sub) {
  refs <- sort(unique(sub$ref_clean)); refs <- refs[vapply(refs, function(r) length(cells[[r]]) > 0, logical(1))]
  sub <- sub[sub$ref_clean %in% refs, ]
  Xs_base <- des$apply(sub)
  do.call(rbind, lapply(KS, function(k) {
    z <- ex[[as.character(k)]]
    Xq <- cbind(des$dom, z = z$z_dom)
    Xs <- lapply(refs, function(r) { i <- sub$ref_clean == r; cbind(Xs_base[i, , drop = FALSE], z = z$apply(sub[i, ])) })
    cl <- cells[refs]
    fit <- fit_cond_ref(Xq, a, cl, Xs)
    base <- data.frame(sample_set = nm, k = k, n_sites = nrow(sub), n_references = length(refs))
    if (fit$bad) return(cbind(base, data.frame(beta_iqr = NA, lo = NA, hi = NA, ratio_q25_vs_q75 = NA, n_boot_ok = 0L, converged = fit$converged)))
    est <- fit$coef[length(fit$coef)]
    if (k %in% KS_BOOT) {
      message("  conditional bootstrap: ", nm, ", k = ", k, " (", B_COND, " resamples of ", length(refs), " references)")
      ci <- ci_row(boot_cond_ref(Xq, a, cl, Xs, B_COND), z$iqr_per_sd, z$sd_log, est)[, c("beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "n_boot_ok")]
    } else ci <- data.frame(beta_iqr = round(est * z$iqr_per_sd, 3), lo = NA, hi = NA, ratio_q25_vs_q75 = round(exp(-est * z$iqr_per_sd), 2), n_boot_ok = NA_integer_)
    cbind(base, ci, converged = fit$converged)
  }))
}
cond <- do.call(rbind, lapply(names(sets), function(nm) cond_one(nm, sets[[nm]])))
write_csv_safe(cond, P("output", "tables", "S15_conditional_profile.csv"))

unc6 <- clu_prof[clu_prof$sample_set == "S1_openair" & clu_prof$clusters == "block 40 km (primary)" & clu_prof$k == 6, ]
c1_6 <- cond[cond$sample_set == "S1_openair" & cond$k == 6, ]
verdict <- if (!is.finite(c1_6$beta_iqr)) "conditional fit failed" else if (c1_6$beta_iqr <= 0) "SIGN FLIPS: depends on between-reference differences" else
  if (c1_6$beta_iqr >= 0.5 * unc6$beta_iqr && is.finite(c1_6$lo) && c1_6$lo > 0) "NOT explained by between-reference effort (same sign, >= half, CI excludes 0)" else
    "INCONCLUSIVE (same sign; CI includes 0 or smaller than half)"
write_csv_safe(data.frame(item = c("B_COND", "KS_BOOT", "BLOCK_KM", "unconditional S1 k6 beta_iqr", "conditional S1 k6 beta_iqr", "verdict (rule fixed before result)"),
                          value = c(B_COND, paste(KS_BOOT, collapse = ","), BLOCK_KM, unc6$beta_iqr, c1_6$beta_iqr, verdict)), P("output", "logs", "S15_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S15_sessionInfo.txt"))

cat("\n========== S15 SUMMARY ==========\n")
cat("beta = log intensity ratio per IQR of log distance (positive = intensity RISES with distance).\n")
cat("\nA. UNCONDITIONAL MODEL, k = 6, by cluster definition (one-step bootstrap):\n")
print(clu_prof[clu_prof$k == 6, c("sample_set", "clusters", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75")], row.names = FALSE)
cat("\nB. FOOTPRINTS (survey area per reference, hull of all its sites + 10 km):\n")
cat("   references:", nrow(fp_sum), "| median footprint cells:", stats::median(fp_sum$n_cells), "| share of domain, min/median/max:",
    paste(range(fp_sum$share_of_domain), collapse = " - "), "/", stats::median(fp_sum$share_of_domain), "\n")
cat("\nB. REFERENCE-CONDITIONAL LIKELIHOOD (reference-bootstrap CI only for k in", paste(KS_BOOT, collapse = ","), "):\n")
print(cond, row.names = FALSE)
cat("\nVERDICT for S1 k = 6 (rule fixed before the result):", verdict, "\n")
cat("\nNOTE: the conditional model uses only within-footprint contrasts and removes reference-specific effort; it has much less information.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S15. Paste this console output back.\n")