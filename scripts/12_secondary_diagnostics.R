# scripts/12_secondary_diagnostics.R
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

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S12 secondary diagnostics (frozen design): sample sets, coordinate tiers, negative control, leave-one-out ==")

BLOCK_KM <- as.numeric(Sys.getenv("S12_BLOCK", "40"))
B_ONE    <- as.integer(Sys.getenv("S12_B_ONE", "2000"))
B_FULL   <- as.integer(Sys.getenv("S12_B_FULL", "1000"))
KS       <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f_q <- P("data", "processed", "quadrature_det.rds"); f_s <- P("data", "processed", "sites_cov.rds")
for (x in c(f_q, f_s)) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid <- sf::st_drop_geometry(readRDS(f_q))
dom  <- grid[grid$in_survey & !grid$in_core, ]
a    <- dom$area_km2
sc   <- readRDS(f_s)
sc   <- sc[sc$in_domain, ]
des  <- design_from_domain(dom)
ex   <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)
Xq_list <- lapply(KS, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)); names(Xq_list) <- KS
s1 <- sc[sc$S1_openair %in% TRUE, ]
message("analysis domain cells: ", nrow(dom), " | sites in domain: ", nrow(sc), " | S1: ", nrow(s1), " | block: ", BLOCK_KM, " km")

# ---- 1. sample sets (S1 = primary) ---------------------------------------------------------------------------
sets <- list(S1_openair = s1, S2_all = sc, S3_riverbed = sc[sc$S3_riverbed %in% TRUE, ], S4_other = sc[sc$S4_other %in% TRUE, ])
prof_sets <- do.call(rbind, lapply(names(sets), function(nm) subset_profile(nm, sets[[nm]], KS, des, ex, Xq_list, a, BLOCK_KM, B_ONE)))
write_csv_safe(prof_sets, P("output", "tables", "S12_sample_sets_profile.csv"))

# ---- 2. coordinate tiers within S1 -----------------------------------------------------------------------------
tier_def <- list(S1_tier_A = s1[s1$coord_tier == "A_field_verified", ],
                 S1_tier_AB = s1[s1$coord_tier %in% c("A_field_verified", "B_author_reported"), ],
                 S1_tier_B = s1[s1$coord_tier == "B_author_reported", ],
                 S1_tier_C = s1[s1$coord_tier == "C_map_based", ])
prof_tier <- do.call(rbind, lapply(names(tier_def), function(nm) subset_profile(nm, tier_def[[nm]], KS, des, ex, Xq_list, a, BLOCK_KM, B_ONE)))
write_csv_safe(prof_tier, P("output", "tables", "S12_tiers_profile.csv"))

# ---- 3. negative control: lineament distance replaces river distance (S1) ---------------------------------------------
exl <- exposure_var(dom, "d_lineament")
Xq_lin <- cbind(int = 1, des$dom, z = exl$z_dom)
Xp_lin <- cbind(int = 1, des$apply(s1), z = exl$apply(s1))
f_lin <- fit_real(Xp_lin, Xq_lin, a)
cl1 <- cluster_ids(as.matrix(s1[, c("x", "y")]), BLOCK_KM * 1000)
message("negative control: lineament, ", B_FULL, " full-refit bootstrap resamples")
bl <- boot_refit(list(lin = Xq_lin), list(lin = Xp_lin), a, cl1, B_FULL)$beta[, 1]
nc_row <- cbind(model = "lineament only (negative control)", ci_row(bl, exl$iqr_per_sd, exl$sd_log, f_lin$coef[length(f_lin$coef)]))
# joint model: river k6 + lineament
Xq_j <- cbind(Xq_list[["6"]], lin = exl$z_dom); Xp_j <- cbind(int = 1, des$apply(s1), z = ex[["6"]]$apply(s1), lin = exl$apply(s1))
f_j <- fit_real(Xp_j, Xq_j, a); pj <- ncol(Xq_j)
bj6 <- boot_onestep_j(f_j, Xp_j, cl1, pj - 1L, B_ONE); bjl <- boot_onestep_j(f_j, Xp_j, cl1, pj, B_ONE)
nc_row <- rbind(nc_row,
                cbind(model = "joint: river k6 coefficient", ci_row(bj6, ex[["6"]]$iqr_per_sd, ex[["6"]]$sd_log, f_j$coef[pj - 1L])),
                cbind(model = "joint: lineament coefficient", ci_row(bjl, exl$iqr_per_sd, exl$sd_log, f_j$coef[pj])))
nc_row$aic <- c(round(f_lin$aic, 2), round(f_j$aic, 2), round(f_j$aic, 2))
write_csv_safe(nc_row, P("output", "tables", "S12_negative_control.csv"))

# ---- 4. leave-one-reference-out and leave-one-block-out (k = 6, S1) --------------------------------------------------------
s1$blk <- cl1
loro <- leave_one_out(s1, "ref_clean", des, ex[["6"]], Xq_list[["6"]], a)
lobo <- leave_one_out(s1, "blk", des, ex[["6"]], Xq_list[["6"]], a)
write_csv_safe(loro, P("output", "tables", "S12_leave_one_reference_out.csv"))
write_csv_safe(lobo, P("output", "tables", "S12_leave_one_block_out.csv"))
write_csv_safe(data.frame(item = c("BLOCK_KM", "B_ONE", "B_FULL", "k_loo", "n_S1", "seed"),
                          value = c(BLOCK_KM, B_ONE, B_FULL, 6, nrow(s1), SEED)), P("output", "logs", "S12_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S12_sessionInfo.txt"))

show <- c("subset", "n_sites", "n_clusters", "k", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "ratio_lo", "ratio_hi", "aic")
cat("\n========== S12 SUMMARY ==========\n")
cat("beta = log intensity ratio per IQR of log distance (positive = intensity RISES with distance); ratio = intensity at q25 / at q75 distance.\n")
cat("CIs here are one-step cluster-bootstrap approximations (", BLOCK_KM, " km blocks), not full refits.\n", sep = "")
cat("\n1. SAMPLE SETS, k = 6:\n"); print(prof_sets[prof_sets$k == 6, show], row.names = FALSE)
cat("\n   Sample sets, full profile (beta per IQR with CI):\n"); print(prof_sets[, c("subset", "k", "n_sites", "beta_iqr", "lo", "hi", "aic")], row.names = FALSE)
cat("\n2. COORDINATE TIERS within S1, k = 6:\n"); print(prof_tier[prof_tier$k == 6, show], row.names = FALSE)
cat("\n   Tiers, full profile:\n"); print(prof_tier[, c("subset", "k", "n_sites", "beta_iqr", "lo", "hi", "aic")], row.names = FALSE)
cat("\n3. NEGATIVE CONTROL (lineament distance), S1:\n"); print(nc_row[, c("model", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "ratio_lo", "ratio_hi", "n_boot_ok", "aic")], row.names = FALSE)
cat("\n4. LEAVE-ONE-REFERENCE-OUT, k = 6 (largest influence first; full-sample beta shown in S11):\n"); print(utils::head(loro, 10), row.names = FALSE)
cat("   range of beta over", nrow(loro), "references:", min(loro$beta_iqr, na.rm = TRUE), "to", max(loro$beta_iqr, na.rm = TRUE), "| failed fits:", sum(!loro$fit_ok), "\n")
cat("\n   LEAVE-ONE-BLOCK-OUT (", BLOCK_KM, " km), k = 6:\n", sep = ""); print(lobo, row.names = FALSE)
cat("   range of beta over", nrow(lobo), "blocks:", min(lobo$beta_iqr, na.rm = TRUE), "to", max(lobo$beta_iqr, na.rm = TRUE), "\n")
cat("\nNOTE: all items are secondary/robustness analyses of the frozen design. They help interpret the S11 result; they do not redefine H1.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S12. Paste this console output back.\n")