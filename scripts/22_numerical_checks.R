# scripts/22_numerical_checks.R
# S22 (EXPLORATORY numerical checks, post-freeze; site-water results HAD been seen):
#  (1) epsilon of the down-weighted Poisson likelihood (frozen 1e-6) at 1e-4, 1e-6, 1e-8, and ridge vs unpenalised ML, S1, k = 3..9;
#  (2) origin of the bootstrap block grid: 40 km blocks with the grid shifted by 0, 10, 20 and 30 km in each axis (16 origins), and the same
#      for 20 km blocks (shifts 0, 5, 10, 15 km); one-step bootstrap; full refit at the worst 40 km origin.
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
source(here::here("R", "numerics_check.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S22 numerical checks (EXPLORATORY, post-freeze) ==")

KS        <- 3:9
K         <- 6L
B_ONE     <- as.integer(Sys.getenv("S22_B_ONESTEP", "2000"))
B_FULL    <- as.integer(Sys.getenv("S22_B_FULL", "1000"))
SKIP_FULL <- Sys.getenv("S22_SKIP_FULL", "0") == "1"
EXPECT_S1 <- as.integer(Sys.getenv("S22_EXPECT_S1", "88"))
EXPECT_B  <- as.numeric(Sys.getenv("S22_EXPECT_BETA", "1.236"))
EXPECT_C40 <- as.integer(Sys.getenv("S22_EXPECT_CL40", "16"))

# reading rules, written BEFORE the run (deviation log A30):
#  (1) epsilon is immaterial if the S1 beta per IQR at k = 3..9 changes by less than 0.005 between eps = 1e-6 and each of 1e-4 and 1e-8;
#  (2) the 40 km result is robust to the block-grid origin if the one-step CI excludes 0 at all 16 origins AND the full-refit CI at the
#      origin with the lowest one-step lower limit also excludes 0.

f <- list(q = P("data", "processed", "quadrature_det.rds"), s = P("data", "processed", "sites_cov.rds"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid   <- sf::st_drop_geometry(readRDS(f$q))
sc_all <- readRDS(f$s)
need_s <- c("x", "y", "S1_openair", "in_domain", "elev_m", "slope_deg", "gm4", "frac_crop", "frac_built", "d_road_major", "d_settle", paste0("d_k", KS))
if (!all(need_s %in% names(sc_all))) stop("sites_cov.rds lacks: ", paste(setdiff(need_s, names(sc_all)), collapse = ", "), call. = FALSE)
dom <- grid[grid$in_survey & !grid$in_core, ]
s1  <- sc_all[(sc_all$S1_openair %in% TRUE) & sc_all$in_domain, ]
if (nrow(s1) != EXPECT_S1) stop("S1 in domain n = ", nrow(s1), " (expected ", EXPECT_S1, ").", call. = FALSE)

des <- design_from_domain(dom); a <- dom$area_km2
ex  <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), as.character(KS))
Xq  <- stats::setNames(lapply(KS, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)), as.character(KS))
Xs1 <- des$apply(s1)
Xp  <- stats::setNames(lapply(KS, function(k) cbind(int = 1, Xs1, z = ex[[as.character(k)]]$apply(s1))), as.character(KS))
iqr_k <- vapply(ex, function(e) e$iqr_per_sd, 0)
fit6  <- fit_real(Xp[["6"]], Xq[["6"]], a)
if (fit6$bad) stop("S1 k = 6 baseline fit failed.", call. = FALSE)
beta6 <- fit6$coef[length(fit6$coef)] * iqr_k[["6"]]
message("domain cells ", nrow(dom), " | S1 ", nrow(s1), " | baseline k = 6 beta per IQR ", round(beta6, 3))
if (abs(beta6 - EXPECT_B) > 0.002) stop("Baseline estimate differs from S11 (", EXPECT_B, ").", call. = FALSE)

# ---- (1) epsilon and ridge vs ML --------------------------------------------------------------------------------------------------
ec <- eps_check(Xp, Xq, a, iqr_k)
ec$beta_iqr <- round(ec$beta_iqr, 6); ec$beta_iqr_ml <- round(ec$beta_iqr_ml, 6); ec$eps_term_share <- signif(ec$eps_term_share, 3)
write_csv_safe(ec, P("output", "tables", "S22_epsilon_check.csv"))
ref <- ec[ec$eps == 1e-6, c("k", "beta_iqr")]; names(ref)[2] <- "beta_ref"
ecm <- merge(ec, ref, by = "k"); ecm$diff_vs_1e6 <- ecm$beta_iqr - ecm$beta_ref
d_up <- max(abs(ecm$diff_vs_1e6[ecm$eps == 1e-4])); d_dn <- max(abs(ecm$diff_vs_1e6[ecm$eps == 1e-8]))
d_ml <- max(abs(ec$beta_iqr[ec$eps == 1e-6] - ec$beta_iqr_ml[ec$eps == 1e-6]))
d_ml6 <- abs(ec$beta_iqr[ec$eps == 1e-6 & ec$k == 6] - ec$beta_iqr_ml[ec$eps == 1e-6 & ec$k == 6])

# ---- (2) block-grid origin ----------------------------------------------------------------------------------------------------------
xy <- as.matrix(s1[, c("x", "y")])
org <- NULL
for (bk in c(40, 20)) {
  offs <- bk * c(0, 0.25, 0.5, 0.75)
  for (ox in offs) for (oy in offs) {
    set.seed(SEED + round(bk * 1000 + ox * 10 + oy))
    org <- rbind(org, origin_row(fit6, Xp[["6"]], xy, bk, ox, oy, ex[["6"]]$iqr_per_sd, ex[["6"]]$sd_log, B_ONE))
  }
}
n_cl0 <- org$n_clusters[org$block_km == 40 & org$ox_km == 0 & org$oy_km == 0]
if (n_cl0 != EXPECT_C40) stop("Origin (0, 0) at 40 km gives ", n_cl0, " clusters (expected ", EXPECT_C40, ").", call. = FALSE)
write_csv_safe(org, P("output", "tables", "S22_block_origin.csv"))
o40 <- org[org$block_km == 40, ]; o20 <- org[org$block_km == 20, ]
worst <- o40[which.min(o40$lo), ]

fw <- NULL
if (!SKIP_FULL) {
  set.seed(SEED + 400L)
  clw <- shift_clusters(xy, 40, worst$ox_km, worst$oy_km)
  bt  <- boot_refit(list(`6` = Xq[["6"]]), list(`6` = Xp[["6"]]), a, clw, B_FULL, progress = 250L)
  fw  <- cbind(origin_ox_km = worst$ox_km, origin_oy_km = worst$oy_km, n_clusters = length(unique(clw)),
               ci_row(bt$beta[, "6"], ex[["6"]]$iqr_per_sd, ex[["6"]]$sd_log, fit6$coef[length(fit6$coef)]))
  write_csv_safe(fw[, c("origin_ox_km", "origin_oy_km", "n_clusters", "beta_iqr", "lo", "hi", "n_boot_ok")], P("output", "tables", "S22_block_origin_worst_fullrefit.csv"))
}
write_csv_safe(data.frame(item = c("KS", "K", "B_ONESTEP", "B_FULL", "PRIOR_SD", "SIM_EPS_frozen", "eps_values", "origins_40km", "origins_20km", "seed", "status"),
                          value = c(paste(KS, collapse = ","), K, B_ONE, if (SKIP_FULL) "skipped" else B_FULL, PRIOR_SD, SIM_EPS, "1e-4,1e-6,1e-8", "0,10,20,30 km per axis", "0,5,10,15 km per axis", SEED,
                                    "EXPLORATORY numerical checks, post-freeze; site-water results had been seen")), P("output", "tables", "S22_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S22_sessionInfo.txt"))

# ---- console summary ----------------------------------------------------------------------------------------------------------------
cat("\n========== S22 SUMMARY (EXPLORATORY numerical checks) ==========\n")
cat("\n(1) S1 beta per IQR by presence weight eps (frozen 1e-6), ridge fit; ML = unpenalised fit at eps 1e-6:\n")
print(reshape(ec[, c("k", "eps", "beta_iqr")], idvar = "k", timevar = "eps", direction = "wide"), row.names = FALSE)
cat("\nmax |beta(1e-4) - beta(1e-6)| over k = 3..9: ", signif(d_up, 3), "\nmax |beta(1e-8) - beta(1e-6)| over k = 3..9: ", signif(d_dn, 3), "\n", sep = "")
cat("max eps-term share of the objective at eps = 1e-4: ", signif(max(ec$eps_term_share[ec$eps == 1e-4]), 3), " | at 1e-6: ", signif(max(ec$eps_term_share[ec$eps == 1e-6]), 3), "\n", sep = "")
cat("ridge vs unpenalised ML at eps 1e-6: k = 6 difference ", signif(d_ml6, 3), " | max over k ", signif(d_ml, 3), "\n", sep = "")
cat("RULE (1): ", if (d_up < 0.005 && d_dn < 0.005) "epsilon immaterial" else "epsilon MATTERS", "\n", sep = "")
cat("\n(2) 40 km blocks, one-step 95% CI, grid origin shifted (ox, oy in km); S1 k = 6, beta per IQR (full sample ", round(beta6, 3), "):\n", sep = "")
print(o40[order(o40$lo), c("ox_km", "oy_km", "n_clusters", "max_cluster_share", "beta_iqr", "lo", "hi", "excludes_zero")], row.names = FALSE)
cat("\n40 km: clusters ", min(o40$n_clusters), " to ", max(o40$n_clusters), " | lower limit ", min(o40$lo), " to ", max(o40$lo), " | upper limit ", min(o40$hi), " to ", max(o40$hi),
    " | CI excludes 0 at ", sum(o40$excludes_zero), " of ", nrow(o40), " origins\n", sep = "")
cat("20 km: clusters ", min(o20$n_clusters), " to ", max(o20$n_clusters), " | lower limit ", min(o20$lo), " to ", max(o20$lo), " | upper limit ", min(o20$hi), " to ", max(o20$hi),
    " | CI excludes 0 at ", sum(o20$excludes_zero), " of ", nrow(o20), " origins\n", sep = "")
if (!is.null(fw)) {
  cat("\nFull refit (", B_FULL, " resamples) at the worst 40 km origin (", worst$ox_km, ", ", worst$oy_km, " km), S11 primary 0.694 to 1.660:\n", sep = "")
  print(fw[, c("origin_ox_km", "origin_oy_km", "n_clusters", "beta_iqr", "lo", "hi", "n_boot_ok")], row.names = FALSE)
  cat("RULE (2): ", if (all(o40$excludes_zero) && fw$lo > 0) "block-origin robust" else "SENSITIVE to the block-grid origin", "\n", sep = "")
}
cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("\nDONE S22. Paste this console output back.\n")