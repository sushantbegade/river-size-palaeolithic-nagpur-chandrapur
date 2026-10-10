# scripts/11_primary_model.R
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

set.seed(SEED)
t0 <- Sys.time()
message("== S11 primary model: size-threshold profile on S1 open-air sites (pre-registered design) ==")

B         <- as.integer(Sys.getenv("S11_B", "1000"))
BLOCKS_KM <- as.numeric(strsplit(Sys.getenv("S11_BLOCKS", "40,20,10"), ",")[[1]])   # first = primary (pre-registration)
KS        <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)
K_PRIMARY <- 6L

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f_q <- P("data", "processed", "quadrature_det.rds"); f_s <- P("data", "processed", "sites_cov.rds")
for (x in c(f_q, f_s)) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)

grid <- sf::st_drop_geometry(readRDS(f_q))
dom  <- grid[grid$in_survey & !grid$in_core, ]
a    <- dom$area_km2
sc   <- readRDS(f_s)
s1   <- sc[sc$S1_openair %in% TRUE & sc$in_domain, ]
stopifnot(!anyNA(s1[, c("elev_m", "slope_deg", "gm4", "frac_crop", "frac_built", "d_road_major", "d_settle", paste0("d_k", KS))]))
message("analysis domain cells: ", nrow(dom), " | S1 sites in domain: ", nrow(s1))

des <- design_from_domain(dom)
Xs  <- des$apply(s1)
ex  <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)
Xq_list <- lapply(KS, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)); names(Xq_list) <- KS
Xp_list <- lapply(KS, function(k) cbind(int = 1, Xs, z = ex[[as.character(k)]]$apply(s1)));  names(Xp_list) <- KS

# ---- full-sample fits -----------------------------------------------------------------------------------
fits <- lapply(as.character(KS), function(k) fit_real(Xp_list[[k]], Xq_list[[k]], a)); names(fits) <- KS
if (any(vapply(fits, function(f) f$bad, logical(1)))) stop("A full-sample fit failed: ", paste(names(fits)[vapply(fits, function(f) f$bad, logical(1))], collapse = ", "), call. = FALSE)
ml6 <- fit_ppm(Xp_list[["6"]], Xq_list[["6"]], a)
message("k = 6: ridge z coefficient ", round(fits[["6"]]$coef[length(fits[["6"]]$coef)], 4),
        " | ML ", round(ml6$coef[length(ml6$coef)], 4), " (ML converged: ", ml6$converged, ")")

cl_of <- function(bk) cluster_ids(as.matrix(s1[, c("x", "y")]), bk * 1000)
prim_bk <- BLOCKS_KM[1]

# ---- bootstrap: all k at the primary block size; k = 6 only at the sensitivity block sizes ---------------------
message("bootstrap, primary block ", prim_bk, " km, all k: ", B, " resamples x ", length(KS), " refits")
bp <- boot_refit(Xq_list, Xp_list, a, cl_of(prim_bk), B)
boots <- list(); boots[[as.character(prim_bk)]] <- bp$beta[, as.character(K_PRIMARY)]
for (bk in BLOCKS_KM[-1]) {
  message("bootstrap, sensitivity block ", bk, " km, k = ", K_PRIMARY)
  boots[[as.character(bk)]] <- boot_refit(Xq_list[as.character(K_PRIMARY)], Xp_list[as.character(K_PRIMARY)], a, cl_of(bk), B)$beta[, 1]
}

# ---- primary result (k = 6) ------------------------------------------------------------------------------------
est6 <- fits[["6"]]$coef[length(fits[["6"]]$coef)]
prim <- do.call(rbind, lapply(names(boots), function(bk) {
  r <- ci_row(boots[[bk]], ex[["6"]]$iqr_per_sd, ex[["6"]]$sd_log, est6)
  cbind(block_km = as.numeric(bk), method = "full-refit bootstrap", r)
}))
for (bk in BLOCKS_KM) {                                   # cross-check: one-step bootstrap (frozen helper)
  Xp <- Xp_list[["6"]]; o <- boot_onestep(fits[["6"]], Xp, cl_of(bk), 2000L)
  prim <- rbind(prim, cbind(block_km = bk, method = "one-step bootstrap (check)", ci_row(o, ex[["6"]]$iqr_per_sd, ex[["6"]]$sd_log, est6)))
}
prim$se_model_beta_iqr <- round(fits[["6"]]$se[length(fits[["6"]]$se)] * ex[["6"]]$iqr_per_sd, 3)
prim$n_sites <- nrow(s1); prim$n_clusters <- vapply(prim$block_km, function(b) length(unique(cl_of(b))), integer(1))
write_csv_safe(prim, P("output", "tables", "S11_primary_k6.csv"))

fb <- prim[prim$method == "full-refit bootstrap", ]
supported <- fb$hi[fb$block_km == prim_bk] < 0
decision <- data.frame(item = c("confirmatory k", "sample", "n sites in analysis domain", "primary block (km)", "beta per IQR",
                                "95% CI low", "95% CI high", "intensity ratio q25 vs q75 (>1 = more sites near drainage)", "ratio 95% CI",
                                "H1 supported (pre-registered rule: 40 km CI excludes 0, expected sign)"),
                       value = c(K_PRIMARY, "S1 open-air", nrow(s1), prim_bk, fb$beta_iqr[fb$block_km == prim_bk],
                                 fb$lo[fb$block_km == prim_bk], fb$hi[fb$block_km == prim_bk],
                                 fb$ratio_q25_vs_q75[fb$block_km == prim_bk],
                                 paste0(fb$ratio_lo[fb$block_km == prim_bk], " to ", fb$ratio_hi[fb$block_km == prim_bk]),
                                 ifelse(supported, "YES", "NO")))
write_csv_safe(decision, P("output", "tables", "S11_decision.csv"))

# ---- profile (all k) at the primary block --------------------------------------------------------------------------
prof <- do.call(rbind, lapply(as.character(KS), function(k) {
  est <- fits[[k]]$coef[length(fits[[k]]$coef)]
  cbind(k = as.integer(k), ci_row(bp$beta[, k], ex[[k]]$iqr_per_sd, ex[[k]]$sd_log, est),
        aic = round(fits[[k]]$aic, 2), converged = fits[[k]]$converged)
}))
prof$delta_aic_descriptive <- round(prof$aic - min(prof$aic), 2)
write_csv_safe(prof, P("output", "tables", "S11_profile.csv"))
okrow <- stats::complete.cases(bp$aic)
khat_b <- as.integer(colnames(bp$aic))[apply(bp$aic[okrow, , drop = FALSE], 1, which.min)]
kd <- data.frame(k = KS, share_of_bootstrap_best_k = round(vapply(KS, function(k) mean(khat_b == k), numeric(1)), 3))
kd$full_sample_best <- kd$k == prof$k[which.min(prof$aic)]
write_csv_safe(kd, P("output", "tables", "S11_khat_distribution.csv"))
kbest <- prof$k[which.min(prof$aic)]
within1 <- sum(kd$share_of_bootstrap_best_k[abs(kd$k - kbest) <= 1])

# ---- nuisance coefficients at k = 6 (model-based SEs, descriptive) -----------------------------------------------
cf <- fits[["6"]]
coefs <- data.frame(term = c("intercept", colnames(des$dom), "log distance k6 (z)"), estimate_std = round(cf$coef, 3),
                    se_model = round(cf$se, 3))
write_csv_safe(coefs, P("output", "tables", "S11_coefficients_k6.csv"))
write_csv_safe(data.frame(item = c("B", "BLOCKS_KM", "PRIOR_SD", "K_PRIMARY", "n_sites", "n_domain_cells", "seed", "bootstrap_failures_primary_block"),
                          value = c(B, paste(BLOCKS_KM, collapse = ","), PRIOR_SD, K_PRIMARY, nrow(s1), nrow(dom), SEED,
                                    sum(!stats::complete.cases(bp$aic)))),
               P("output", "logs", "S11_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S11_sessionInfo.txt"))

cat("\n========== S11 SUMMARY ==========\n")
cat("S1 sites in analysis domain:", nrow(s1), "| domain cells:", nrow(dom), "| clusters at ", prim_bk, " km: ",
    length(unique(cl_of(prim_bk))), "\n", sep = "")
cat("\nPRIMARY (k = 6). beta = log intensity ratio per IQR of log distance; ratio = intensity at q25 distance / at q75 distance:\n")
print(prim[, c("block_km", "method", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "ratio_lo", "ratio_hi", "per_doubling", "n_boot_ok", "se_model_beta_iqr")], row.names = FALSE)
cat("\nPre-registered rule for H1: 40 km full-refit CI excludes 0 in the expected direction (beta < 0) ->", ifelse(supported, "SUPPORTED", "NOT SUPPORTED"), "\n")
cat("(Absence of support is not evidence of absence; minimum detectable ratio is about 1.9 per IQR.)\n")
cat("\nPROFILE (k = 3..9; ", prim_bk, " km blocks; delta-AIC is descriptive only, no thresholds):\n", sep = "")
print(prof[, c("k", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "ratio_lo", "ratio_hi", "per_doubling", "aic", "delta_aic_descriptive", "n_boot_ok")], row.names = FALSE)
cat("\nBootstrap distribution of the best k (min AIC per resample); full-sample best k =", kbest,
    "| bootstrap mass within +-1 of it:", round(within1, 3), "\n")
print(kd, row.names = FALSE)
cat("\nCoefficients at k = 6 (standardised; model-based SE; descriptive):\n"); print(coefs, row.names = FALSE)
cat("\nNOTE: trunk-river thresholds (k8-k9) are exploratory. Negative control, leave-one-reference-out, spatial CV and robustness follow in S12-S14.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S11. Paste this console output back.\n")