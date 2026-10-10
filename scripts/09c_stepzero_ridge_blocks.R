# scripts/09c_stepzero_ridge_blocks.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S09c Step 0: ridge-stabilised fit and cluster block size (NO real site is used) ==")

SIM_REPS     <- as.integer(Sys.getenv("S09_REPS", "200"))
SIM_NS       <- c(93L, 48L, 27L)
SIM_KS_POWER <- c(6L, 8L)
SIM_RATIOS   <- c(1, 1.5, 2, 3)
SIM_SIGMAS   <- c(0, 0.7)
BLOCKS_KM    <- c(10, 20, 40)
SIM_RANGE_KM <- 10
BOOT_B       <- 300L

f_q <- P("data", "processed", "quadrature_det.rds")
if (!file.exists(f_q)) stop("Missing: ", f_q, "\n  Run S08 first.", call. = FALSE)
dd   <- sf::st_drop_geometry(readRDS(f_q))
need <- c("x", "y", "in_survey", "in_core", "area_km2", "elev_m", "slope_deg", "gm4",
          "frac_crop", "frac_built", "d_road_major", "d_settle", paste0("d_k", SIZE_K))
if (!all(need %in% names(dd))) stop("quadrature_det.rds lacks: ", paste(setdiff(need, names(dd)), collapse = ", "), call. = FALSE)
dom <- dd[dd$in_survey & !dd$in_core, ]
if (anyNA(dom[, need])) stop("NA in analysis-domain covariates.", call. = FALSE)
a  <- dom$area_km2
xy <- as.matrix(dom[, c("x", "y")])
nblk <- vapply(BLOCKS_KM, function(b) length(unique(cluster_ids(xy, b * 1000))), integer(1))
message("analysis domain cells: ", nrow(dom), " | blocks: ", paste0(BLOCKS_KM, " km = ", nblk, collapse = ", "))

X  <- nuisance_design_c(dom)
zs <- stats::setNames(lapply(SIZE_K, function(k) z_size(dom, k)), SIZE_K)

# ---- A. ridge vs maximum likelihood: shrinkage and stability on easy and hard patterns ----------------
message("A. ridge vs ML on simulated patterns (n = 93 and 27)")
Xq <- cbind(int = 1, X, z = zs[["6"]]$z); ik <- zs[["6"]]$iqr_per_sd
eta0 <- drop(X %*% SIM_TRUTH_C[colnames(X)]) + (log(2) / ik) * zs[["6"]]$z
cmp <- do.call(rbind, lapply(c(93L, 27L), function(nn) {
  do.call(rbind, lapply(seq_len(30L), function(r) {
    cnt <- simulate_counts(eta0, a, nn); idx <- rep(seq_along(cnt), cnt)
    if (length(idx) < 8L) return(NULL)
    Xp <- Xq[idx, , drop = FALSE]
    m <- fit_ppm(Xp, Xq, a); rg <- fit_ppm_ridge(Xp, Xq, a)
    data.frame(n_target = nn, ml_ok = m$converged && all(abs(m$coef) < 8), ridge_ok = rg$converged && all(abs(rg$coef) < 8),
               beta_ml = m$coef[ncol(Xq)] * ik, beta_ridge = rg$coef[ncol(Xq)] * ik)
  }))
}))
cmp_s <- do.call(rbind, lapply(split(cmp, cmp$n_target), function(d) {
  both <- d[d$ml_ok & d$ridge_ok, ]
  data.frame(n_target = d$n_target[1], n_patterns = nrow(d), ml_usable = round(mean(d$ml_ok), 3),
             ridge_usable = round(mean(d$ridge_ok), 3),
             mean_beta_ml = round(mean(both$beta_ml), 3), mean_beta_ridge = round(mean(both$beta_ridge), 3),
             mean_abs_diff = round(mean(abs(both$beta_ml - both$beta_ridge)), 3), true_beta = round(log(2), 3))
}))
write_csv_safe(cmp_s, P("output", "tables", "S09c_ridge_vs_ml.csv"))

# ---- B. grid -------------------------------------------------------------------------------------------
n_fits <- SIM_REPS * length(SIM_SIGMAS) * length(SIM_NS) * length(SIM_KS_POWER) * length(SIM_RATIOS)
message("B. simulation: ", n_fits, " ridge fits, each with ", length(BLOCKS_KM), " bootstraps (B = ", BOOT_B, ")")
tt <- Sys.time()
res <- run_power_c(X, zs, a, xy, SIM_KS_POWER, SIM_NS, SIM_RATIOS, SIM_SIGMAS, BLOCKS_KM, SIM_REPS, SIM_RANGE_KM * 1000, B = BOOT_B)
message("simulation minutes: ", round(as.numeric(difftime(Sys.time(), tt, units = "mins")), 2))
ps  <- summarise_power_c(res)
mde <- mde_table_c(ps)
gat <- evaluate_gates_c(ps, mde)
write_csv_safe(ps,  P("output", "tables", "S09c_power_summary.csv"))
write_csv_safe(mde, P("output", "tables", "S09c_mde.csv"))
write_csv_safe(gat, P("output", "tables", "S09c_gates.csv"))
write_csv_safe(data.frame(item = c("SIM_REPS", "SIM_NS", "SIM_KS_POWER", "SIM_RATIOS", "SIM_SIGMAS", "BLOCKS_KM", "PRIOR_SD", "BOOT_B",
                                   "GATE_TYPE1_MAX", "GATE_COVER_MIN", "GATE_MDE_MAX", "GATE_POWER", "GATE_FAIL_MAX", "seed"),
                          value = c(SIM_REPS, paste(SIM_NS, collapse = ","), paste(SIM_KS_POWER, collapse = ","),
                                    paste(SIM_RATIOS, collapse = ","), paste(SIM_SIGMAS, collapse = ","), paste(BLOCKS_KM, collapse = ","),
                                    PRIOR_SD, BOOT_B, GATE_TYPE1_MAX, GATE_COVER_MIN, GATE_MDE_MAX, GATE_POWER, GATE_FAIL_MAX, SEED)),
               P("output", "logs", "S09c_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S09c_sessionInfo.txt"))

cat("\n========== S09c SUMMARY ==========\n")
cat("Blocks:", paste0(BLOCKS_KM, " km = ", nblk, collapse = ", "), "| prior SD on standardised slopes:", PRIOR_SD, "\n")
cat("\nRidge vs ML (usable = converged, no diverging coefficient):\n"); print(cmp_s, row.names = FALSE)
cat("\nPOWER / SIZE / COVERAGE by block size (rej_uncond counts failed fits as non-rejections; ratio 1 = null):\n")
print(as.data.frame(ps[order(ps$sigma, ps$n_target, ps$k, ps$block_km, ps$ratio), ]), row.names = FALSE)
cat("\nMinimum detectable ratio per IQR (80% power, unconditional):\n")
print(as.data.frame(mde[order(mde$sigma, mde$n_target, mde$k, mde$block_km), ]), row.names = FALSE)
cat("\nGATES (type I <= ", GATE_TYPE1_MAX, ", coverage >= ", GATE_COVER_MIN, ", failed fits <= ", GATE_FAIL_MAX,
    ", MDE <= ", GATE_MDE_MAX, " at n = 93):\n", sep = "")
print(gat, row.names = FALSE)
cat("\nNOTE: simulated clusters are spatial blocks only; the real analysis uses reference x block. All truth values are assumptions.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S09c. Paste this console output back.\n")