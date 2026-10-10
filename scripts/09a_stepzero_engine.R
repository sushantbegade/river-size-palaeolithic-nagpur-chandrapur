# scripts/09a_stepzero_engine.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S09a Step 0: engine check and null/alternative simulation (NO real site is used) ==")

# ---- simulation settings (all are assumptions; change only via the deviation log) ------------------
SIM_REPS     <- as.integer(Sys.getenv("S09_REPS", "200"))
SIM_N        <- 93L            # expected number of sites (S1 open-air)
SIM_SIGMAS   <- c(0, 0.7)      # SD of log-Gaussian field: 0 = Poisson, 0.7 = clustered
SIM_RANGE_KM <- 10             # field correlation range (exponential kernel)
SIM_RATIO    <- 1.5            # alternative: intensity ratio per IQR of log distance
SIM_K_TRUE   <- 6L
SIM_KS       <- SIZE_K[SIZE_K >= 3L]   # k = 3,4,...,9 (k = 1 is below coordinate error)

f_q <- P("data", "processed", "quadrature_det.rds")
if (!file.exists(f_q)) stop("Missing: ", f_q, "\n  Run S08 first.", call. = FALSE)
grid <- readRDS(f_q)
dd   <- sf::st_drop_geometry(grid)
need <- c("x", "y", "in_survey", "in_core", "area_km2", "elev_m", "slope_deg", "gm4",
          "frac_crop", "frac_built", "d_road_major", "d_settle", paste0("d_k", SIZE_K))
if (!all(need %in% names(dd))) stop("quadrature_det.rds lacks: ", paste(setdiff(need, names(dd)), collapse = ", "), call. = FALSE)
dom <- dd[dd$in_survey & !dd$in_core, ]
if (anyNA(dom[, need])) stop("NA in analysis-domain covariates.", call. = FALSE)
a   <- dom$area_km2
xy  <- as.matrix(dom[, c("x", "y")])
message("analysis domain cells: ", nrow(dom), " | area km2: ", round(sum(a), 1))

X  <- nuisance_design(dom)
zs <- stats::setNames(lapply(SIZE_K, function(k) z_size(dom, k)), SIZE_K)
cm <- stats::cor(cbind(X, sapply(zs, function(z) z$z)))
nuis_vs_z <- round(cm[colnames(X), as.character(SIZE_K)], 2)
write_csv_safe(cbind(nuisance = rownames(nuis_vs_z), as.data.frame(nuis_vs_z)),
               P("output", "tables", "S09a_nuisance_vs_size_correlations.csv"))
nz <- round(cm[colnames(X), colnames(X)], 2)
write_csv_safe(cbind(nuisance = rownames(nz), as.data.frame(nz)), P("output", "tables", "S09a_nuisance_correlations.csv"))

# ---- A. engine check: Newton fit vs glm lattice and vs down-weighted glm.fit -----------------------
message("A. engine check")
Xq <- cbind(int = 1, X, z = zs[["6"]]$z)
eta <- drop(X %*% SIM_NUIS_TRUTH[colnames(X)]) + 0.4 * zs[["6"]]$z
cnt <- simulate_counts(eta, a, 400)
idx <- rep(seq_along(cnt), cnt)
f   <- fit_ppm(Xq[idx, , drop = FALSE], Xq, a)
gl  <- stats::glm(cnt ~ Xq - 1 + offset(log(a)), family = stats::poisson())
y   <- c(rep(1 / SIM_EPS, length(idx)), rep(0, nrow(Xq)))
w   <- c(rep(SIM_EPS, length(idx)), a)
gf  <- stats::glm.fit(rbind(Xq[idx, , drop = FALSE], Xq), y, weights = w, family = stats::poisson(), start = f$coef)
eng <- data.frame(check = c("newton_converged", "newton_iterations", "max_abs_coef_diff_vs_lattice_glm",
                            "max_abs_coef_diff_vs_dwpr_glm_fit", "max_abs_se_diff_vs_lattice_glm"),
                  value = c(f$converged, f$iter, signif(max(abs(f$coef - stats::coef(gl))), 3),
                            signif(max(abs(f$coef - gf$coefficients)), 3),
                            signif(max(abs(f$se - sqrt(diag(stats::vcov(gl))))), 3)))
write_csv_safe(eng, P("output", "tables", "S09a_engine_check.csv"))

# ---- B. field check ---------------------------------------------------------------------------------
message("B. field check (approximate Gaussian field vs exp(-d/range))")
fc <- field_check(xy, 1, SIM_RANGE_KM * 1000, nfield = 300L, npair = 2000L)
write_csv_safe(fc$by_bin, P("output", "tables", "S09a_field_check.csv"))

# ---- C. null and alternative simulation -------------------------------------------------------------
message("C. simulation: ", SIM_REPS, " replicates x ", length(SIM_SIGMAS), " sigma x 2 scenarios x ", length(SIM_KS), " k")
tt <- Sys.time()
sim <- run_sim(X, zs, a, xy, SIM_KS, SIM_SIGMAS, SIM_RATIO, SIM_K_TRUE, SIM_N, SIM_REPS, SIM_RANGE_KM * 1000)
message("simulation minutes: ", round(as.numeric(difftime(Sys.time(), tt, units = "mins")), 2))
summ <- summarise_sim(sim)
summ$true_beta_iqr <- ifelse(summ$scenario == "alt" & summ$k == SIM_K_TRUE, round(log(SIM_RATIO), 3), 0)
write_csv_safe(summ, P("output", "tables", "S09a_sim_summary.csv"))
conv <- stats::aggregate(cbind(converged, separation) ~ scenario + sigma, sim, mean)
conv$converged <- round(conv$converged, 4); conv$separation <- round(conv$separation, 4)
write_csv_safe(conv, P("output", "tables", "S09a_sim_convergence.csv"))
kh <- khat_table(sim)
write_csv_safe(kh, P("output", "tables", "S09a_khat_preview.csv"))
write_csv_safe(data.frame(item = c("SIM_REPS", "SIM_N", "SIM_SIGMAS", "SIM_RANGE_KM", "SIM_RATIO", "SIM_K_TRUE", "SIM_KS",
                                   "SIM_EPS", "nuisance_truth", "domain_cells", "seed"),
                          value = c(SIM_REPS, SIM_N, paste(SIM_SIGMAS, collapse = ","), SIM_RANGE_KM, SIM_RATIO, SIM_K_TRUE,
                                    paste(SIM_KS, collapse = ","), SIM_EPS,
                                    paste(paste0(names(SIM_NUIS_TRUTH), "=", SIM_NUIS_TRUTH), collapse = ";"), nrow(dom), SEED)),
               P("output", "logs", "S09a_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S09a_sessionInfo.txt"))

cat("\n========== S09a SUMMARY ==========\n")
cat("Analysis domain cells:", nrow(dom), "| nuisance columns:", paste(colnames(X), collapse = ", "), "\n")
cat("\nNuisance vs log-distance (standardised) correlations:\n"); print(nuis_vs_z)
cat("\nLargest |r| among nuisance columns:\n")
nzz <- nz; diag(nzz) <- 0; print(round(max(abs(nzz)), 2)); print(which(abs(nzz) == max(abs(nzz)), arr.ind = TRUE)[1, , drop = FALSE])
cat("\nEngine check:\n"); print(eng, row.names = FALSE)
cat("\nField check (target SD ", fc$target_sd, " | empirical SD ", fc$marginal_sd, "):\n", sep = ""); print(fc$by_bin, row.names = FALSE)
cat("\nConvergence share and share of fits with quasi-separation (a coefficient diverges, usually an empty class):\n"); print(conv, row.names = FALSE)
cat("\nNULL scenario (true effect 0): bias, SE calibration, rejection at 5%:\n")
print(summ[summ$scenario == "null", ], row.names = FALSE)
cat("\nALTERNATIVE (k_true = ", SIM_K_TRUE, ", ratio ", SIM_RATIO, " per IQR; log ratio ", round(log(SIM_RATIO), 3), "):\n", sep = "")
print(summ[summ$scenario == "alt", ], row.names = FALSE)
cat("\nPreview: k chosen by minimum AIC under the alternative (counts of replicates):\n"); print(kh, row.names = FALSE)
cat("\nNOTE: nuisance truth, SIGMA, range and effect size are ASSUMPTIONS. Model-based SEs are NOT the inference method;\n",
    "the real analysis uses the cluster bootstrap. se_ratio < 1 under sigma > 0 shows why.\n", sep = "")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S09a. Paste this console output back.\n")