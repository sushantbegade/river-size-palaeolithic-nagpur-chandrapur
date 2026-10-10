# scripts/09b_stepzero_power.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S09b Step 0: power, size and CI coverage with cluster bootstrap (NO real site is used) ==")

# ---- settings (assumptions; any change goes to the deviation log) ----------------------------------
SIM_REPS     <- as.integer(Sys.getenv("S09_REPS", "200"))
SIM_NS       <- c(93L, 48L, 27L)       # S1 open-air; ~MP-only; ~UP-only
SIM_KS_POWER <- c(6L, 8L)              # mid-profile and trunk-river threshold
SIM_RATIOS   <- c(1, 1.25, 1.5, 2, 3)  # intensity ratio per IQR of log distance (1 = null)
SIM_SIGMAS   <- c(0, 0.7)
SIM_RANGE_KM <- 10
BOOT_B       <- 300L
VAL_REPS     <- as.integer(Sys.getenv("S09_VAL_REPS", "5"))
VAL_B        <- 60L

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
cl <- cluster_ids(xy)
message("analysis domain cells: ", nrow(dom), " | 20 km blocks: ", length(unique(cl)))

X  <- nuisance_design_c(dom)
zs <- stats::setNames(lapply(SIZE_K, function(k) z_size(dom, k)), SIZE_K)
message("nuisance columns: ", paste(colnames(X), collapse = ", "))

# ---- A. validate the one-step bootstrap against a full-refit bootstrap ------------------------------
message("A. one-step vs full-refit cluster bootstrap (", VAL_REPS, " reps x 2 settings x B = ", VAL_B, " refits)")
val <- list()
for (set in list(list(n = 93L, sg = 0.7), list(n = 27L, sg = 0.7))) {
  k <- 6L; Xq <- cbind(int = 1, X, z = zs[["6"]]$z); ik <- zs[["6"]]$iqr_per_sd
  for (r in seq_len(VAL_REPS)) {
    eta <- drop(X %*% SIM_TRUTH_C[colnames(X)]) + rff_field(set$sg, SIM_RANGE_KM * 1000)(xy) +
      (log(1.5) / ik) * zs[["6"]]$z
    cnt <- simulate_counts(eta, a, set$n); idx <- rep(seq_along(cnt), cnt)
    if (length(idx) < 8L) next
    Xp <- Xq[idx, , drop = FALSE]; f <- fit_ppm(Xp, Xq, a)
    if (!f$converged) next
    b1 <- boot_onestep(f, Xp, cl[idx], 1000L) * ik
    b2 <- boot_full(Xq, a, Xp, cl[idx], VAL_B) * ik
    val[[length(val) + 1L]] <- data.frame(n_target = set$n, n_sites = length(idx), n_clusters = length(unique(cl[idx])),
                                          se_onestep = sd(b1), se_full = sd(b2), se_model = f$se[ncol(Xq)] * ik,
                                          ratio_onestep_over_full = sd(b1) / sd(b2))
  }
}
val <- do.call(rbind, val)
val[, 4:7] <- round(val[, 4:7], 3)
write_csv_safe(val, P("output", "tables", "S09b_bootstrap_validation.csv"))

# ---- B. power / size / coverage grid -------------------------------------------------------------------
n_fits <- SIM_REPS * length(SIM_SIGMAS) * length(SIM_NS) * length(SIM_KS_POWER) * length(SIM_RATIOS)
message("B. simulation: ", SIM_REPS, " reps x ", length(SIM_SIGMAS), " sigma x ", length(SIM_NS), " n x ",
        length(SIM_KS_POWER), " k x ", length(SIM_RATIOS), " ratios = ", n_fits, " fits (B = ", BOOT_B, " each)")
tt <- Sys.time()
res <- run_power(X, zs, a, xy, SIM_KS_POWER, SIM_NS, SIM_RATIOS, SIM_SIGMAS, SIM_REPS, SIM_RANGE_KM * 1000, B = BOOT_B)
message("simulation minutes: ", round(as.numeric(difftime(Sys.time(), tt, units = "mins")), 2))
ps  <- summarise_power(res)
mde <- mde_table(ps)
gat <- evaluate_gates(ps, mde)
write_csv_safe(ps,  P("output", "tables", "S09b_power_summary.csv"))
write_csv_safe(mde, P("output", "tables", "S09b_mde.csv"))
write_csv_safe(gat, P("output", "tables", "S09b_gates.csv"))
write_csv_safe(data.frame(item = c("SIM_REPS", "SIM_NS", "SIM_KS_POWER", "SIM_RATIOS", "SIM_SIGMAS", "SIM_RANGE_KM", "BOOT_B",
                                   "cluster_block_m", "GATE_TYPE1_MAX", "GATE_COVER_MIN", "GATE_MDE_MAX", "GATE_POWER",
                                   "nuisance_truth", "alluvial", "seed"),
                          value = c(SIM_REPS, paste(SIM_NS, collapse = ","), paste(SIM_KS_POWER, collapse = ","),
                                    paste(SIM_RATIOS, collapse = ","), paste(SIM_SIGMAS, collapse = ","), SIM_RANGE_KM, BOOT_B,
                                    20000, GATE_TYPE1_MAX, GATE_COVER_MIN, GATE_MDE_MAX, GATE_POWER,
                                    paste(paste0(names(SIM_TRUTH_C), "=", SIM_TRUTH_C), collapse = ";"),
                                    "collapsed into pediplain", SEED)),
               P("output", "logs", "S09b_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S09b_sessionInfo.txt"))

cat("\n========== S09b SUMMARY ==========\n")
cat("Analysis domain cells:", nrow(dom), "| nuisance columns:", paste(colnames(X), collapse = ", "), "\n")
cat("\nBootstrap validation (one-step vs full refit; ratio should be ~1):\n"); print(val, row.names = FALSE)
cat("\nPOWER / SIZE / COVERAGE (ratio = 1 is the null; rej_ci = bootstrap CI excludes 0; rej_z = bootstrap-SE z test):\n")
print(as.data.frame(ps[order(ps$sigma, ps$n_target, ps$k, ps$ratio), ]), row.names = FALSE)
cat("\nMinimum detectable ratio per IQR at 80% power (NA = not reached within the grid):\n")
print(as.data.frame(mde[order(mde$sigma, mde$n_target, mde$k), ]), row.names = FALSE)
cat("\nGATES (thresholds: type I <= ", GATE_TYPE1_MAX, ", coverage >= ", GATE_COVER_MIN, ", MDE <= ", GATE_MDE_MAX, " at n = 93):\n", sep = "")
print(gat, row.names = FALSE)
cat("\nNOTE: clusters here are 20 km blocks only; the real analysis uses reference x block. Nuisance truth, field SD, range are assumptions.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S09b. Paste this console output back.\n")