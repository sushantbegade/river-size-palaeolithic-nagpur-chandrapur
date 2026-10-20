# scripts/17_exploratory.R
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
source(here::here("R", "exploratory.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S17 EXPLORATORY items (frozen-design document, section 5): distance shape, single-period models, period contrast ==")

BLOCK_KM <- as.numeric(Sys.getenv("S17_BLOCK", "40"))
B_ONE    <- as.integer(Sys.getenv("S17_B_ONE", "2000"))
B_PER    <- as.integer(Sys.getenv("S17_B_PERIOD", "300"))
KS_SHAPE <- c(4L, 6L, 8L)
KS_PER   <- c(5L, 6L, 7L)

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f_q <- P("data", "processed", "quadrature_det.rds"); f_s <- P("data", "processed", "sites_cov.rds")
for (x in c(f_q, f_s)) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid <- sf::st_drop_geometry(readRDS(f_q)); dom <- grid[grid$in_survey & !grid$in_core, ]; a <- dom$area_km2
sc <- readRDS(f_s); sc <- sc[sc$in_domain, ]
des <- design_from_domain(dom)
KS <- sort(unique(c(KS_SHAPE, KS_PER)))
ex <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)
Xq_list <- lapply(KS, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)); names(Xq_list) <- KS
s1 <- sc[sc$S1_openair %in% TRUE, ]; s3 <- sc[sc$S3_riverbed %in% TRUE, ]

# ---- A. distance shape: five equal-area distance classes (class 3 = reference) --------------------------------------------------------
shape_sets <- list(S1_openair = s1, S2_all = sc, S3_riverbed = s3)
shape <- do.call(rbind, lapply(names(shape_sets), function(nm)
  do.call(rbind, lapply(KS_SHAPE, function(k) binned_fit(nm, k, shape_sets[[nm]], ex[[as.character(k)]], des, dom, a, BLOCK_KM, B_ONE)))))
write_csv_safe(shape, P("output", "tables", "S17_distance_shape.csv"))

# ---- B. single-period models and period contrast --------------------------------------------------------------------------------------
per_def <- function(d) list(LP = d[d$period_class == "LP", ], MP = d[d$period_class == "MP", ], UP = d[d$period_class == "UP", ])
per_counts <- do.call(rbind, lapply(c("S1_openair", "S2_all"), function(nm) {
  d <- if (nm == "S1_openair") s1 else sc
  data.frame(sample_set = nm, period = c("LP", "MP", "UP", "multi"), n_sites = c(sum(d$period_class == "LP"), sum(d$period_class == "MP"),
                                                                                 sum(d$period_class == "UP"), sum(d$period_class == "multi")))
}))
write_csv_safe(per_counts, P("output", "tables", "S17_period_counts.csv"))
per_prof <- do.call(rbind, lapply(c("S1_openair", "S2_all"), function(nm) {
  d <- if (nm == "S1_openair") s1 else sc
  do.call(rbind, lapply(names(per_def(d)), function(p) {
    r <- subset_profile(paste(nm, p, sep = "_"), per_def(d)[[p]], KS_PER, des, ex, Xq_list, a, BLOCK_KM, B_ONE); cbind(sample_set = nm, period = p, r)
  }))
}))
write_csv_safe(per_prof, P("output", "tables", "S17_period_profiles.csv"))

contrast <- do.call(rbind, lapply(c("S1_openair", "S2_all"), function(nm) {
  d <- if (nm == "S1_openair") s1 else sc; subs <- per_def(d)
  message("period contrast bootstrap, ", nm, " (", B_PER, " resamples x 3 periods, k = 6)")
  bm <- boot_period_contrast(subs, des, ex[["6"]], Xq_list[["6"]], a, BLOCK_KM, B_PER)
  pr <- list(c("LP", "MP"), c("MP", "UP"), c("LP", "UP"))
  do.call(rbind, lapply(pr, function(p) {
    Xp <- function(s) cbind(int = 1, des$apply(s), z = ex[["6"]]$apply(s))
    est <- function(s) { f <- fit_real(Xp(s), Xq_list[["6"]], a); if (f$bad) NA_real_ else f$coef[length(f$coef)] * ex[["6"]]$iqr_per_sd }
    e1 <- est(subs[[p[1]]]); e2 <- est(subs[[p[2]]]); dd <- bm[, p[1]] - bm[, p[2]]; dd <- dd[is.finite(dd)]
    q <- stats::quantile(dd, c(.025, .975), names = FALSE)
    data.frame(sample_set = nm, contrast = paste(p, collapse = " - "), beta_first = round(e1, 3), beta_second = round(e2, 3),
               difference = round(e1 - e2, 3), lo = round(q[1], 3), hi = round(q[2], 3), n_boot_ok = length(dd))
  }))
}))
write_csv_safe(contrast, P("output", "tables", "S17_period_contrast.csv"))
write_csv_safe(data.frame(item = c("BLOCK_KM", "B_ONE", "B_PERIOD", "KS_SHAPE", "KS_PER", "label"),
                          value = c(BLOCK_KM, B_ONE, B_PER, paste(KS_SHAPE, collapse = ","), paste(KS_PER, collapse = ","), "EXPLORATORY (declared exploratory in the frozen design)")),
               P("output", "logs", "S17_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S17_sessionInfo.txt"))

cat("\n========== S17 SUMMARY (EXPLORATORY) ==========\n")
cat("A. DISTANCE SHAPE: five equal-area classes of distance to size-k rivers (class 1 nearest). rel_intensity = intensity relative to class 3, adjusted for nuisance terms; CI = one-step cluster bootstrap.\n")
print(shape[, c("subset", "k", "bin", "km_from", "km_to", "n_sites", "sites_per_100km2", "rel_intensity", "lo", "hi")], row.names = FALSE)
cat("\nB. PERIOD COUNTS (in analysis domain):\n"); print(per_counts, row.names = FALSE)
cat("\n   Single-period profiles (beta per IQR; UP samples are descriptive only):\n")
print(per_prof[per_prof$k == 6, c("sample_set", "period", "n_sites", "n_clusters", "k", "beta_iqr", "lo", "hi")], row.names = FALSE)
cat("\n   Period contrasts at k = 6 (difference of beta per IQR; block bootstrap):\n"); print(contrast, row.names = FALSE)
cat("\nNOTE: exploratory. Single-period samples are small (minimum detectable ratio >= 2.75 for n = 27); do not read non-significance as absence.\n")
cat("Not run: visibility and preservation offsets (no usable independent covariate), trunk-only (covered by the k8-k9 profile in S11-S16).\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S17. Paste this console output back.\n")