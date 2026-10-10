# scripts/13_site_type_composition.R
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
source(here::here("R", "site_type.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S13 EXPLORATORY (post-hoc, addendum A10): does the site-type label (open-air vs riverbed) depend on distance to large rivers? ==")

BLOCK_KM <- as.numeric(Sys.getenv("S13_BLOCK", "40"))
B        <- as.integer(Sys.getenv("S13_B", "1000"))
KS       <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f_q <- P("data", "processed", "quadrature_det.rds"); f_s <- P("data", "processed", "sites_cov.rds")
for (x in c(f_q, f_s)) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid <- sf::st_drop_geometry(readRDS(f_q))
dom  <- grid[grid$in_survey & !grid$in_core, ]
a    <- dom$area_km2
sc   <- readRDS(f_s); sc <- sc[sc$in_domain, ]
ty   <- sc[(sc$S1_openair %in% TRUE) | (sc$S3_riverbed %in% TRUE), ]
ty$is_riverbed <- as.numeric(ty$S3_riverbed %in% TRUE)          # 1 = riverbed, 0 = open-air (own column: "y" is the site northing)
cl   <- cluster_ids(as.matrix(ty[, c("x", "y")]), BLOCK_KM * 1000)
cl_s1 <- cluster_ids(as.matrix(ty[ty$S1_openair %in% TRUE, c("x", "y")]), BLOCK_KM * 1000)
if (length(unique(cl)) < length(unique(cl_s1))) stop("Cluster ids look wrong: fewer clusters than the S1 subset alone.", call. = FALSE)
message("sites with an open-air or riverbed label in the analysis domain: ", nrow(ty), " (riverbed ", sum(ty$is_riverbed), ", open-air ", sum(1 - ty$is_riverbed),
        ") | clusters at ", BLOCK_KM, " km: ", length(unique(cl)), " | S4 'other' sites excluded: ", sum(sc$S4_other %in% TRUE))

des <- design_from_domain(dom)
Xn  <- des$apply(ty)
ex  <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)

rowfun <- function(k, with_nuisance) {
  z <- ex[[as.character(k)]]
  X <- if (with_nuisance) cbind(int = 1, Xn, z = z$apply(ty)) else cbind(int = 1, z = z$apply(ty))
  f <- fit_logit_ridge(X, ty$is_riverbed)
  if (f$bad) return(data.frame(model = ifelse(with_nuisance, "distance + nuisance", "distance only"), k = k, fit_ok = FALSE))
  bs <- boot_logit(X, ty$is_riverbed, cl, B)
  ci <- ci_row(bs, z$iqr_per_sd, z$sd_log, f$coef[length(f$coef)])
  zq <- stats::quantile(z$z_dom, c(.25, .75), names = FALSE)
  eta0 <- f$coef[1]; bz <- f$coef[length(f$coef)]
  base <- if (with_nuisance) 0 else 0
  pq <- stats::plogis(eta0 + bz * zq + base)
  cbind(data.frame(model = ifelse(with_nuisance, "distance + nuisance", "distance only"), k = k, n = nrow(ty), n_riverbed = sum(ty$is_riverbed)),
        ci, data.frame(p_riverbed_at_q25 = round(pq[1], 3), p_riverbed_at_q75 = round(pq[2], 3), aic = round(f$aic, 2), fit_ok = TRUE))
}
message("fitting ", length(KS), " thresholds x 2 models, each with ", B, " cluster-bootstrap refits")
prof <- do.call(rbind, c(lapply(KS, rowfun, with_nuisance = FALSE), lapply(KS, rowfun, with_nuisance = TRUE)))
prof$delta_aic_descriptive <- NA_real_
for (m in unique(prof$model)) { i <- prof$model == m & prof$fit_ok; prof$delta_aic_descriptive[i] <- round(prof$aic[i] - min(prof$aic[i]), 2) }
write_csv_safe(prof, P("output", "tables", "S13_site_type_profile.csv"))

# consistency check at k = 6: intensity effect in S1 vs S2 (pooled) vs log-odds effect of the type label
Xq6 <- cbind(int = 1, des$dom, z = ex[["6"]]$z_dom)
fit_int <- function(sub) { Xp <- cbind(int = 1, des$apply(sub), z = ex[["6"]]$apply(sub)); f <- fit_real(Xp, Xq6, a); f$coef[length(f$coef)] * ex[["6"]]$iqr_per_sd }
b1 <- fit_int(sc[sc$S1_openair %in% TRUE, ]); b2 <- fit_int(sc); b3 <- fit_int(sc[sc$S3_riverbed %in% TRUE, ])
lg <- prof[prof$k == 6 & prof$model == "distance + nuisance", ]
dec <- data.frame(quantity = c("intensity beta per IQR, S1 open-air", "intensity beta per IQR, S2 all", "intensity beta per IQR, S3 riverbed",
                               "S1 minus S2", "log-odds of riverbed per IQR (type model, with nuisance)", "minus the log-odds value (-beta), for comparison with S1 minus S2"),
                  value = round(c(b1, b2, b3, b1 - b2, lg$beta_iqr, -lg$beta_iqr), 3))
write_csv_safe(dec, P("output", "tables", "S13_decomposition_k6.csv"))
write_csv_safe(data.frame(item = c("BLOCK_KM", "B", "n_sites", "n_riverbed", "n_openair", "n_clusters", "prior_sd", "label"),
                          value = c(BLOCK_KM, B, nrow(ty), sum(ty$is_riverbed), sum(1 - ty$is_riverbed), length(unique(cl)), PRIOR_SD, "EXPLORATORY post-hoc (A10)")),
               P("output", "logs", "S13_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S13_sessionInfo.txt"))

show <- c("model", "k", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "ratio_lo", "ratio_hi", "p_riverbed_at_q25", "p_riverbed_at_q75", "aic", "n_boot_ok")
cat("\n========== S13 SUMMARY (EXPLORATORY) ==========\n")
cat("Outcome: riverbed (1) vs open-air (0) among labelled sites. beta = log-odds of 'riverbed' per IQR increase in log distance;\n")
cat("ratio_q25_vs_q75 = odds of riverbed at the near (q25) distance relative to the far (q75) distance; p_* = predicted P(riverbed) at domain q25 / q75 distance.\n")
cat("\nSites:", nrow(ty), "(riverbed", sum(ty$is_riverbed), "| open-air", sum(1 - ty$is_riverbed), ") | clusters:", length(unique(cl)), "\n")
cat("\nTYPE MODEL PROFILE:\n"); print(prof[prof$fit_ok, show], row.names = FALSE)
cat("\nCONSISTENCY CHECK, k = 6 (approximate: intensity ratio of types is not exactly a log-odds):\n"); print(dec, row.names = FALSE)
cat("\nNOTE: exploratory and post-hoc. A strong dependence of the label on river distance supports the classification-artefact reading of S11-S12;\n",
    "it does not show which explanation (label, burial, behaviour) is true. Label definitions are not independent of position.\n", sep = "")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S13. Paste this console output back.\n")