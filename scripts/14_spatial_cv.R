# scripts/14_spatial_cv.R
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
source(here::here("R", "spatial_cv.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S14 spatial cross-validation (frozen design): leave-one-block-out predictive log-likelihood ==")

BLOCKS_KM <- as.numeric(strsplit(Sys.getenv("S14_BLOCKS", "40,20"), ",")[[1]])
KS        <- c(3L, 4L, 5L, 6L, 7L, 8L, 9L)

f_man <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
if (!file.exists(f_man)) stop("Design freeze manifest not found. Run scripts/freeze_manifest.R first.", call. = FALSE)
f_q <- P("data", "processed", "quadrature_det.rds"); f_s <- P("data", "processed", "sites_cov.rds")
for (x in c(f_q, f_s)) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid <- sf::st_drop_geometry(readRDS(f_q))
dom  <- grid[grid$in_survey & !grid$in_core, ]
a    <- dom$area_km2
sc   <- readRDS(f_s); sc <- sc[sc$in_domain, ]
des  <- design_from_domain(dom)
ex   <- stats::setNames(lapply(KS, function(k) exposure_domain(dom, k)), KS)
exl  <- exposure_var(dom, "d_lineament")
sets <- list(S1_openair = sc[sc$S1_openair %in% TRUE, ], S2_all = sc, S3_riverbed = sc[sc$S3_riverbed %in% TRUE, ])
message("analysis domain cells: ", nrow(dom), " | sets: ", paste(names(sets), vapply(sets, nrow, integer(1)), sep = "=", collapse = ", "))

summ <- list(); folds_all <- list()
for (nm in names(sets)) for (bk in BLOCKS_KM) {
  sub <- sets[[nm]]
  dsn <- cv_design(dom, sub, des, ex, exl, KS)
  cell_blk <- cluster_ids(as.matrix(dom[, c("x", "y")]), bk * 1000)
  site_blk <- cluster_ids(as.matrix(sub[, c("x", "y")]), bk * 1000)
  message(nm, ", ", bk, " km blocks: ", length(unique(cell_blk)), " folds x ", length(dsn$q), " models")
  fd <- cv_run(dsn$q, dsn$p, a, cell_blk, site_blk)
  fd$set <- nm; fd$block_km <- bk; folds_all[[paste(nm, bk)]] <- fd
  s <- summarise_cv(fd, nrow(sub)); s$set <- nm; s$block_km <- bk; s$n_sites <- nrow(sub)
  s$n_folds_total <- length(unique(cell_blk)); s$n_folds_with_sites <- length(unique(site_blk))
  summ[[paste(nm, bk)]] <- s
}
summ <- do.call(rbind, summ); folds_all <- do.call(rbind, folds_all)
summ <- summ[, c("set", "block_km", "n_sites", "n_folds_total", "n_folds_with_sites", "model", "n_folds_used", "n_folds_dropped",
                 "cv_loglik", "delta_vs_nuisance", "delta_per_site", "folds_improved", "folds_worse", "se_delta_sum")]
write_csv_safe(summ, P("output", "tables", "S14_cv_summary.csv"))
write_csv_safe(folds_all[, c("set", "block_km", "fold", "model", "n_test_sites", "n_test_cells", "loglik", "fit_ok")], P("output", "tables", "S14_cv_folds.csv"))
write_csv_safe(data.frame(item = c("BLOCKS_KM", "models", "prior_sd", "score", "seed"),
                          value = c(paste(BLOCKS_KM, collapse = ","), "intercept-only; nuisance; nuisance + river k3..k9; nuisance + lineament", PRIOR_SD,
                                    "held-out Poisson-process log-likelihood, leave-one-block-out, no leakage", SEED)),
               P("output", "logs", "S14_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S14_sessionInfo.txt"))

cat("\n========== S14 SUMMARY ==========\n")
cat("Score = held-out Poisson-process log-likelihood summed over leave-one-block-out folds (higher is better).\n")
cat("delta_vs_nuisance > 0 means the model predicts held-out blocks better than nuisance terms alone; folds_improved counts blocks where it does.\n")
for (nm in names(sets)) for (bk in BLOCKS_KM) {
  cat("\n", nm, " | ", bk, " km blocks:\n", sep = "")
  print(summ[summ$set == nm & summ$block_km == bk, c("model", "n_folds_used", "n_folds_dropped", "cv_loglik", "delta_vs_nuisance", "delta_per_site", "folds_improved", "folds_worse", "se_delta_sum")], row.names = FALSE)
}
cat("\nNOTE: delta is a sum over folds; se_delta_sum is a rough fold-level scale (sd of fold differences x sqrt(folds)), not a formal test.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S14. Paste this console output back.\n")