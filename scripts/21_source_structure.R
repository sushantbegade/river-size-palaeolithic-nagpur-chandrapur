# scripts/21_source_structure.R
# S21 (EXPLORATORY, post-freeze; site-water results HAD been seen). Source structure of (A) the S1 k = 6 result and (B) the S1/S3 labels.
#  A1 author-family counts; A2 S1 k = 6 full-refit bootstrap with clusters = author family and family x 40 km block;
#  A3 leave-one-family-out and a full-precision leave-one-block-out (check of the 1.118 tie in S12);
#  B  within-reference test: among references that report BOTH labels, does the riverbed label depend on distance within the reference?
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
if (!requireNamespace("survival", quietly = TRUE)) stop("Package survival is needed (it ships with R as a recommended package).", call. = FALSE)
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "model_primary.R"))
source(here::here("R", "diagnostics.R"))
source(here::here("R", "domain_sens.R"))
source(here::here("R", "source_struct.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S21 source structure (EXPLORATORY, post-freeze) ==")

K        <- 6L
BLOCK_KM <- 40
B_FULL   <- as.integer(Sys.getenv("S21_B_FULL", "1000"))
B_LABEL  <- as.integer(Sys.getenv("S21_B_LABEL", "1000"))
SKIP_FULL <- Sys.getenv("S21_SKIP_FULL", "0") == "1"
EXPECT_S1 <- as.integer(Sys.getenv("S21_EXPECT_S1", "88"))
EXPECT_B  <- as.numeric(Sys.getenv("S21_EXPECT_BETA", "1.236"))   # S11 baseline estimate; a mismatch means the data or code changed

# reading rule for B, written BEFORE the run (deviation log A28): the riverbed label depends on distance WITHIN references if the
# fixed-effect estimate (references reporting both labels) is negative with a 95 % CI that excludes 0; the label is mainly a source convention
# if that CI includes 0. If it excludes 0, its size is compared with the pooled estimate on all S1 + S3 sites; no causal reading either way.

f <- list(q = P("data", "processed", "quadrature_det.rds"), s = P("data", "processed", "sites_cov.rds"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S08 and S10 first.", call. = FALSE)
grid   <- sf::st_drop_geometry(readRDS(f$q))
sc_all <- readRDS(f$s)
need_s <- c("site_id", "ref_id", "ref_clean", "x", "y", "S1_openair", "S3_riverbed", "in_domain", "elev_m", "slope_deg", "gm4", "frac_crop", "frac_built",
            "d_road_major", "d_settle", paste0("d_k", K))
if (!all(need_s %in% names(sc_all))) stop("sites_cov.rds lacks: ", paste(setdiff(need_s, names(sc_all)), collapse = ", "), call. = FALSE)
dom <- grid[grid$in_survey & !grid$in_core, ]
sc  <- sc_all[sc_all$in_domain, ]
sc$fam <- ref_family(sc$ref_clean)
s1  <- sc[sc$S1_openair %in% TRUE, ]
if (nrow(s1) != EXPECT_S1) stop("S1 in domain n = ", nrow(s1), " (expected ", EXPECT_S1, ").", call. = FALSE)
message("domain cells: ", nrow(dom), " | S1 ", nrow(s1), " | S3 ", sum(sc$S3_riverbed %in% TRUE))

des <- design_from_domain(dom); ex <- exposure_domain(dom, K); a <- dom$area_km2
Xq  <- cbind(int = 1, des$dom, z = ex$z_dom)
Xp1 <- cbind(int = 1, des$apply(s1), z = ex$apply(s1))
full <- fit_real(Xp1, Xq, a)
if (full$bad) stop("S1 k = 6 baseline fit failed.", call. = FALSE)
beta_full <- full$coef[length(full$coef)] * ex$iqr_per_sd
message("baseline S1 k = 6 beta per IQR: ", round(beta_full, 3), " (S11: ", EXPECT_B, ")")
if (abs(beta_full - EXPECT_B) > 0.002) stop("Baseline estimate differs from S11 (", EXPECT_B, ").", call. = FALSE)

# ---- A1: author-family structure -------------------------------------------------------------------------------------------------
ft <- family_table(sc, sc$fam)
write_csv_safe(ft, P("output", "tables", "S21_family_table.csv"))
fs <- ft[ft$n_S1 > 0, ]
fam_sum <- data.frame(item = c("references (strings) with S1 sites", "author families with S1 sites", "top-1 family share of S1", "top-3 family share of S1", "families with exactly 1 S1 site",
                               "references with both labels (S1 and S3)", "references with S1 only", "references with S3 only"),
                      value = c(length(unique(s1$ref_clean)), nrow(fs), round(max(fs$share_of_S1), 3), round(sum(sort(fs$share_of_S1, decreasing = TRUE)[1:min(3, nrow(fs))]), 3), sum(fs$n_S1 == 1),
                                NA, NA, NA))
rt <- table(sc$ref_clean[sc$S1_openair %in% TRUE | sc$S3_riverbed %in% TRUE], ifelse(sc$S3_riverbed[sc$S1_openair %in% TRUE | sc$S3_riverbed %in% TRUE] %in% TRUE, "S3", "S1"))
both <- rownames(rt)[rt[, "S1"] > 0 & rt[, "S3"] > 0]
fam_sum$value[6] <- length(both); fam_sum$value[7] <- sum(rt[, "S1"] > 0 & rt[, "S3"] == 0); fam_sum$value[8] <- sum(rt[, "S3"] > 0 & rt[, "S1"] == 0)
write_csv_safe(fam_sum, P("output", "tables", "S21_family_summary.csv"))

# ---- A2: S1 k = 6 bootstrap with family clusters --------------------------------------------------------------------------------
blk <- cluster_ids(as.matrix(s1[, c("x", "y")]), BLOCK_KM * 1000)
a2 <- NULL
if (!SKIP_FULL) {
  one <- function(label, cl) {
    set.seed(SEED + nchar(label))
    bt <- boot_refit(list(`6` = Xq), list(`6` = Xp1), a, cl, B_FULL, progress = 250L)
    cbind(clusters = label, n_clusters = length(unique(cl)), ci_row(bt$beta[, "6"], ex$iqr_per_sd, ex$sd_log, full$coef[length(full$coef)]))
  }
  a2 <- rbind(one("author family", s1$fam), one("author family x 40 km block", paste(s1$fam, blk)))
  write_csv_safe(a2, P("output", "tables", "S21_family_bootstrap_k6.csv"))
}

# ---- A3: leave-one-family-out and full-precision leave-one-block-out ---------------------------------------------------------------
s1$blk <- blk
lof <- loo_full(s1, "fam", des, ex, Xq, a)
lob <- loo_full(s1, "blk", des, ex, Xq, a)
write_csv_safe(lof, P("output", "tables", "S21_leave_one_family_out_k6.csv"))
write_csv_safe(lob, P("output", "tables", "S21_leave_one_block_out_k6_fullprecision.csv"))

# ---- B: within-reference label test ----------------------------------------------------------------------------------------------
sb <- sc[(sc$S1_openair %in% TRUE) | (sc$S3_riverbed %in% TRUE), ]
yb <- as.numeric(sb$S3_riverbed %in% TRUE)
zb <- ex$apply(sb); Xn <- des$apply(sb)
mixed <- sb$ref_clean %in% both
fam_mixed <- sb$fam %in% names(which(tapply(yb, sb$fam, function(v) length(unique(v)) > 1)))
set.seed(SEED + 500L)
rows <- list(
  label_row("all S1+S3 sites, pooled (no reference effect)", yb, zb, sb$ref_clean, sb$ref_clean, ex$iqr_per_sd, fe = FALSE, B = B_LABEL),
  label_row("references with both labels, pooled", yb[mixed], zb[mixed], sb$ref_clean[mixed], sb$ref_clean[mixed], ex$iqr_per_sd, fe = FALSE, B = B_LABEL),
  label_row("references with both labels, fixed effect per reference", yb[mixed], zb[mixed], sb$ref_clean[mixed], sb$ref_clean[mixed], ex$iqr_per_sd, fe = TRUE, B = B_LABEL),
  label_row("author families with both labels, fixed effect per family", yb[fam_mixed], zb[fam_mixed], sb$fam[fam_mixed], sb$fam[fam_mixed], ex$iqr_per_sd, fe = TRUE, B = B_LABEL),
  label_row("references with both labels, pooled + nuisance", yb[mixed], zb[mixed], sb$ref_clean[mixed], sb$ref_clean[mixed], ex$iqr_per_sd, X = Xn[mixed, , drop = FALSE], fe = FALSE, B = B_LABEL),
  label_row("references with both labels, fixed effect + nuisance", yb[mixed], zb[mixed], sb$ref_clean[mixed], sb$ref_clean[mixed], ex$iqr_per_sd, X = Xn[mixed, , drop = FALSE], fe = TRUE, B = B_LABEL))
lab <- do.call(rbind, rows)
write_csv_safe(lab, P("output", "tables", "S21_label_within_reference_k6.csv"))
write_csv_safe(data.frame(item = c("K", "BLOCK_KM", "B_FULL", "B_LABEL", "seed", "family_rule", "status"),
                          value = c(K, BLOCK_KM, if (SKIP_FULL) "skipped" else B_FULL, B_LABEL, SEED, "first-named author of ref_clean (letters only)",
                                    "EXPLORATORY, post-freeze; site-water results had been seen")), P("output", "tables", "S21_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S21_sessionInfo.txt"))

# ---- console summary ----------------------------------------------------------------------------------------------------------------
cat("\n========== S21 SUMMARY (EXPLORATORY, post-freeze) ==========\n")
cat("\nA1 family structure (S1 sites in domain):\n"); print(fam_sum, row.names = FALSE)
cat("\nTop families by S1 sites:\n"); print(utils::head(ft, 8), row.names = FALSE)
if (!is.null(a2)) { cat("\nA2 S1 k = 6, full-refit CI by cluster definition (S11 block 40 km: 0.694 to 1.660):\n"); print(a2[, c("clusters", "n_clusters", "beta_iqr", "lo", "hi", "n_boot_ok")], row.names = FALSE) }
cat("\nA3 leave-one-family-out, k = 6 (full beta ", round(beta_full, 3), "); five largest changes:\n", sep = ""); print(utils::head(lof[, c("left_out", "n_left_out", "beta_iqr", "change_vs_full")], 5), row.names = FALSE)
cat("range over families: ", round(min(lof$beta_iqr, na.rm = TRUE), 3), " to ", round(max(lof$beta_iqr, na.rm = TRUE), 3), " (", nrow(lof), " families)\n", sep = "")
cat("\nA3 leave-one-block-out, k = 6, six decimals, five largest changes:\n"); print(utils::head(transform(lob[, c("left_out", "n_left_out", "n_remaining", "beta_iqr", "loglik", "tied_beta_6dp")], beta_iqr = round(beta_iqr, 6), loglik = round(loglik, 4)), 5), row.names = FALSE)
cat("blocks with identical beta (6 dp): ", sum(lob$tied_beta_6dp), " | range ", round(min(lob$beta_iqr, na.rm = TRUE), 3), " to ", round(max(lob$beta_iqr, na.rm = TRUE), 3), "\n", sep = "")
if (any(lob$tied_beta_6dp)) { cat("TIED rows:\n"); print(lob[lob$tied_beta_6dp, c("left_out", "n_left_out", "n_remaining", "beta_iqr", "loglik")], row.names = FALSE) }
cat("\nB riverbed label vs log distance (k = 6), log-odds per IQR; negative = riverbed label more likely NEAR large rivers (S13: -1.371 with nuisance, 40 km blocks):\n")
print(lab[, c("model", "n_sites", "n_S1", "n_S3", "n_strata", "log_odds_per_IQR", "lo", "hi", "odds_ratio_near_vs_far", "n_boot_ok")], row.names = FALSE)
fe_row <- lab[lab$model == "references with both labels, fixed effect per reference", ]
cat("\nReading rule (A28): within-reference estimate ", fe_row$log_odds_per_IQR, " (", fe_row$lo, " to ", fe_row$hi, "): ",
    if (fe_row$hi < 0) "label depends on distance WITHIN references (CI excludes 0)" else "CI includes 0 or is positive: not distinguishable from a source convention", "\n", sep = "")
cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("\nDONE S21. Paste this console output back.\n")