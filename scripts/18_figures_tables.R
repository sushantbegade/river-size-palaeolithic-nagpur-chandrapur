# scripts/18_figures_tables.R   (S18 final; REPLACES the earlier version)
# Reads SAVED outputs only (no refits). Every plotted value comes from a saved table; FigSource_*.csv
# lists the plotted values; fig_audit() re-checks every plotted value against the cited table/row/column.
# Map points: data/public/sites_public_1km.csv only (grid verified).
suppressPackageStartupMessages({ library(here); library(sf) })
source(here::here("R", "config.R")); source(here::here("R", "utils.R")); source(here::here("R", "figures.R"))
set.seed(SEED); t0 <- Sys.time()
message("== S18 figures and paper tables (final) ==")
OUT_F <- P("output", "figures"); OUT_T <- P("output", "tables", "paper")
for (d in c(OUT_F, OUT_T, P("output", "logs"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, name) write_csv_safe(x, file.path(OUT_T, name))

# ---- 1. load every required saved table; stop (do not skip silently) if any is missing --------------------------
need <- c("S10_site_domain", "S11_profile", "S11_primary_k6", "S11_decision", "S11_khat_distribution", "S12_sample_sets_profile", "S12_tiers_profile",
          "S12_negative_control", "S13_site_type_profile", "S13_decomposition_k6", "S14_cv_summary", "S15_cluster_definitions", "S15_conditional_profile",
          "S16_variants_profile", "S16_thinning_k6", "S17_distance_shape", "S17_period_profiles", "S17_period_contrast")
tb <- setNames(lapply(need, function(n) rd_csv(paste0(n, ".csv"))), need)
miss <- names(tb)[vapply(tb, is.null, TRUE)]
if (length(miss)) stop("Saved tables not found by rd_csv(): ", paste0(miss, ".csv", collapse = ", "), call. = FALSE)

# ---- 2. build FigSource tables (= what is drawn) and write them ---------------------------------------------------
m1 <- fig1_load(CRS_ANALYSIS)
src <- list(Fig1 = fig1_data(m1), Fig2 = fig2_data(tb), Fig3 = fig3_data(tb), Fig4 = fig4_data(tb), Fig5 = fig5_data(tb), Fig6 = fig6_data(tb))
fs_files <- character(0)
for (nm in names(src)) {                      # x/y are plotted statistics (k, estimate), written as x_pos / y_val so no header looks like a coordinate
  f <- file.path(OUT_T, sprintf("FigSource_%s.csv", nm)); d <- src[[nm]]
  names(d)[names(d) == "x"] <- "x_pos"; names(d)[names(d) == "y"] <- "y_val"
  write_csv_safe(d, f); fs_files <- c(fs_files, f)
}

# ---- 3. audit: written FigSource vs original saved tables ---------------------------------------------------------
aud <- fig_audit(fs_files[-1], reader = rd_csv)          # Fig1 rows are counts of the public file, not table values
if (aud$n_bad > 0) { print(utils::head(aud$bad, 20)); stop("FigSource audit FAILED: ", aud$n_bad, " mismatches. No figure written.", call. = FALSE) }
message(sprintf("Audit passed: %d plotted values (y/lo/hi) match their saved tables exactly.", aud$n_checked))

# ---- 4. draw ------------------------------------------------------------------------------------------------------
bb <- sf::st_bbox(m1$bnd); pad <- 0.04
ext <- list(xlim = c(bb[["xmin"]] - pad * (bb[["xmax"]] - bb[["xmin"]]), bb[["xmax"]] + pad * (bb[["xmax"]] - bb[["xmin"]])),
            ylim = c(bb[["ymin"]] - pad * (bb[["ymax"]] - bb[["ymin"]]), bb[["ymax"]] + pad * (bb[["ymax"]] - bb[["ymin"]])))
outs <- c(fig_save("Fig1_study_area", 7.0, 8.8, function(fam) fig1_draw(m1, fam, ext, CRS_ANALYSIS, inset_pos = "bottomleft", scalebar_pos = "bottomright"), OUT_F),
          fig_save("Fig2_profile", 7.4, 3.9, function(fam) fig2_draw(src$Fig2, fam), OUT_F),
          fig_save("Fig3_distance_shape", 7.2, 8.4, function(fam) fig3_draw(src$Fig3, fam), OUT_F),
          fig_save("Fig4_spatial_cv", 7.6, 3.7, function(fam) fig4_draw(src$Fig4, fam), OUT_F),
          fig_save("Fig5_robustness", 8.6, 6.8, function(fam) fig5_draw(src$Fig5, fam), OUT_F),
          fig_save("Fig6_site_type_model", 9.0, 4.2, function(fam) fig6_draw(src$Fig6, fam), OUT_F))

# ---- 5. paper tables (CSV) + one markdown bundle -------------------------------------------------------------------
Ec <- fig_est_ci; Fm <- fig_fmt
t1 <- tb$S10_site_domain
pr <- tb$S11_primary_k6
t2 <- data.frame(block_km = pr$block_km, method = pr$method, beta_per_IQR = Ec(pr$beta_iqr, pr$lo, pr$hi), intensity_ratio_q25_vs_q75 = Ec(pr$ratio_q25_vs_q75, pr$ratio_lo, pr$ratio_hi),
                 ratio_per_doubling = Ec(pr$per_doubling, pr$per_doubling_lo, pr$per_doubling_hi), n_sites = pr$n_sites, clusters = pr$n_clusters, usable_resamples = pr$n_boot_ok)
p1 <- tb$S11_profile; s12 <- tb$S12_sample_sets_profile; kh <- tb$S11_khat_distribution
g <- function(s) { d <- s12[s12$subset == s, ]; d[order(d$k), ] }
t3 <- data.frame(k = p1$k, S1_full_refit_40km = Ec(p1$beta_iqr, p1$lo, p1$hi), S2_all_one_step_40km = Ec(g("S2_all")$beta_iqr, g("S2_all")$lo, g("S2_all")$hi),
                 S3_riverbed_one_step_40km = Ec(g("S3_riverbed")$beta_iqr, g("S3_riverbed")$lo, g("S3_riverbed")$hi),
                 S1_share_of_bootstrap_best_k = sprintf("%.1f%%", 100 * kh$share_of_bootstrap_best_k[match(p1$k, kh$k)]),
                 status = ifelse(p1$k == 6, "confirmatory", ifelse(p1$k >= 8, "exploratory", "")))
f5 <- src$Fig5
t4 <- data.frame(group = f5$note, set = f5$series, analysis = f5$x_label, beta_per_IQR = Ec(f5$y, f5$lo, f5$hi), interval_method = f5$interval_type)
t4 <- t4[order(match(t4$group, unique(f5$note)), match(t4$analysis, unique(f5$x_label)), t4$set), ]
th <- tb$S16_thinning_k6
t4b <- data.frame(subset = th$subset, thin_km = th$thin_km, replicates = th$replicates, usable = th$usable, mean_n_kept = th$mean_n_kept, median_beta_per_IQR = Fm(th$median_beta_iqr, 3),
                  replicate_range_2.5_97.5 = paste0(Fm(th$q025, 3), " to ", Fm(th$q975, 3)), share_positive = th$share_positive)
nc <- tb$S12_negative_control
t5 <- data.frame(model = nc$model, beta_per_IQR = Ec(nc$beta_iqr, nc$lo, nc$hi), intensity_ratio_q25_vs_q75 = Ec(nc$ratio_q25_vs_q75, nc$ratio_lo, nc$ratio_hi), n_boot_ok = nc$n_boot_ok)
ty <- tb$S13_site_type_profile; ty <- ty[ty$fit_ok %in% TRUE, ]
t6 <- data.frame(model = ty$model, k = ty$k, log_odds_per_IQR = Ec(ty$beta_iqr, ty$lo, ty$hi), odds_ratio_near_vs_far = Ec(ty$ratio_q25_vs_q75, ty$ratio_lo, ty$ratio_hi),
                 p_riverbed_at_q25 = Fm(ty$p_riverbed_at_q25, 3), p_riverbed_at_q75 = Fm(ty$p_riverbed_at_q75, 3), n = ty$n, n_riverbed = ty$n_riverbed)
cv <- tb$S14_cv_summary; cv <- cv[order(cv$set, -cv$block_km, match(cv$model, unique(cv$model))), ]
t7 <- data.frame(set = cv$set, block_km = cv$block_km, model = cv$model, n_folds_used = cv$n_folds_used, delta_loglik_vs_nuisance = Fm(cv$delta_vs_nuisance, 2),
                 se_delta_sum = Fm(cv$se_delta_sum, 2), delta_per_site = Fm(cv$delta_per_site, 3), folds_improved = cv$folds_improved, folds_worse = cv$folds_worse)
pp <- tb$S17_period_profiles; pp <- pp[pp$k == 6, ]; pc <- tb$S17_period_contrast
t8 <- data.frame(sample_set = pp$sample_set, period = pp$period, n_sites = pp$n_sites, n_clusters = pp$n_clusters, beta_per_IQR = Ec(pp$beta_iqr, pp$lo, pp$hi))
t8b <- data.frame(sample_set = pc$sample_set, contrast = pc$contrast, difference = Ec(pc$difference, pc$lo, pc$hi))
sh <- tb$S17_distance_shape
t9 <- data.frame(subset = sh$subset, k = sh$k, class = sh$bin, km_from = sh$km_from, km_to = sh$km_to, n_sites = sh$n_sites, sites_per_100km2 = sh$sites_per_100km2,
                 rel_intensity = Fm(sh$rel_intensity, 2), ci = ifelse(is.na(sh$lo), "reference", paste0(Fm(sh$lo, 2), " to ", Fm(sh$hi, 2))))
tabs <- list(T1_samples_in_analysis_domain = t1, T2_primary_k6_S1 = t2, T3_profile_by_sample_set = t3, T4_robustness_k6 = t4, T4b_thinning_k6 = t4b,
             T5_negative_control = t5, T6_site_type_model_EXPLORATORY = t6, T7_spatial_cv = t7, T8_single_period_k6_EXPLORATORY = t8,
             T8b_period_contrasts_EXPLORATORY = t8b, T9_distance_shape_EXPLORATORY = t9)
for (n in names(tabs)) wr(tabs[[n]], paste0(n, ".csv"))
h1 <- tb$S11_decision$value[grepl("^H1 supported", tb$S11_decision$item)]
ttl <- c(T1_samples_in_analysis_domain = "Table 1. Samples in the analysis domain.",
         T2_primary_k6_S1 = paste0("Table 2. Confirmatory test, k = 6, S1 open-air. Pre-registered rule result (H1 supported): ", h1, ". beta > 0 means intensity RISES with distance (opposite to H1)."),
         T3_profile_by_sample_set = "Table 3. Threshold profile. S1: full-refit bootstrap; S2, S3: one-step bootstrap (40 km blocks). Profile claims at +-1 k only; k8-k9 exploratory.",
         T4_robustness_k6 = "Table 4. Robustness at k = 6. Interval methods differ by row (last column); 'extra' rows were not pre-registered. S2 has no tier or full-refit rows because none were computed.",
         T4b_thinning_k6 = "Table 4b. Spatial thinning, k = 6 (100 replicates). Range = 2.5-97.5% of replicate estimates, not a confidence interval.",
         T5_negative_control = "Table 5. Negative control (lineament distance), S1, 40 km blocks.",
         T6_site_type_model_EXPLORATORY = "Table 6. EXPLORATORY. Log-odds of riverbed vs open-air per IQR of log distance (cluster bootstrap).",
         T7_spatial_cv = "Table 7. Spatial cross-validation: held-out log-likelihood gain over the nuisance model; SE is of the summed fold differences (not a CI).",
         T8_single_period_k6_EXPLORATORY = "Table 8. EXPLORATORY. Single-period estimates, k = 6. Do not compare AIC across periods.",
         T8b_period_contrasts_EXPLORATORY = "Table 8b. EXPLORATORY. Period contrasts; all intervals include 0, so no period claim.",
         T9_distance_shape_EXPLORATORY = "Table 9. EXPLORATORY. Intensity relative to the reference class (class 3) in five equal-area distance classes.")
md <- c("# Paper tables (generated by S18; all values from saved tables)", "")
for (n in names(tabs)) md <- c(md, paste0("**", ttl[[n]], "**"), "", fig_md_table(tabs[[n]]), "")
writeLines(md, file.path(OUT_T, "paper_tables.md"))

# ---- 6. captions (counts are computed from the plotted data, not typed) ---------------------------------------------
a5 <- f5[f5$panel == "a" & f5$x > 0, ]; b5 <- f5[f5$panel == "b", ]
n1 <- sum(a5$lo > 0 | a5$hi < 0); n2 <- sum(b5$lo <= 0 & b5$hi >= 0)
cap <- c("# Figure captions (draft; every number is in FigSource_FigN.csv)", "",
         "**Fig. 1.** Study area (Nagpur and Chandrapur districts) in true shape and scale (equal-aspect projected map, EPSG:32644; graticule in degrees). WRIS rivers by order (< 5, 5, >= 6), removed areas (ponds, reservoirs, quarries + 100 m) and findspots by label. Findspot positions are the 1 km coarsened public version. Inset: location in India (boundary file supplied by the authors).", "",
         "**Fig. 2.** Log intensity ratio per IQR of log distance to drainage of WRIS order >= k (positive = intensity rises with distance). S1 (open-air): full-refit cluster bootstrap, 40 km blocks. S2 and S3: one-step approximation. H1 predicted a negative value in S1; the opposite was found. S3 is at water by construction. k = 8-9 are exploratory.", "",
         "**Fig. 3.** EXPLORATORY. Intensity relative to the middle distance class (five equal-area classes of distance to rivers of order >= k, nuisance-adjusted; 95% one-step cluster bootstrap intervals as saved; log axis). Open symbol = reference class.", "",
         "**Fig. 4.** Spatial cross-validation: held-out log-likelihood gain of each river model (and the lineament negative control) over the nuisance-only model, 40 km (filled) and 20 km (open) blocks. Whiskers are +-1 SE of the summed fold differences, a rough scale and not a confidence interval.", "",
         "**Fig. 5.** Robustness of the k = 6 estimate. Left, S1; right, S2. Interval methods differ by group (see Table 4). Variants labelled extra were not pre-registered. S2 has no coordinate-tier or full-refit rows because none were computed.",
         sprintf("Of %d S1 intervals shown, %d exclude zero; of %d S2 intervals shown, %d include zero.", nrow(a5), n1, nrow(b5), n2), "",
         "**Fig. 6.** EXPLORATORY and post hoc. (a) Log-odds that a site is riverbed rather than open-air per IQR of log distance (negative = riverbed sites nearer rivers). (b) Decomposition at k = 6. Because the label depends on river distance, S1 and S2 are not classification-independent tests of river-size preference.")
writeLines(cap, file.path(OUT_T, "captions.md"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S18_sessionInfo.txt"))

# ---- 7. manifest + console summary ----------------------------------------------------------------------------------
allout <- c(outs, list.files(OUT_T, full.names = TRUE)); allout <- allout[!grepl("S18_outputs_manifest", allout)]
wr(data.frame(file = basename(allout), bytes = file.size(allout), md5 = unname(tools::md5sum(allout))), "S18_outputs_manifest.csv")
kv <- function(df, cond) df[cond, ]
cat("\n========== S18 SUMMARY ==========\n")
cat("Figures written:", paste(basename(outs), collapse = ", "), "\n")
cat("Paper tables:", length(tabs), "CSV + paper_tables.md + captions.md in output/tables/paper\n")
cat(sprintf("Audit: %d plotted values checked against saved tables, %d mismatches.\n", aud$n_checked, aud$n_bad))
cat(sprintf("Fig1: %d public 1 km points drawn (%s); rivers %d features; India inset %s.\n", nrow(m1$pts), paste(names(table(m1$pts$cls)), table(m1$pts$cls), collapse = ", "), nrow(m1$rivers), if (is.null(m1$india)) "NOT drawn (file missing)" else "drawn"))
cat("\nVALUES TO CHECK BY EYE ON THE FIGURES (from FigSource):\n")
f2 <- src$Fig2; f4 <- src$Fig4; f3 <- src$Fig3; f6 <- src$Fig6
pk <- function(d, cond, lab) { r <- d[cond, ][1, ]; cat(sprintf("  %-46s %s (%s to %s)\n", lab, Fm(r$y, 3), Fm(r$lo, 3), Fm(r$hi, 3))) }
pk(f2, f2$panel == "a" & f2$x == 6, "Fig2 S1 k6 beta (full refit)")
pk(f2, f2$panel == "a" & f2$x == 8, "Fig2 S1 k8 beta (lower limit must be -0.002)")
pk(f2, f2$series == "S2" & f2$x == 3, "Fig2 S2 k3")
pk(f2, f2$series == "S3" & f2$x == 3, "Fig2 S3 k3")
pk(f3, f3$panel == "b" & f3$x == 1, "Fig3 S1 k6 class 1 (rel. intensity)")
pk(f4, f4$series == "S1 40 km" & f4$x_label == "k6", "Fig4 S1 k6 gain (+-1 SE)")
pk(f4, f4$series == "S1 40 km" & f4$x_label == "Lin", "Fig4 S1 lineament gain (+-1 SE)")
for (lab in c("baseline (one-step comparator)", "5 km survey domain", "reference-conditional likelihood", "tier B (n = 32)")) pk(f5, f5$panel == "a" & f5$x_label == lab, paste("Fig5 S1", lab))
pk(f5, f5$panel == "b" & f5$x_label == "baseline (one-step comparator)", "Fig5 S2 baseline")
pk(f6, f6$series == "distance + nuisance" & f6$x == 6, "Fig6 nuisance k6 log-odds")
pk(f6, f6$series == "distance only" & f6$x == 3, "Fig6 distance-only k3 log-odds")
cat("\nData safety: map uses data/public/sites_public_1km.csv only (grid verified); FigSource_Fig1 holds counts only; no exact coordinate is drawn or tabulated.\n")
cat("Elapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S18. Paste this console output back, and open each figure to check it before use.\n")