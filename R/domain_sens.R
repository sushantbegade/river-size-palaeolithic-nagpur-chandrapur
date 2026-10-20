# R/domain_sens.R -- functions for script 20 (EXPLORATORY, post-freeze): sensitivity of the S1 result to domain construction, plus k = 2.
# Uses the frozen helpers unchanged. Needs: config.R, quadrature.R (survey_domain, nearest_dist, size_line_sets, log_d),
# stepzero.R, stepzero_b.R (cluster_ids), stepzero_c.R, model_primary.R (design_from_domain, exposure_domain, fit_real, boot_refit, ci_row, nuisance_raw),
# diagnostics.R (subset_profile).

# first-named author (letters only) of a reference string; composite strings ("A, 1961; B, 2018") are assigned to their FIRST author
ref_family <- function(ref_clean) {
  x <- trimws(as.character(ref_clean))
  f <- sub("^([A-Za-z]+).*$", "\\1", x)
  ifelse(grepl("^[A-Za-z]", x), f, x)
}

in_poly <- function(pts, poly) lengths(sf::st_intersects(pts, poly)) > 0

# survey-domain polygon from a SUBSET of sites; same rule as the baseline (hull if >= SURVEY_MIN_N sites of a reference, else points; buffer; clip)
subset_domain_poly <- function(sites_sf, keep, study) {
  stopifnot(is.logical(keep), length(keep) == nrow(sites_sf), !anyNA(keep), any(keep))
  survey_domain(sites_sf[keep, ], study)
}

# design for a domain: identical to design_from_domain() unless the domain has no alluvial cells
# (the frozen nuisance_design() stops on a constant column; the alluvial column is dropped from the model anyway)
design_for_domain <- function(dom) {
  if (any(dom$gm4 == "alluvial", na.rm = TRUE)) return(design_from_domain(dom))
  keep <- c("elev", "slope", "gm_plateau", "crop", "built", "road", "settle")
  R <- nuisance_raw(dom)[, keep, drop = FALSE]
  sds <- apply(R, 2, stats::sd)
  if (any(!is.finite(sds) | sds == 0)) stop("Constant or NA nuisance column: ", paste(names(sds)[!is.finite(sds) | sds == 0], collapse = ", "), call. = FALSE)
  cen <- colMeans(R)
  list(dom = sweep(sweep(R, 2, cen, "-"), 2, sds, "/"), center = cen, scale = sds,
       apply = function(d) sweep(sweep(nuisance_raw(d)[, keep, drop = FALSE], 2, cen, "-"), 2, sds, "/"))
}

# IQR of log distance (over the domain) for each k
iqr_log_by_k <- function(dom, ks) {
  stats::setNames(vapply(ks, function(k) stats::IQR(log_d(dom[[paste0("d_k", k)]])), 0), as.character(ks))
}

# one-step profile for one domain and one site set (S1): ridge fit, one-step 40 km cluster-bootstrap CI, AIC; effect per IQR of THIS domain
profile_domain <- function(label, dm, s1, ks, des, block_km = 40, B = 2000L) {
  ex <- stats::setNames(lapply(ks, function(k) exposure_domain(dm, k)), as.character(ks))
  Xq <- stats::setNames(lapply(ks, function(k) cbind(int = 1, des$dom, z = ex[[as.character(k)]]$z_dom)), as.character(ks))
  r <- subset_profile(label, s1, ks, des, ex, Xq, dm$area_km2, block_km, B)
  r$iqr_log <- unname(vapply(as.character(ks), function(k) ex[[k]]$iqr_per_sd * ex[[k]]$sd_log, 0)[as.character(r$k)])
  r
}

# beta, lo and hi are linear in the IQR used: beta on the BASELINE-domain IQR = beta_own x IQR_base / IQR_domain
to_base_scale <- function(prof, iqr_base_by_k) {
  f <- unname(iqr_base_by_k[as.character(prof$k)]) / prof$iqr_log
  prof$beta_base <- round(prof$beta_iqr * f, 3); prof$lo_base <- round(prof$lo * f, 3); prof$hi_base <- round(prof$hi * f, 3)
  prof$ratio_q25_vs_q75_base <- round(exp(-prof$beta_base), 2)
  prof
}

# full-refit cluster bootstrap at one k (same estimator and bootstrap as S11), with the CI on the domain's own IQR and on the baseline IQR
fullrefit_k <- function(dm, s1, k, des, iqr_base, block_km = 40, B = 1000L) {
  kk <- as.character(k)
  ex <- exposure_domain(dm, k); a <- dm$area_km2
  Xq <- cbind(int = 1, des$dom, z = ex$z_dom)
  Xp <- cbind(int = 1, des$apply(s1), z = ex$apply(s1))
  f <- fit_real(Xp, Xq, a)
  if (f$bad) return(NULL)
  cl <- cluster_ids(as.matrix(s1[, c("x", "y")]), block_km * 1000)
  bt <- boot_refit(stats::setNames(list(Xq), kk), stats::setNames(list(Xp), kk), a, cl, B, progress = 250L)
  est <- f$coef[length(f$coef)]
  own <- ci_row(bt$beta[, kk], ex$iqr_per_sd, ex$sd_log, est)
  bas <- ci_row(bt$beta[, kk], iqr_base / ex$sd_log, ex$sd_log, est)
  data.frame(k = k, n_sites = nrow(s1), n_clusters = length(unique(cl)), beta_iqr = own$beta_iqr, lo = own$lo, hi = own$hi,
             beta_base = bas$beta_iqr, lo_base = bas$lo, hi_base = bas$hi, ratio_q25_vs_q75_base = bas$ratio_q25_vs_q75,
             n_boot_ok = own$n_boot_ok, B = B, iqr_log_domain = round(ex$iqr_per_sd * ex$sd_log, 3), iqr_log_base = round(iqr_base, 3))
}

domain_composition <- function(label, dm, s1, base_area, k = 6L, near_m = 5000, block_km = 40) {
  dk <- paste0("d_k", k); v <- log_d(dm[[dk]])
  data.frame(domain = label, n_cells = nrow(dm), area_km2 = round(sum(dm$area_km2), 1), share_of_baseline_area = round(sum(dm$area_km2) / base_area, 3),
             n_S1 = nrow(s1), n_ref_S1 = length(unique(s1$ref_id)),
             n_clusters_40km = length(unique(cluster_ids(as.matrix(s1[, c("x", "y")]), block_km * 1000))),
             share_area_within_5km_k6 = round(mean(dm[[dk]] <= near_m), 3), share_S1_within_5km_k6 = round(mean(s1[[dk]] <= near_m), 3),
             sd_log_d_k6 = round(stats::sd(v), 3), iqr_log_d_k6 = round(stats::IQR(v), 3))
}

# data safety: refuse to write a table with coordinate-like columns
write_s20 <- function(x, path) {
  bad <- grepl("^(x|y|x_?c|y_?c|lon|lat|long|longitude|latitude|easting|northing|coord.*|utm.*)$", tolower(names(x)))
  if (any(bad)) stop("Refusing to write coordinate-like column(s): ", paste(names(x)[bad], collapse = ", "), call. = FALSE)
  write_csv_safe(x, path)
}