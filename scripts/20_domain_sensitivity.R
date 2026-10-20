# scripts/20_domain_sensitivity.R
# S20 (EXPLORATORY, post-freeze; site-water results HAD been seen). Two questions:
#  (1) Does the S1 k = 6 association depend on how the survey domain was built? The baseline domain uses hulls around ALL sites of
#      every reference (S1, S3, S4). Same 88 S1 sites, different quadrature domains, all subsets of the baseline domain.
#  (2) What does k = 2 look like (never analysed; no recorded reason)?
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
source(here::here("R", "domain_sens.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S20 domain-construction sensitivity and k = 2 (EXPLORATORY, post-freeze) ==")

KS        <- 3:9
KS2       <- 2:9
BLOCK_KM  <- 40
B_ONE     <- as.integer(Sys.getenv("S20_B_ONESTEP", "2000"))
B_FULL    <- as.integer(Sys.getenv("S20_B_FULL", "1000"))
SKIP_FULL <- Sys.getenv("S20_SKIP_FULL", "0") == "1"
INT_KM    <- as.numeric(Sys.getenv("S20_INTERIOR_KM", "15"))        # D5 only; not a frozen value
EXPECT_S1 <- as.integer(Sys.getenv("S20_EXPECT_S1", "88"))
EXPECT_CL <- stats::setNames(as.integer(strsplit(Sys.getenv("S20_EXPECT_CL", "16,32,52"), ",")[[1]]), c("40", "20", "10"))  # S1 occupied blocks (S10_cluster_structure.csv)

# reading rule for (1), written BEFORE the run (deviation log A25): the S1 k = 6 association is "domain-robust" if, in each of
# D1, D2, D2f, D3 and D4, the full-refit 40 km CI on the BASELINE IQR scale excludes 0 with beta > 0. D5 is reported, not used in the rule.
RULE_DOMAINS <- c("D1", "D2", "D2f", "D3", "D4")

f <- list(q = P("data", "processed", "quadrature_det.rds"), s = P("data", "processed", "sites_cov.rds"),
          h = P("data", "processed", "hydro_prep.gpkg"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S02, S08 and S10 first.", call. = FALSE)

grid   <- sf::st_drop_geometry(readRDS(f$q))
sc_all <- readRDS(f$s)
bnd    <- sf::st_read(f$h, layer = "boundary", quiet = TRUE)
study  <- sf::st_union(sf::st_geometry(bnd))
rv     <- sf::st_read(f$h, layer = "rivers", quiet = TRUE)
chan   <- if ("removed_channel" %in% sf::st_layers(f$h)$name)
  sf::st_union(sf::st_geometry(sf::st_read(f$h, layer = "removed_channel", quiet = TRUE))) else sf::st_sfc(sf::st_polygon(), crs = CRS_ANALYSIS)

need_g <- c("x", "y", "in_survey", "in_core", "area_km2", "elev_m", "slope_deg", "gm4", "frac_crop", "frac_built", "d_road_major", "d_settle",
            paste0("d_k", SIZE_K))
need_s <- c("site_id", "ref_id", "ref_clean", "x", "y", "S1_openair", "in_core", "in_domain", "gm4", "elev_m", "slope_deg", "frac_crop", "frac_built",
            "d_road_major", "d_settle", paste0("d_k", SIZE_K))
if (!all(need_g %in% names(grid))) stop("quadrature_det.rds lacks: ", paste(setdiff(need_g, names(grid)), collapse = ", "), call. = FALSE)
if (!all(need_s %in% names(sc_all))) stop("sites_cov.rds lacks: ", paste(setdiff(need_s, names(sc_all)), collapse = ", "), call. = FALSE)

grid_pts  <- sf::st_as_sf(grid[, c("x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
sites_sf  <- sf::st_as_sf(sc_all[, c("site_id", "ref_id", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
has_chan  <- length(chan) && !all(sf::st_is_empty(chan))
if (!"in_channel" %in% names(grid)) { grid$in_channel <- if (has_chan) in_poly(grid_pts, chan) else FALSE; message("in_channel rebuilt from removed_channel layer") }
s1flag    <- sc_all$S1_openair %in% TRUE
if (anyNA(sc_all$ref_id)) stop("NA ref_id among sites.", call. = FALSE)

# ---- (0) baseline reconstruction checks -------------------------------------------------------------------------------------
P0 <- survey_domain(sites_sf, study)
chk_cells <- in_poly(grid_pts, P0)
if (!identical(unname(chk_cells), unname(grid$in_survey))) stop("Baseline survey domain NOT reproduced: ", sum(chk_cells != grid$in_survey), " cells differ.", call. = FALSE)
chk_sites <- in_poly(sites_sf, P0) & !sc_all$in_core
if (!identical(unname(chk_sites), unname(sc_all$in_domain))) stop("Site domain membership NOT reproduced: ", sum(chk_sites != sc_all$in_domain), " sites differ.", call. = FALSE)
base <- grid$in_survey & !grid$in_core
message("baseline reproduced: ", sum(base), " domain cells | S1 in domain: ", sum(s1flag & sc_all$in_domain))
if (sum(s1flag & sc_all$in_domain) != EXPECT_S1) stop("Baseline S1 n differs from expected (", EXPECT_S1, ").", call. = FALSE)
s1_0 <- sc_all[s1flag & sc_all$in_domain, ]
for (b in names(EXPECT_CL)) {
  n_cl <- length(unique(cluster_ids(as.matrix(s1_0[, c("x", "y")]), as.numeric(b) * 1000)))
  if (n_cl != EXPECT_CL[[b]]) stop("S1 clusters at ", b, " km: ", n_cl, " (expected ", EXPECT_CL[[b]], "). Cluster ids or data changed.", call. = FALSE)
}

# ---- (1) domain variants: subsets of the baseline ---------------------------------------------------------------------------
fam      <- ref_family(sc_all$ref_clean)
refs_s1  <- unique(sc_all$ref_id[s1flag]); fam_s1 <- unique(fam[s1flag])
keep_set <- list(D1 = s1flag, D2 = sc_all$ref_id %in% refs_s1, D2f = fam %in% fam_s1)
polys    <- lapply(keep_set, function(k) subset_domain_poly(sites_sf, k, study))
cells <- list(D0 = base)
site_in <- list(D0 = sc_all$in_domain)
for (nm in names(polys)) {
  cells[[nm]]   <- in_poly(grid_pts, polys[[nm]]) & !grid$in_core
  site_in[[nm]] <- in_poly(sites_sf, polys[[nm]]) & !sc_all$in_core
}
site_chan <- if (has_chan) in_poly(sites_sf, chan) else rep(FALSE, nrow(sc_all))
cells$D3   <- base & !grid$in_channel;                  site_in$D3 <- sc_all$in_domain & !site_chan
cells$D4   <- base & !(grid$gm4 %in% "alluvial");       site_in$D4 <- sc_all$in_domain & !(sc_all$gm4 %in% "alluvial")
bl         <- sf::st_boundary(study)
d_bnd      <- rep(NA_real_, nrow(grid)); d_bnd[base] <- as.numeric(sf::st_distance(grid_pts[base, ], bl))
d_bnd_site <- as.numeric(sf::st_distance(sites_sf, bl))
cells$D5   <- base & !is.na(d_bnd) & d_bnd >= INT_KM * 1000; site_in$D5 <- sc_all$in_domain & d_bnd_site >= INT_KM * 1000
cells <- cells[c("D0", "D1", "D2", "D2f", "D3", "D4", "D5")]

# buffers of different hulls are polygon approximations of arcs (about 3 m sagitta at 10 km), so a cell centre within metres of an edge can fall
# inside a subset polygon but outside the baseline polygon. Each variant is therefore intersected with the baseline; a larger discrepancy stops the run.
for (nm in setdiff(names(cells), "D0")) {
  n_out <- sum(cells[[nm]] & !base)
  if (n_out > 0.001 * sum(cells[[nm]])) stop(nm, " is not a subset of the baseline domain (", n_out, " cells outside).", call. = FALSE)
  if (n_out > 0) message("  ", nm, ": ", n_out, " edge cell(s) outside the baseline polygon removed (arc approximation)")
  cells[[nm]]   <- cells[[nm]] & base
  if (any(site_in[[nm]] & !sc_all$in_domain)) stop(nm, ": a site is in the domain but not in the baseline domain.", call. = FALSE)
}
n_s1 <- vapply(names(cells), function(nm) sum(s1flag & site_in[[nm]]), 0L)
for (nm in c("D1", "D2", "D2f")) if (n_s1[[nm]] != EXPECT_S1) stop(nm, ": S1 n = ", n_s1[[nm]], " (expected ", EXPECT_S1, ").", call. = FALSE)
message("domains: ", paste(sprintf("%s %d cells / S1 %d", names(cells), vapply(cells, sum, 0L), n_s1), collapse = " | "))

dom_of <- function(nm) grid[cells[[nm]], ]
s1_of  <- function(nm) sc_all[s1flag & site_in[[nm]], ]
base_area <- sum(grid$area_km2[base])
iqr_base  <- iqr_log_by_k(dom_of("D0"), KS)

comp <- do.call(rbind, lapply(names(cells), function(nm) domain_composition(nm, dom_of(nm), s1_of(nm), base_area)))
write_s20(comp, P("output", "tables", "S20_domain_composition.csv"))

# ---- one-step profile k = 3..9 in every domain --------------------------------------------------------------------------------
des_of <- list()
prof <- list(); full <- list()
for (i in seq_along(cells)) {
  nm <- names(cells)[i]; dm <- dom_of(nm); s1 <- s1_of(nm)
  des_of[[nm]] <- design_for_domain(dm)
  set.seed(SEED + i)
  p <- profile_domain(nm, dm, s1, KS, des_of[[nm]], BLOCK_KM, B_ONE)
  prof[[nm]] <- to_base_scale(p, iqr_base)
  if (!SKIP_FULL) {
    set.seed(SEED + 100L + i)
    fr <- fullrefit_k(dm, s1, 6L, des_of[[nm]], iqr_base[["6"]], BLOCK_KM, B_FULL)
    if (!is.null(fr)) full[[nm]] <- cbind(domain = nm, fr)
  }
  message("  done: ", nm, "  (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
}
prof <- do.call(rbind, prof); rownames(prof) <- NULL
prof <- prof[, c("subset", "k", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "beta_base", "lo_base", "hi_base", "ratio_q25_vs_q75_base", "aic", "fit_ok", "n_boot_ok")]
names(prof)[1] <- "domain"
write_s20(prof, P("output", "tables", "S20_domain_profile.csv"))
if (length(full)) { full <- do.call(rbind, full); rownames(full) <- NULL; write_s20(full, P("output", "tables", "S20_domain_k6_fullrefit.csv")) }

# ---- (2) k = 2 ------------------------------------------------------------------------------------------------------------------
rv$length_km <- as.numeric(sf::st_length(rv)) / 1000
rin <- if ("in_study" %in% names(rv)) rv[rv$in_study %in% TRUE, ] else rv
ord <- sf::st_drop_geometry(rin) |> dplyr::filter(!is.na(ordsh)) |> dplyr::group_by(ordsh) |>
  dplyr::summarise(n_segments = dplyr::n(), length_km = round(sum(length_km), 1), .groups = "drop") |> dplyr::arrange(ordsh)
ord$share_of_length <- round(ord$length_km / sum(ord$length_km), 4)
ord$length_km_order_ge <- round(rev(cumsum(rev(ord$length_km))), 1)
write_s20(as.data.frame(ord), P("output", "tables", "S20_k2_orders.csv"))

empty_hr <- sf::st_sf(HYRIV_ID = integer(0), UPLAND_SKM = numeric(0), geometry = sf::st_sfc(crs = CRS_ANALYSIS))
lines_k  <- function(k) size_line_sets(rv, empty_hr, numeric(0), ks = as.integer(k))$with_ext[[as.character(k)]]
dm0 <- dom_of("D0"); pts0 <- grid_pts[cells$D0, ]
d3  <- nearest_dist(pts0, lines_k(3L))
if (max(abs(d3 - dm0$d_k3)) > 1e-6) stop("Rebuilt d_k3 differs from stored d_k3 (max ", max(abs(d3 - dm0$d_k3)), " m): river layer or logic changed.", call. = FALSE)
message("rebuilt d_k3 equals stored d_k3 over the domain (check passed)")
dm0$d_k2 <- nearest_dist(pts0, lines_k(2L))
sc_all$d_k2 <- nearest_dist(sites_sf, lines_k(2L))
s1_k2 <- sc_all[s1flag & site_in$D0, ]
lg <- function(v) log_d(v)
cor3 <- function(a, b) round(stats::cor(lg(a), lg(b)), 3)
k2sum <- data.frame(
  item = c("domain cells", "median d_k1 (km)", "median d_k2 (km)", "median d_k3 (km)", "IQR log d_k1", "IQR log d_k2", "IQR log d_k3",
           "cor log d_k1 vs d_k2", "cor log d_k2 vs d_k3", "cor log d_k1 vs d_k3",
           "share of cells with d_k2 = d_k1 (within 1 m)", "share of cells with d_k2 = d_k3 (within 1 m)",
           "S1 sites: median d_k2 (km)", "S1 sites: share within 5 km of d_k2"),
  value = c(nrow(dm0), round(stats::median(dm0$d_k1) / 1000, 3), round(stats::median(dm0$d_k2) / 1000, 3), round(stats::median(dm0$d_k3) / 1000, 3),
            round(stats::IQR(lg(dm0$d_k1)), 3), round(stats::IQR(lg(dm0$d_k2)), 3), round(stats::IQR(lg(dm0$d_k3)), 3),
            cor3(dm0$d_k1, dm0$d_k2), cor3(dm0$d_k2, dm0$d_k3), cor3(dm0$d_k1, dm0$d_k3),
            round(mean(abs(dm0$d_k2 - dm0$d_k1) < 1), 3), round(mean(abs(dm0$d_k2 - dm0$d_k3) < 1), 3),
            round(stats::median(s1_k2$d_k2) / 1000, 3), round(mean(s1_k2$d_k2 <= 5000), 3)))
write_s20(k2sum, P("output", "tables", "S20_k2_summary.csv"))

set.seed(SEED + 200L)
des0 <- des_of$D0
p2 <- profile_domain("D0_with_k2", dm0, s1_k2, KS2, des0, BLOCK_KM, B_ONE)
p2 <- p2[, c("subset", "k", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "ratio_q25_vs_q75", "aic", "fit_ok", "n_boot_ok")]
names(p2)[1] <- "domain"
write_s20(p2, P("output", "tables", "S20_k2_profile.csv"))
fr2 <- NULL
if (!SKIP_FULL) {
  set.seed(SEED + 300L)
  fr2 <- fullrefit_k(dm0, s1_k2, 2L, des0, stats::IQR(lg(dm0$d_k2)), BLOCK_KM, B_FULL)
  if (!is.null(fr2)) write_s20(fr2[, c("k", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "n_boot_ok", "B")], P("output", "tables", "S20_k2_fullrefit.csv"))
}

write_s20(data.frame(item = c("KS", "KS_with_k2", "BLOCK_KM", "B_ONESTEP", "B_FULL", "INTERIOR_KM", "PRIOR_SD", "seed", "family_rule", "RULE_DOMAINS", "status"),
                     value = c(paste(KS, collapse = ","), paste(KS2, collapse = ","), BLOCK_KM, B_ONE, if (SKIP_FULL) "skipped" else B_FULL, INT_KM, PRIOR_SD, SEED,
                               "first-named author of ref_clean (letters only)", paste(RULE_DOMAINS, collapse = ","),
                               "EXPLORATORY, post-freeze; site-water results had been seen")),
          P("output", "tables", "S20_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S20_sessionInfo.txt"))

# ---- console summary ------------------------------------------------------------------------------------------------------------
cat("\n========== S20 SUMMARY (EXPLORATORY, post-freeze) ==========\n")
cat("\nDomain composition (the S1 site set is identical in D0, D1, D2 and D2f; D3 to D5 may drop sites):\n"); print(comp, row.names = FALSE)
cat("\nS1 profile, one-step 40 km CI. beta_iqr = per IQR of the domain's own log distance; *_base = per BASELINE IQR:\n")
print(prof[prof$k %in% c(4, 6, 8), c("domain", "k", "n_sites", "beta_iqr", "lo", "hi", "beta_base", "lo_base", "hi_base", "aic")], row.names = FALSE)
if (length(full)) {
  cat("\nk = 6 full refit (", B_FULL, " resamples, 40 km blocks):\n", sep = ""); print(full[, c("domain", "n_sites", "n_clusters", "beta_iqr", "lo", "hi", "beta_base", "lo_base", "hi_base", "n_boot_ok")], row.names = FALSE)
  ok <- full[full$domain %in% RULE_DOMAINS, ]
  cat("\nReading rule (A25): domain-robust iff beta_base > 0 and lo_base > 0 in each of ", paste(RULE_DOMAINS, collapse = ", "), "\n", sep = "")
  rr <- data.frame(domain = ok$domain, beta_base = ok$beta_base, lo_base = ok$lo_base, met = ok$beta_base > 0 & ok$lo_base > 0)
  print(rr, row.names = FALSE)
  cat("RULE RESULT:", if (all(rr$met) && nrow(rr) == length(RULE_DOMAINS)) "domain-robust" else "SENSITIVE to domain construction", "\n")
  s11 <- P("output", "tables", "S11_primary_k6.csv")
  if (file.exists(s11)) { t11 <- utils::read.csv(s11); cat("\nS11 primary table (for comparison with D0; random draws differ):\n"); print(t11[, seq_len(min(6, ncol(t11)))], row.names = FALSE) }
}
cat("\nk = 2 (WRIS order >= 2) vs k = 1 and k = 3 over the baseline domain:\n"); print(k2sum, row.names = FALSE)
cat("\nWRIS length by ordsh inside the study:\n"); print(as.data.frame(ord), row.names = FALSE)
cat("\nS1 profile k = 2..9 in D0 (one-step 40 km CI):\n"); print(p2[, c("k", "n_sites", "beta_iqr", "lo", "hi", "aic")], row.names = FALSE)
if (!is.null(fr2)) { cat("\nk = 2 full refit (own IQR scale):\n"); print(fr2[, c("k", "beta_iqr", "lo", "hi", "n_boot_ok")], row.names = FALSE) }
cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("\nDONE S20. Paste this console output back.\n")