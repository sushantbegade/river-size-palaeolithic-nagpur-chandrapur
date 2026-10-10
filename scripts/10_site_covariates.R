# scripts/10_site_covariates.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf); library(terra)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hydrography_prep.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "covariates.R"))
source(here::here("R", "detectability.R"))
source(here::here("R", "osm_xml.R"))
source(here::here("R", "external_data.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))   # cluster_ids
source(here::here("R", "site_covariates.R"))

set.seed(SEED)
t0 <- Sys.time()
message("== S10 site covariates: extraction + QA + structure counts (BLIND: no covariate-site relationship is computed) ==")

BLOCKS_KM <- c(10, 20, 40)
QA_N      <- 1500L

# ---- inputs ----------------------------------------------------------------------------------------
f <- list(sites = P("data", "processed", "sites_clean.rds"), h = P("data", "processed", "hydro_prep.gpkg"),
          e = P("data", "processed", "external_region.gpkg"), q = P("data", "processed", "quadrature_det.rds"),
          d = P("data", "processed", "domains.gpkg"), t = P("data", "processed", "terrain_510m.tif"),
          w = P("data", "processed", "worldcover_frac_100m.tif"), cal = P("output", "tables", "S05_wris_order_vs_hr_upland.csv"))
for (x in f) if (!file.exists(x)) stop("Missing: ", x, "\n  Run S01-S08 first.", call. = FALSE)

sites <- readRDS(f$sites)
bnd   <- sf::st_read(f$h, layer = "boundary", quiet = TRUE)
study <- sf::st_union(sf::st_geometry(bnd))
rv    <- sf::st_read(f$h, layer = "rivers", quiet = TRUE)
gm    <- sf::st_read(f$h, layer = "geomorph", quiet = TRUE)
core  <- sf::st_union(sf::st_geometry(sf::st_read(f$h, layer = "removed_core", quiet = TRUE)))
hr    <- sf::st_read(f$e, layer = "hydrorivers", quiet = TRUE)
dsurv <- sf::st_union(sf::st_geometry(sf::st_read(f$d, layer = "survey_domain", quiet = TRUE)))
grid  <- readRDS(f$q)
cal   <- readr::read_csv(f$cal, show_col_types = FALSE)
upl   <- stats::setNames(cal$upland_median, as.character(cal$ordsh))

# ---- context: the SAME surfaces used for the grid -------------------------------------------------------
hro  <- hr_outside(hr, study)
sets <- size_line_sets(rv, hro, upl)
lin_path <- DATA("Linement", "Linement.shp")
if (!file.exists(lin_path)) stop("Lineament file not found: ", lin_path, call. = FALSE)
lin <- sf::st_geometry(read_vector_utm(lin_path, "Lineaments"))
message("reading OSM tiles with the direct XML parser...")
rfiles <- list.files(EXT("osm"), "^roads_.*\\.osm$", full.names = TRUE)
pfiles <- list.files(EXT("osm"), "^places_.*\\.osm$", full.names = TRUE)
rd <- read_osm_xml_tiles(rfiles, parse_osm_roads_xml); pl <- read_osm_xml_tiles(pfiles, parse_osm_places_xml)
rd <- sf::st_transform(rd[!is.na(classify_highway(rd$highway)), ], CRS_ANALYSIS)
pl <- sf::st_transform(pl[pl$place %in% PLACE_SETTLE, ], CRS_ANALYSIS)
ctx <- list(terrain = terra::rast(f$t), geomorph = gm, sets = sets, lineament = lin,
            roads_major = sf::st_geometry(rd[rd$highway %in% HWY_MAJOR, ]), roads_any = sf::st_geometry(rd),
            places_all = sf::st_geometry(pl), places_town = sf::st_geometry(pl[pl$place %in% PLACE_TOWN, ]),
            landcover = terra::rast(f$w))

# ---- QA: recompute covariates at random grid cell centres and compare with the stored grid ------------
message("QA: ", QA_N, " random grid cells recomputed with the site code path")
gd  <- grid[grid$in_survey & !grid$in_core, ]
qi  <- sample.int(nrow(gd), min(QA_N, nrow(gd)))
qpts <- sf::st_as_sf(sf::st_drop_geometry(gd)[qi, c("x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
cq   <- covariates_at(qpts, ctx)
qa   <- qa_against_grid(cq, sf::st_drop_geometry(gd)[qi, ])
write_csv_safe(qa, P("output", "tables", "S10_extraction_qa.csv"))
if (any(qa$n_mismatch > 0)) {
  print(qa, row.names = FALSE)
  stop("QA FAILED: the site code path does not reproduce the grid values. No site covariates were written.", call. = FALSE)
}
message("QA passed: all covariates reproduce the stored grid values")

# ---- sites ----------------------------------------------------------------------------------------------
spts <- sf::st_as_sf(sites[, c("site_id", "x", "y")], coords = c("x", "y"), crs = CRS_ANALYSIS, remove = FALSE)
cs   <- covariates_at(spts, ctx)
sc   <- cbind(sites, cs)
sc$in_survey <- lengths(sf::st_intersects(spts, dsurv)) > 0
sc$in_core   <- lengths(sf::st_intersects(spts, core)) > 0
sc$in_domain <- sc$in_survey & !sc$in_core
cov_vars <- c("elev_m", "slope_deg", "gm4", paste0("d_k", SIZE_K), "d_lineament", "d_road_major", "d_road_any",
              "d_settle", "d_town", paste0("frac_", names(WC_GROUPS)))
sc$any_na <- !stats::complete.cases(sc[, cov_vars])
for (b in BLOCKS_KM) sc[[paste0("blk", b)]] <- cluster_ids(as.matrix(sc[, c("x", "y")]), b * 1000)
saveRDS(sc, P("data", "processed", "sites_cov.rds"))          # exact coordinates: stays in data/ (git-ignored)

# ---- count-only tables (no coordinates, no covariate values) -----------------------------------------------
sets_def <- c(S1_openair = "S1_openair", S2_all = "S2_all", S3_riverbed = "S3_riverbed", S4_other = "S4_other")
dom_tab <- do.call(rbind, lapply(names(sets_def), function(s) {
  m <- sc[[sets_def[[s]]]] %in% TRUE
  data.frame(sample_set = s, n_total = sum(m), n_outside_survey_domain = sum(m & !sc$in_survey),
             n_in_removed_core = sum(m & sc$in_core), n_in_analysis_domain = sum(m & sc$in_domain),
             n_in_domain_with_na_covariate = sum(m & sc$in_domain & sc$any_na))
}))
write_csv_safe(dom_tab, P("output", "tables", "S10_site_domain.csv"))
na_tab <- data.frame(covariate = cov_vars, n_na_sites_in_domain = vapply(cov_vars, function(v) sum(is.na(sc[[v]]) & sc$in_domain), numeric(1)))
write_csv_safe(na_tab, P("output", "tables", "S10_na_by_covariate.csv"))
dxy <- as.matrix(sf::st_drop_geometry(grid[grid$in_survey & !grid$in_core, c("x", "y")]))
dblk <- stats::setNames(lapply(BLOCKS_KM, function(b) length(unique(cluster_ids(dxy, b * 1000)))), BLOCKS_KM)
st1 <- cluster_structure(as.matrix(sc[sc$S1_openair %in% TRUE & sc$in_domain, c("x", "y")]), BLOCKS_KM, dblk); st1$sample_set <- "S1_openair"
st2 <- cluster_structure(as.matrix(sc[sc$in_domain, c("x", "y")]), BLOCKS_KM, dblk); st2$sample_set <- "S2_all"
clu <- rbind(st1, st2)
write_csv_safe(clu, P("output", "tables", "S10_cluster_structure.csv"))
s1 <- sc[sc$S1_openair %in% TRUE & sc$in_domain, ]
reft <- sort(table(s1$ref_clean), decreasing = TRUE)
ref_tab <- data.frame(sample_set = "S1_openair", n_references = length(reft), top1_share = round(reft[[1]] / sum(reft), 3),
                      top3_share = round(sum(reft[1:min(3, length(reft))]) / sum(reft), 3), singleton_references = sum(reft == 1L))
write_csv_safe(ref_tab, P("output", "tables", "S10_reference_structure.csv"))
write_csv_safe(file_checksums(f$sites), P("output", "logs", "S10_input_checksums.csv"))
write_csv_safe(data.frame(item = c("BLOCKS_KM", "QA_N", "QA_tolerance_m", "QA_tolerance_fraction", "n_sites", "blind_rule"),
                          value = c(paste(BLOCKS_KM, collapse = ","), QA_N, "1e-6 (distances); 1e-3 (elevation, slope)", 1e-6, nrow(sc),
                                    "no covariate value of any site is printed or tabulated in S10")),
               P("output", "logs", "S10_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S10_sessionInfo.txt"))

cat("\n========== S10 SUMMARY ==========\n")
cat("QA vs grid (", length(qi), " random cells): all covariates reproduce stored values; max abs diff per variable:\n", sep = "")
print(qa[, c("variable", "n_compared", "max_abs_diff", "n_mismatch")], row.names = FALSE)
cat("\nSites by sample set and analysis-domain membership:\n"); print(dom_tab, row.names = FALSE)
cat("\nNA count by covariate (sites in domain):\n"); print(na_tab[na_tab$n_na_sites_in_domain > 0, ], row.names = FALSE)
if (!any(na_tab$n_na_sites_in_domain > 0)) cat("  none\n")
cat("\nCluster structure (counts only):\n"); print(clu, row.names = FALSE)
cat("\nReference structure of S1:\n"); print(ref_tab, row.names = FALSE)
cat("\nBLIND: no covariate value of any site was printed or tabulated. data/processed/sites_cov.rds holds exact coordinates (git-ignored).\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S10. Paste this console output back.\n")