# scripts/09e_stepzero_profile.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))
source(here::here("R", "stepzero_e.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S09e Step 0: how well does the AIC profile recover the true size threshold k? (NO real site is used) ==")

SIM_REPS <- as.integer(Sys.getenv("S09_REPS", "80"))
KS       <- SIZE_K[SIZE_K >= 3L]          # 3,4,5,6,7,8,9
K_TRUES  <- c(5L, 6L, 7L, 8L, 9L)
RATIOS   <- c(2, 3)
NS       <- c(93L, 48L)
SETTINGS <- list(A = list(sigma = 0.7, range_km = 10), B = list(sigma = 1.0, range_km = 30))   # B = stress

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
X  <- nuisance_design_c(dom)
zs <- stats::setNames(lapply(SIZE_K, function(k) z_size(dom, k)), SIZE_K)

n_fits <- SIM_REPS * length(SETTINGS) * length(NS) * length(K_TRUES) * length(RATIOS) * length(KS)
message("simulation: ", n_fits, " ridge fits (", SIM_REPS, " reps x ", length(SETTINGS), " settings x ", length(NS), " n x ",
        length(K_TRUES), " true k x ", length(RATIOS), " ratios x ", length(KS), " fitted k)")
tt <- Sys.time()
sums <- list(); conf <- list()
for (nm in names(SETTINGS)) {
  st <- SETTINGS[[nm]]
  message("setting ", nm, ": sigma ", st$sigma, ", range ", st$range_km, " km")
  res <- run_profile(X, zs, a, xy, KS, K_TRUES, RATIOS, NS, st$sigma, st$range_km * 1000, SIM_REPS)
  s <- summarise_profile(res); s$setting <- nm; s$sigma <- st$sigma; s$range_km <- st$range_km
  sums[[nm]] <- s
  for (n_t in NS) for (ra in RATIOS) {
    cf <- confusion_profile(res, n_t, ra, KS, K_TRUES)
    conf[[paste(nm, n_t, ra)]] <- list(setting = nm, n = n_t, ratio = ra, table = cf)
  }
}
message("simulation minutes: ", round(as.numeric(difftime(Sys.time(), tt, units = "mins")), 2))
summ <- do.call(rbind, sums)
write_csv_safe(summ, P("output", "tables", "S09e_profile_summary.csv"))
conf_flat <- do.call(rbind, lapply(conf, function(z) cbind(setting = z$setting, n_target = z$n, ratio = z$ratio,
                                                           k_true = rownames(z$table), z$table)))
write_csv_safe(conf_flat, P("output", "tables", "S09e_profile_confusion.csv"))
write_csv_safe(data.frame(item = c("SIM_REPS", "KS", "K_TRUES", "RATIOS", "NS", "settings", "PRIOR_SD", "k_classes", "seed"),
                          value = c(SIM_REPS, paste(KS, collapse = ","), paste(K_TRUES, collapse = ","), paste(RATIOS, collapse = ","),
                                    paste(NS, collapse = ","),
                                    paste(vapply(names(SETTINGS), function(n) paste0(n, ": sigma ", SETTINGS[[n]]$sigma, ", range ", SETTINGS[[n]]$range_km, " km"), character(1)), collapse = "; "),
                                    PRIOR_SD, "low = k<=5, mid = 6-7, high = 8-9", SEED)),
               P("output", "logs", "S09e_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S09e_sessionInfo.txt"))

cat("\n========== S09e SUMMARY ==========\n")
cat("Profile recovery: p_exact / p_within1 / p_class (low k<=5, mid 6-7, high 8-9); median and 90th percentile of AIC(true k) - min AIC:\n")
print(summ[order(summ$setting, summ$n_target, summ$ratio, summ$k_true), c("setting", "sigma", "range_km", "n_target", "ratio", "k_true",
                                                                          "n_ok", "p_exact", "p_within1", "p_class", "mean_abs_err",
                                                                          "median_daic_true", "p90_daic_true")], row.names = FALSE)
cat("\nPooled over true k (mean):\n")
print(stats::aggregate(cbind(p_exact, p_within1, p_class) ~ setting + n_target + ratio, summ, function(x) round(mean(x), 3)), row.names = FALSE)
cat("\nConfusion matrices (rows = true k, columns = k chosen by min AIC), ratio 2, n = 93:\n")
for (nm in names(SETTINGS)) { cat("\nsetting", nm, "\n"); print(conf[[paste(nm, 93L, 2)]]$table) }
cat("\nNOTE: truth values (nuisance, field SD, range, effect size) are assumptions. Adjacent k are correlated (S06: r up to 0.80), so exact k is expected to be hard.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S09e. Paste this console output back.\n")