# scripts/09d_stepzero_block_rule.R
suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "quadrature.R"))
source(here::here("R", "stepzero.R"))
source(here::here("R", "stepzero_b.R"))
source(here::here("R", "stepzero_c.R"))

options(max.print = 1e5)
set.seed(SEED)
t0 <- Sys.time()
message("== S09d Step 0: is the cluster block-size rule robust to the (unknown) spatial range? (NO real site is used) ==")

SIM_REPS   <- as.integer(Sys.getenv("S09_REPS", "200"))
RANGES_KM  <- c(5, 10, 30, 60)     # correlation range of the simulated clustering field
SIGMAS     <- c(0.7, 1.0)          # field SD (moderate, strong)
NS         <- c(93L, 48L)
KS         <- c(6L, 8L)
RATIOS     <- c(1, 2)
BLOCKS_KM  <- c(10, 20, 40)
BOOT_B     <- 300L

# block-size rule, proposed BEFORE seeing S09d results (goes to the deviation log):
#   pick the SMALLEST block size with  mean type I <= 0.06, worst-cell type I <= 0.09,
#   mean coverage >= 0.93 and worst-cell coverage >= 0.90 over all cells; if none qualifies use the LARGEST.
RULE <- c(t1_mean = 0.06, t1_worst = 0.09, cov_mean = 0.93, cov_worst = 0.90)

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

n_fits <- SIM_REPS * length(RANGES_KM) * length(SIGMAS) * length(NS) * length(KS) * length(RATIOS)
message("simulation: ", n_fits, " ridge fits, each with ", length(BLOCKS_KM), " bootstraps (B = ", BOOT_B, ")")
tt <- Sys.time()
ps_all <- list()
for (rg in RANGES_KM) {
  message("field range ", rg, " km")
  res <- run_power_c(X, zs, a, xy, KS, NS, RATIOS, SIGMAS, BLOCKS_KM, SIM_REPS, rg * 1000, B = BOOT_B, progress = 25L)
  ps <- summarise_power_c(res); ps$range_km <- rg
  ps_all[[as.character(rg)]] <- ps
}
ps <- do.call(rbind, ps_all)
message("simulation minutes: ", round(as.numeric(difftime(Sys.time(), tt, units = "mins")), 2))
write_csv_safe(ps, P("output", "tables", "S09d_power_summary.csv"))

null <- ps[ps$ratio == 1, ]; alt <- ps[ps$ratio > 1, ]
keys <- unique(ps[, c("range_km", "sigma", "n_target", "block_km")])
cells <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  s <- function(d) d$range_km == keys$range_km[i] & d$sigma == keys$sigma[i] &
    d$n_target == keys$n_target[i] & d$block_km == keys$block_km[i]
  data.frame(keys[i, ],
             type1 = round(mean(null$rej_uncond[s(null)]), 3),
             coverage = round(mean(c(null$coverage[s(null)], alt$coverage[s(alt)])), 3),
             power_ratio2 = round(mean(alt$rej_uncond[s(alt)]), 3),
             se_boot_over_sd = round(mean(c(null$se_boot_over_sd[s(null)], alt$se_boot_over_sd[s(alt)])), 2),
             share_failed = round(mean(ps$share_failed[s(ps)]), 3))
}))
write_csv_safe(cells, P("output", "tables", "S09d_cells.csv"))

rule <- do.call(rbind, lapply(BLOCKS_KM, function(b) {
  d <- cells[cells$block_km == b, ]
  data.frame(block_km = b, t1_mean = round(mean(d$type1), 3), t1_worst = max(d$type1),
             cov_mean = round(mean(d$coverage), 3), cov_worst = min(d$coverage),
             power_ratio2_mean = round(mean(d$power_ratio2), 3),
             power_ratio2_n93 = round(mean(d$power_ratio2[d$n_target == 93]), 3),
             fail_mean = round(mean(d$share_failed), 3))
}))
rule$passes <- rule$t1_mean <= RULE["t1_mean"] & rule$t1_worst <= RULE["t1_worst"] &
  rule$cov_mean >= RULE["cov_mean"] & rule$cov_worst >= RULE["cov_worst"]
chosen <- if (any(rule$passes)) min(rule$block_km[rule$passes]) else max(rule$block_km)
write_csv_safe(rule, P("output", "tables", "S09d_block_rule.csv"))
write_csv_safe(data.frame(item = c("SIM_REPS", "RANGES_KM", "SIGMAS", "NS", "KS", "RATIOS", "BLOCKS_KM", "PRIOR_SD", "BOOT_B",
                                   names(RULE), "chosen_block_km", "seed"),
                          value = c(SIM_REPS, paste(RANGES_KM, collapse = ","), paste(SIGMAS, collapse = ","), paste(NS, collapse = ","),
                                    paste(KS, collapse = ","), paste(RATIOS, collapse = ","), paste(BLOCKS_KM, collapse = ","),
                                    PRIOR_SD, BOOT_B, RULE, chosen, SEED)),
               P("output", "logs", "S09d_parameters.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S09d_sessionInfo.txt"))

cat("\n========== S09d SUMMARY ==========\n")
cat("Blocks in domain (10/20/40 km): ", paste(vapply(BLOCKS_KM, function(b) length(unique(cluster_ids(xy, b * 1000))), integer(1)), collapse = " / "), "\n", sep = "")
cat("\nPer cell (type1 = null rejection pooled over k; coverage pooled; power at ratio 2; se_boot_over_sd = median bootstrap SE / true SD):\n")
print(cells[order(cells$block_km, cells$range_km, cells$sigma, cells$n_target), ], row.names = FALSE)
cat("\nBlock-size rule (mean type I <= ", RULE["t1_mean"], ", worst <= ", RULE["t1_worst"],
    "; mean coverage >= ", RULE["cov_mean"], ", worst >= ", RULE["cov_worst"], "):\n", sep = "")
print(rule, row.names = FALSE)
cat("\nBlock size selected by the rule:", chosen, "km\n")
cat("\nNOTE: simulated clusters are spatial blocks only; the real analysis uses reference x block. All truth values are assumptions.\n")
cat("\nElapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S09d. Paste this console output back.\n")