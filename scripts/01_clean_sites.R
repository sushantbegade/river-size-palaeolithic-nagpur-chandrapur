suppressPackageStartupMessages({
  library(here); library(readxl); library(readr); library(dplyr)
  library(tidyr); library(stringr); library(sf)
})
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "clean_sites.R"))

set.seed(SEED)
message("== S01 clean sites ==")

xlsx      <- site_file_path()
crosswalk <- readr::read_csv(P("docs", "crosswalk_period.csv"), show_col_types = FALSE)

raw   <- read_sites_raw(xlsx)
sites <- clean_sites(raw, crosswalk)
sites <- add_utm(sites)
validate_sites(sites)
message("rows: ", nrow(sites), " | unique ids: ", dplyr::n_distinct(sites$site_id))

dir.create(P("data", "processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(sites, P("data", "processed", "sites_clean.rds"))
write_csv_safe(sites, P("data", "processed", "sites_clean.csv"))
write_csv_safe(coarsen_sites(sites, 1000), P("data", "public", "sites_public_1km.csv"))

tabs <- summarise_sites(sites)
for (nm in names(tabs)) write_csv_safe(tabs[[nm]], P("output", "tables", paste0(nm, ".csv")))

thin <- thinning_table(sites, dists_m = c(2000, 5000), n_rep = 50L, seed = SEED)
write_csv_safe(thin, P("output", "tables", "S01_thinning_counts.csv"))

write_csv_safe(file_checksums(xlsx), P("output", "logs", "S01_input_checksums.csv"))
writeLines(capture.output(sessionInfo()), P("output", "logs", "S01_sessionInfo.txt"))

cat("\n========== S01 SUMMARY ==========\n")
cat("N sites:", nrow(sites), "\n")
cat("\nSite type:\n");        print(table(sites$site_type))
cat("\nCoordinate tier:\n");  print(table(sites$coord_tier))
cat("\nPeriod class (single-period vs multi):\n"); print(table(sites$period_class))
cat("\nPeriod membership (overlapping): LP", sum(sites$LP), "| MP", sum(sites$MP),
    "| UP", sum(sites$UP), "\n")
cat("\nSets: S1_openair", sum(sites$S1_openair), "| S2_all", nrow(sites),
    "| S3_riverbed", sum(sites$S3_riverbed), "| S4_other", sum(sites$S4_other), "\n")
cat("\nMegalith/cairn flagged:", sum(sites$flag_megalith), "\n")
cat("\nDistinct reference strings:", dplyr::n_distinct(sites$ref_clean), "\n")
cat("\nThinning counts (median of 50 random orders):\n"); print(as.data.frame(thin), row.names = FALSE)
cat("\nDONE S01. Files written. Paste this console output back.\n")