# scripts/run_all.R   (S19 version; REPLACES the earlier stub)
# Entry point for reproducibility. Default is the quick "smoke" check (minutes).
# Full rebuild of every analysis step in fresh R sessions (hours; recomputes outputs):
#   RUN_MODE <- "full"; source(here::here("scripts", "run_all.R"))
# S04 downloads external data and runs only if RUN_DOWNLOADS <- TRUE is set first.
# Order of scripts: REPRO_ORDER in R/repro.R.
if (!exists("RUN_MODE")) RUN_MODE <- "smoke"
source(here::here("scripts", "19_reproducibility.R"))