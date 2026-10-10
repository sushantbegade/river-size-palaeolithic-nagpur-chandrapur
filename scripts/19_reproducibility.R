# scripts/19_reproducibility.R   (S19)
# Checks that the project is rebuildable and safe to publish. Reads saved outputs; never touches raw data.
#  RUN_MODE "smoke" (default): environment record, design-freeze check, data-safety scan, tests, S18 run twice in
#                              fresh R sessions (outputs must be identical), comparison with the stored baseline.
#  RUN_MODE "full": additionally re-runs S01-S18 in order, each in a fresh R session (hours; recomputes outputs).
if (!exists("RUN_MODE")) RUN_MODE <- "smoke"
if (!exists("RUN_DOWNLOADS")) RUN_DOWNLOADS <- FALSE        # S04 downloads external data; only with TRUE
suppressPackageStartupMessages(library(here))
source(here::here("R", "config.R")); source(here::here("R", "utils.R")); source(here::here("R", "repro.R"))
t0 <- Sys.time(); root <- normalizePath(here::here(), winslash = "/"); logd <- P("output", "logs")
dir.create(logd, recursive = TRUE, showWarnings = FALSE)
message("== S19 reproducibility check, mode: ", RUN_MODE, " ==")
res <- list()   # check name -> "PASS" / "WARN" / "FAIL"

# ---- 1. environment record ---------------------------------------------------------------------------------------
pk <- repro_pkgs(root); pv <- repro_pkg_versions(pk)
git_v <- tryCatch(system2("git", "--version", stdout = TRUE), error = function(e) "git not found")
env <- c(R.version.string, paste("platform:", R.version$platform), paste("running:", paste(utils::sessionInfo()$running, collapse = " ")),
         paste("project root contains '&' or space:", grepl("[& ]", root)), paste("PALEO_DATA set:", nzchar(Sys.getenv("PALEO_DATA"))),
         paste("renv.lock present:", file.exists(file.path(root, "renv.lock"))), paste("renv loaded:", "renv" %in% loadedNamespaces()), git_v[1], "",
         "Packages used by project code (scan of library(), require(), pkg::):", utils::capture.output(print(pv, row.names = FALSE)))
writeLines(env, file.path(logd, "S19_environment.txt"))
res$packages <- if (any(!pv$installed)) "FAIL" else "PASS"
res$path <- if (grepl("[& ]", root)) "WARN" else "PASS"
res$renv <- if (file.exists(file.path(root, "renv.lock"))) "PASS" else "WARN"

# ---- 2. design freeze --------------------------------------------------------------------------------------------
fz <- repro_verify_freeze(root)
utils::write.csv(fz$table, file.path(logd, "S19_freeze_check.csv"), row.names = FALSE)
res$freeze_manifest <- if (fz$md5_ok) "PASS" else "FAIL"
chg <- fz$table$file[fz$table$status == "changed"]; mis <- fz$table$file[fz$table$status == "missing"]
design_chg <- chg[!grepl("^output/", chg)]

# ---- 3. data safety ----------------------------------------------------------------------------------------------
scan_files <- unlist(lapply(c("R", "scripts", "tests", "docs", "output/tables", "output/logs"), function(d) list.files(file.path(root, d), recursive = TRUE, full.names = TRUE)))
hits <- repro_scan_coords(scan_files, root)
res$coord_scan <- if (length(hits)) "FAIL" else "PASS"
pub <- repro_public_grid(file.path(root, REPRO_PUBLIC)); res$public_grid <- if (isTRUE(pub$ok)) "PASS" else if (is.na(pub$ok)) "WARN" else "FAIL"
gi <- repro_gitignore_check(root); res$gitignore <- if (length(gi)) "WARN" else "PASS"

# ---- 4. run ------------------------------------------------------------------------------------------------------
paper <- P("output", "tables", "paper"); figs <- P("output", "figures"); s18 <- P("scripts", "18_figures_tables.R")
snap <- function() c(repro_md5_set(paper, "\\.(csv|md)$", exclude = "S18_outputs_manifest.csv"), repro_md5_set(figs, "\\.png$"))
runlog <- NULL; tests <- c(passed = NA, failed = NA, skipped = NA, errors = NA)
if (RUN_MODE == "full") {
  before <- repro_md5_set(P("output", "tables"), "\\.csv$")
  for (s in REPRO_ORDER) {
    if (s == "04_download_external.R" && !isTRUE(RUN_DOWNLOADS)) { runlog <- rbind(runlog, data.frame(script = s, status = "skipped (RUN_DOWNLOADS = FALSE)", seconds = 0, log = "")); next }
    message("running ", s, " ..."); r <- repro_run_script(P("scripts", s), root, logd); runlog <- rbind(runlog, r)
    if (r$status != "ok") { message("STOP: ", s, " failed; see output/logs/", r$log); break }
  }
  after <- repro_md5_set(P("output", "tables"), "\\.csv$")
  common <- intersect(names(before), names(after)); differ <- common[before[common] != after[common]]
  res$full_rebuild <- if (all(runlog$status %in% c("ok") | grepl("^skipped", runlog$status)) && nrow(runlog) == length(REPRO_ORDER)) "PASS" else "FAIL"
  res$full_outputs_identical <- if (length(differ)) "WARN" else "PASS"
}
tests <- repro_run_tests(root, logd); res$tests <- if (anyNA(tests)) "FAIL" else if (tests[["failed"]] + tests[["errors"]] > 0) "FAIL" else "PASS"
r1 <- repro_run_script(s18, root, logd, "_run1"); s1 <- snap()
r2 <- repro_run_script(s18, root, logd, "_run2"); s2 <- snap()
runlog <- rbind(runlog, r1, r2)
same12 <- identical(s1, s2); diff12 <- names(s2)[s1[names(s2)] != s2 | is.na(s1[names(s2)])]
res$s18_runs <- if (r1$status == "ok" && r2$status == "ok") "PASS" else "FAIL"
res$s18_deterministic <- if (same12) "PASS" else "FAIL"
bf <- file.path(logd, "S19_baseline_md5.csv"); base_msg <- ""
if (!file.exists(bf)) { utils::write.csv(data.frame(file = names(s2), md5 = unname(s2)), bf, row.names = FALSE); base_msg <- "baseline CREATED from this run (first S19 run)"; res$baseline <- "PASS"
} else { b <- utils::read.csv(bf, stringsAsFactors = FALSE); bb <- setNames(b$md5, b$file); com <- intersect(names(bb), names(s2))
dif <- c(com[bb[com] != s2[com]], setdiff(names(bb), names(s2)), setdiff(names(s2), names(bb)))
base_msg <- if (length(dif)) paste0("baseline DIFFERS in ", length(dif), " file(s): ", paste(utils::head(dif, 8), collapse = ", ")) else "baseline identical"
res$baseline <- if (length(dif)) "WARN" else "PASS" }
utils::write.csv(runlog, file.path(logd, "S19_run_log.csv"), row.names = FALSE)

# ---- 5. console summary ------------------------------------------------------------------------------------------
cat("\n========== S19 SUMMARY ==========\n")
cat(sprintf("Mode: %s | R: %s | project root has '&' or space: %s\n", RUN_MODE, R.version.string, grepl("[& ]", root)))
cat(sprintf("Packages used by code: %d | not installed: %s\n", nrow(pv), if (any(!pv$installed)) paste(pv$package[!pv$installed], collapse = ", ") else "none"))
cat(sprintf("Design freeze manifest MD5 matches the registered value: %s (%s)\n", fz$md5_ok, fz$manifest_md5))
cat(sprintf("Frozen files: %d unchanged, %d changed, %d missing | design files changed after freeze: %d | new design files: %d\n",
            sum(fz$table$status == "unchanged"), length(chg), length(mis), length(design_chg), length(fz$new_design_files)))
if (length(design_chg)) cat("  changed design files (each must be explained in docs/deviation_log*.md):", paste(design_chg, collapse = ", "), "\n")
if (length(mis)) cat("  MISSING frozen files:", paste(mis, collapse = ", "), "\n")
cat("  new design files since freeze:", if (length(fz$new_design_files)) paste(fz$new_design_files, collapse = ", ") else "none", "\n")
cat(sprintf("Coordinate scan (%d files): %s\n", length(scan_files), if (length(hits)) paste0(length(hits), " file(s) flagged (file:line only) -> ", paste(hits, collapse = " | ")) else "no coordinate pairs found"))
cat(sprintf("Public 1 km file: %s\n", pub$msg))
cat(sprintf(".gitignore lacks: %s\n", if (length(gi)) paste(gi, collapse = " ") else "nothing required missing"))
cat(sprintf("Tests (fresh session): passed=%s failed=%s skipped=%s errors=%s\n", tests[["passed"]], tests[["failed"]], tests[["skipped"]], tests[["errors"]]))
cat(sprintf("S18 run 1: %s (%.1f s) | run 2: %s (%.1f s) | outputs identical between runs: %s", r1$status, r1$seconds, r2$status, r2$seconds, same12))
if (!same12) cat(" | differing:", paste(utils::head(diff12, 8), collapse = ", ")); cat("\n")
cat("Baseline:", base_msg, "\n")
if (RUN_MODE == "full") { cat("Full rebuild:\n"); print(runlog[seq_len(nrow(runlog) - 2), ], row.names = FALSE) }
cat("\nCHECKS:\n"); for (n in names(res)) cat(sprintf("  %-24s %s\n", n, res[[n]]))
cat(sprintf("\nOVERALL: %s\n", if (any(unlist(res) == "FAIL")) "FAIL" else if (any(unlist(res) == "WARN")) "PASS WITH WARNINGS" else "PASS"))
cat("Files: output/logs/S19_environment.txt, S19_freeze_check.csv, S19_run_log.csv, S19_baseline_md5.csv, S19_run_*.log\n")
cat("Elapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("DONE S19. Paste this console output back.\n")