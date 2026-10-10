# scripts/freeze_manifest.R -- record checksums of the frozen design (no Git needed)
# Run AFTER docs/deviation_log.md and docs/preregistration.md are final and BEFORE script S10.
suppressPackageStartupMessages(library(here))
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

ls_rel <- function(dir, pattern) {
  f <- list.files(P(dir), pattern = pattern, recursive = TRUE, full.names = TRUE)
  f[!grepl("/(data|renv|\\.git)/", f)]
}
files <- c(ls_rel("R", "\\.R$"), ls_rel("scripts", "\\.R$"), ls_rel("tests", "\\.R$"),
           ls_rel("docs", "\\.(md|csv|txt)$"), ls_rel("output/tables", "\\.csv$"))
files <- sort(unique(files))
if (!any(grepl("deviation_log\\.md$", files)) || !any(grepl("preregistration\\.md$", files))) {
  stop("docs/deviation_log.md and docs/preregistration.md must exist before freezing.", call. = FALSE)
}
rel <- sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", here::here()), "/"), "", files)
man <- data.frame(file = rel, md5 = unname(tools::md5sum(files)), bytes = file.size(files),
                  modified = format(file.mtime(files), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE)
man$frozen_on <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
out <- P("output", "logs", "DESIGN_FREEZE_manifest.csv")
write_csv_safe(man, out)
cat("\n========== FREEZE MANIFEST ==========\n")
cat("Files recorded:", nrow(man), "\n")
cat("By folder:\n"); print(table(sub("/.*$", "", man$file)))
cat("\nMD5 of the manifest file (paste this into the OSF registration):", unname(tools::md5sum(out)), "\n")
cat("Frozen on:", man$frozen_on[1], "\n")
cat("\nCheck: no data, raster or coordinate files are listed above.\n")
cat("DONE freeze_manifest. Paste this console output back.\n")