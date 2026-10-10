# scripts/check_repo_safe.R   (S19 version; REPLACES the earlier one)
# Run AFTER `git add .` and BEFORE `git commit`. Stops if anything risky is staged:
# data/spatial file types, secrets, large files, anything under data/ except the 1 km public file,
# CSV columns that look like coordinates, and exact coordinate pairs hidden in any text file.
suppressPackageStartupMessages(library(here))
source(here::here("R", "repro.R"))
r <- repro_check_repo(here::here())
if (length(r$notes)) cat(paste0(r$notes, "\n"), sep = "")
if (length(r$problems)) {
  cat("NOT SAFE TO COMMIT:\n"); cat(paste0(" - ", r$problems), sep = "\n")
  stop("Fix the issues above (unstage with: git rm --cached <file>).", call. = FALSE)
}
cat("OK: ", r$n_files, " tracked files, none flagged.\n", sep = "")