# scripts/00_setup_project.R -- run ONCE after opening the .Rproj and running renv::init()
# Checks the environment and creates folders. Reads no data.

stopifnot("Open waterbody-palaeolithic.Rproj first" = file.exists("waterbody-palaeolithic.Rproj"))
cat("R version:", R.version.string, "\n")
if (getRversion() < "4.1.0") stop("R >= 4.1 required (native pipe).")

core <- c("here", "readxl", "readr", "dplyr", "tidyr", "stringr", "purrr",
          "tibble", "sf", "testthat")
missing <- core[!vapply(core, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  if (requireNamespace("renv", quietly = TRUE) && file.exists("renv.lock")) {
    renv::install(missing)
  } else install.packages(missing)
}

for (d in c("data/raw", "data/processed", "data/public", "data/external",
            "output/tables", "output/logs", "output/figures", "docs"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

root <- Sys.getenv("PALEO_DATA", "")
cat("PALEO_DATA =", if (nzchar(root)) root else "<NOT SET>", "\n")
if (!nzchar(root)) {
  cat("-> copy .Renviron.template to .Renviron, edit the path, then Session > Restart R.\n")
} else if (!dir.exists(root)) {
  cat("-> folder does not exist. Fix the path in .Renviron.\n")
} else cat("data folder found.\n")

cat("\nPackage versions:\n")
for (p in core) cat(sprintf("  %-10s %s\n", p, as.character(packageVersion(p))))
if (requireNamespace("renv", quietly = TRUE) && file.exists("renv.lock")) {
  renv::snapshot(prompt = FALSE)
  cat("renv.lock updated.\n")
} else cat("renv.lock not found: run renv::init() in the console first.\n")
cat("\nSetup check finished.\n")
