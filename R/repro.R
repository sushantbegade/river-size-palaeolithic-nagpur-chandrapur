# R/repro.R  (S19 reproducibility helpers). No site coordinates are ever printed or written.
REPRO_FREEZE_MD5 <- "88ceb20829091aa0328217c460341466"     # MD5 of output/logs/DESIGN_FREEZE_manifest.csv (design freeze, 4 Oct 2026)
REPRO_ORDER <- c("01_clean_sites.R", "02_prepare_hydrography.R", "03_persistence_index.R", "04_download_external.R", "05_regulation_flags.R",
                 "05b_cumulative_regulation.R", "06_quadrature_predictors.R", "07_terrain_geomorph.R", "08_detectability.R", "08b_osm_diagnostic.R",
                 "09a_stepzero_engine.R", "09b_stepzero_power.R", "09c_stepzero_ridge_blocks.R", "09d_stepzero_block_rule.R", "09e_stepzero_profile.R",
                 "10_site_covariates.R", "11_primary_model.R", "12_secondary_diagnostics.R", "13_site_type_composition.R", "14_spatial_cv.R",
                 "15_survey_bias.R", "16_robustness.R", "17_exploratory.R", "18_figures_tables.R")
REPRO_PUBLIC <- "data/public/sites_public_1km.csv"
REPRO_SCAN_SKIP <- c(REPRO_PUBLIC, "output/tables/S04_gdw_barriers_within100km.csv")   # third-party dam positions, git-ignored
REPRO_COORD_NAMES <- "^(x|y|x_?c|y_?c|x_?utm|y_?utm|lon|lat|long|longitude|latitude|easting|northing|coord.*|utm.*)$"
REPRO_TEXT_EXT <- "\\.(csv|tsv|txt|md|r|rmd|json|ya?ml|log)$"

repro_rscript <- function() file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# ---- packages used by the project code (for the environment record / renv) --------------------------------------
repro_pkgs <- function(root) {
  f <- unlist(lapply(c("R", "scripts", "tests"), function(d) list.files(file.path(root, d), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)))
  pk <- character(0)
  for (x in f) {
    l <- sub("#.*$", "", readLines(x, warn = FALSE))
    pk <- c(pk, unlist(regmatches(l, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=::)", l, perl = TRUE))))
    for (fn in c("library", "require", "requireNamespace")) {
      m <- unlist(regmatches(l, gregexpr(paste0("(?<=", fn, "\\()\\s*[\"']?[A-Za-z][A-Za-z0-9.]*"), l, perl = TRUE)))
      pk <- c(pk, trimws(gsub("[\"']", "", m)))
    }
  }
  base <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "grid", "parallel", "datasets", "compiler", "splines", "stats4")
  sort(setdiff(unique(pk), c(base, "", "pkg")))   # "pkg" = placeholder in comments/strings
}
repro_pkg_versions <- function(pk) {
  data.frame(package = pk, installed = nzchar(vapply(pk, function(p) system.file(package = p), "")),
             version = vapply(pk, function(p) if (nzchar(system.file(package = p))) as.character(utils::packageVersion(p)) else NA_character_, ""),
             row.names = NULL, stringsAsFactors = FALSE)
}

# ---- design freeze verification ----------------------------------------------------------------------------------
repro_verify_freeze <- function(root, manifest = file.path(root, "output", "logs", "DESIGN_FREEZE_manifest.csv"), expect_md5 = REPRO_FREEZE_MD5) {
  if (!file.exists(manifest)) stop("Freeze manifest not found: ", manifest, call. = FALSE)
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE); cur <- file.path(root, m$file); ex <- file.exists(cur)
  now <- rep(NA_character_, nrow(m)); now[ex] <- unname(tools::md5sum(cur[ex]))
  m$status <- ifelse(!ex, "missing", ifelse(now == m$md5, "unchanged", "changed"))
  m$status <- repro_mark_documented(root, m, now)
  design <- unlist(lapply(c("R", "scripts", "tests", "docs"), function(d)
    paste0(d, "/", list.files(file.path(root, d), pattern = "\\.(R|md|csv|txt)$", recursive = TRUE))))
  list(manifest_md5 = unname(tools::md5sum(manifest)), md5_ok = identical(unname(tools::md5sum(manifest)), expect_md5),
       table = m[, c("file", "status")], new_design_files = sort(setdiff(design, m$file)))
}

# ---- data-safety scans (report file and line number only; never print the matching text) -----------------------------
repro_public_grid <- function(path) {
  if (!file.exists(path)) return(list(ok = NA, msg = "public file absent"))
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("x_c", "y_c") %in% names(d))) return(list(ok = FALSE, msg = "x_c / y_c columns missing"))
  ok <- all(abs((d$x_c %% 1000) - 500) < 1e-6) && all(abs((d$y_c %% 1000) - 500) < 1e-6)
  list(ok = ok, msg = if (ok) sprintf("%d rows, all on 1 km cell centres", nrow(d)) else "NOT on the 1 km cell-centre grid")
}
repro_scan_coords <- function(files, root = "", allow = REPRO_PUBLIC) {
  hits <- character(0)
  utm <- "\\b[1-5][0-9]{5}(\\.[0-9]+)?\\b[^0-9\\n]{1,8}\\b2[0-9]{6}(\\.[0-9]+)?\\b"
  lonp <- "\\b(7[89]|80)\\.[0-9]{2}(?=[0-9]*[1-9])[0-9]{2,}\\b"; latp <- "\\b(19|20|21)\\.[0-9]{2}(?=[0-9]*[1-9])[0-9]{2,}\\b"   # needs a non-zero digit after the 2nd decimal: round fixture numbers are not data
  for (f in files) {
    rel <- if (nzchar(root)) sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", normalizePath(root, winslash = "/", mustWork = FALSE)), "/"), "", normalizePath(f, winslash = "/", mustWork = FALSE)) else f
    if (rel %in% allow || !file.exists(f) || !grepl(REPRO_TEXT_EXT, f, ignore.case = TRUE) || file.size(f) > 2e7) next
    l <- readLines(f, warn = FALSE, encoding = "UTF-8"); l <- iconv(l, "UTF-8", "UTF-8", sub = "")
    i <- which(grepl(utm, l, perl = TRUE) | (grepl(lonp, l, perl = TRUE) & grepl(latp, l, perl = TRUE)))
    if (length(i)) hits <- c(hits, sprintf("%s: line(s) %s", rel, paste(utils::head(i, 5), collapse = ",")))
  }
  hits
}
repro_gitignore_check <- function(root) {
  need <- c("data/", "*.tif", "*.shp", "*.gpkg", "*.rds", "*.xlsx", ".Renviron", ".Rhistory", ".RData", "renv/library/")
  f <- file.path(root, ".gitignore"); have <- if (file.exists(f)) trimws(readLines(f, warn = FALSE)) else character(0)
  norm <- function(z) sub("^/", "", z)
  need[!(norm(need) %in% norm(have) | norm(sub("/$", "", need)) %in% norm(have))]
}

# ---- repository check (used by scripts/check_repo_safe.R; works on `git ls-files`) -------------------------------------
repro_check_repo <- function(root) {
  old <- setwd(root); on.exit(setwd(old))
  files <- system2("git", c("-c", "core.quotepath=off", "ls-files"), stdout = TRUE)
  if (!length(files)) stop("No tracked files found. Run `git add .` first.", call. = FALSE)
  bad_ext <- "\\.(tif|tiff|shp|shx|dbf|prj|cpg|gpkg|geojson|zip|xlsx|xls|rds|rda|rdata)$"
  bad_name <- "(^|/)(\\.Renviron|\\.RData|\\.Rhistory)$|(^|/)renv/library/"
  pr <- character(0); notes <- character(0)
  for (f in files) {
    if (!file.exists(f)) next
    if (grepl(bad_ext, f, ignore.case = TRUE)) pr <- c(pr, paste("forbidden type:", f))
    if (grepl(bad_name, f)) pr <- c(pr, paste("must not be tracked:", f))
    if (file.size(f) > 5e6) pr <- c(pr, sprintf("large file (%.1f MB): %s", file.size(f) / 1e6, f))
    if (grepl("^data/", f) && f != REPRO_PUBLIC) pr <- c(pr, paste("tracked file under data/ (only the public 1 km file may be tracked):", f))
    if (grepl("\\.csv$", f, ignore.case = TRUE) && f != REPRO_PUBLIC && file.size(f) > 0) {
      hdr <- gsub('"', "", tolower(trimws(strsplit(readLines(f, n = 1, warn = FALSE), ",")[[1]])))
      hit <- hdr[grepl(REPRO_COORD_NAMES, hdr)]
      if (length(hit)) pr <- c(pr, paste0("coordinate-like column(s) [", paste(hit, collapse = ","), "] in ", f))
    }
  }
  if (REPRO_PUBLIC %in% files) {
    g <- repro_public_grid(REPRO_PUBLIC)
    if (!isTRUE(g$ok)) pr <- c(pr, paste("public coordinate file is not verified on the 1 km grid:", g$msg))
    notes <- c(notes, "NOTE: the 1 km public coordinate file is tracked. The coordinate-sharing policy (ethics/heritage) must be decided before any push.")
  }
  sc <- repro_scan_coords(file.path(root, files), root)
  if (length(sc)) pr <- c(pr, paste("possible exact coordinates (UTM pair or lon/lat pair) in", sc))
  list(n_files = length(files), problems = pr, notes = notes)
}

# ---- running scripts in fresh R sessions -------------------------------------------------------------------------------
repro_run_script <- function(script, root, log_dir, tag = "") {
  log <- file.path(log_dir, paste0("S19_run_", sub("\\.R$", "", basename(script)), tag, ".log"))
  old <- setwd(root); on.exit(setwd(old)); t0 <- Sys.time()
  st <- suppressWarnings(system2(repro_rscript(), shQuote(script), stdout = log, stderr = log))
  data.frame(script = basename(script), status = if (identical(as.integer(st), 0L)) "ok" else paste("FAILED", st), seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
             log = basename(log), stringsAsFactors = FALSE)
}
repro_run_tests <- function(root, log_dir) {
  code <- file.path(log_dir, "S19_run_tests.R"); log <- file.path(log_dir, "S19_run_tests.log")
  writeLines(c("res <- testthat::test_dir(file.path(here::here(), 'tests', 'testthat'), reporter = 'silent', stop_on_failure = FALSE)",
               "d <- as.data.frame(res)",
               "cat(sprintf('S19_TESTS passed=%d failed=%d skipped=%d errors=%d\\n', sum(d$passed), sum(d$failed), sum(d$skipped), sum(d$error)))"), code)
  old <- setwd(root); on.exit(setwd(old))
  suppressWarnings(system2(repro_rscript(), shQuote(code), stdout = log, stderr = log))
  l <- readLines(log, warn = FALSE); l <- l[grepl("^S19_TESTS", l)]
  if (!length(l)) return(c(passed = NA, failed = NA, skipped = NA, errors = NA))
  l <- sub("^S19_TESTS ", "", l[1]); v <- as.integer(regmatches(l, gregexpr("[0-9]+", l))[[1]]); names(v) <- c("passed", "failed", "skipped", "errors"); v
}
repro_md5_set <- function(dir, pattern, exclude = character(0)) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE); f <- f[!basename(f) %in% exclude]
  if (!length(f)) return(setNames(character(0), character(0)))
  setNames(unname(tools::md5sum(f)), basename(f))
}

# a changed frozen file is "changed (documented)" if docs/post_freeze_changes.csv lists it with the current MD5 (column md5_after); otherwise it stays "changed"
repro_mark_documented <- function(root, m, now) {
  pf <- file.path(root, "docs", "post_freeze_changes.csv")
  if (!file.exists(pf)) return(m$status)
  pc <- utils::read.csv(pf, stringsAsFactors = FALSE)
  if (!all(c("file", "md5_after") %in% names(pc))) return(m$status)
  ok <- m$status == "changed" & !is.na(now) & vapply(seq_along(now), function(i) isTRUE(any(pc$file == m$file[i] & pc$md5_after == now[i])), NA)
  m$status[ok] <- "changed (documented)"
  m$status
}
