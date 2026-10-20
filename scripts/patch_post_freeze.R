# scripts/patch_post_freeze.R
# One-off, documented post-freeze edit (deviation log A32). It changes WORDING ONLY in code, one document title and two output labels:
# "pre-registered / pre-registration / OSF registration" -> "frozen design" wording, because the design was frozen and checksummed internally on
# 4 Oct 2026 but never registered externally. It also teaches the freeze check about documented edits (R/repro.R, scripts/19_reproducibility.R).
# Default is a DRY RUN. To apply:  Sys.setenv(PATCH_APPLY = "1"); source(here::here("scripts", "patch_post_freeze.R"))
# Originals are copied to output/logs/post_freeze_backup/ first; every changed file is logged in docs/post_freeze_changes.csv (old and new MD5).
suppressPackageStartupMessages(library(here))
APPLY <- Sys.getenv("PATCH_APPLY", "0") == "1"
root  <- here::here()
message("== post-freeze wording patch (", if (APPLY) "APPLY" else "DRY RUN", ") ==")

# fixed-string replacements, applied in this order to every target file
REPL <- list(
  c("(pre-registered design)", "(frozen design; not externally registered)"),
  c("(pre-registered secondary analyses on the real sites)", "(secondary analyses of the frozen design on the real sites)"),
  c("(pre-registered robustness matrix", "(robustness matrix of the frozen design"),
  c("(EXPLORATORY items from the pre-registration", "(EXPLORATORY items declared in the frozen design"),
  c("(pre-registration section 5)", "(frozen-design document, section 5)"),
  c("(pre-registered as exploratory)", "(declared exploratory in the frozen design)"),
  c("all items are pre-registered secondary/robustness analyses", "all items are secondary/robustness analyses of the frozen design"),
  c("NOT in the pre-registration (labelled)", "NOT in the frozen design (labelled)"),
  c("extra (not pre-registered)", "extra (not in the frozen design)"),
  c("were not pre-registered", "are not in the frozen design"),
  c("# --- pre-registered:", "# --- frozen design:"),
  c("# --- pre-registered ", "# --- frozen design "),
  c("(paste this into the OSF registration)", "(record this value in Online Resource 1; no external registration exists)"),
  c("Pre-registered rule", "Frozen-design rule"),
  c("pre-registered rule", "frozen-design rule"),
  c("(pre-registered)", "(frozen design)"),
  c("(pre-registration)", "(frozen design)"),
  c("(pre-registered", "(frozen design"),
  c("# Pre-registration (draft for OSF)", "# Design-freeze document (internal draft; never filed with any external registry; file name kept for checksum continuity)"))
ALLOWED <- c("docs/preregistration.md", "preregistration\\\\.md", "preregistration\\.md", "preregistration.md")   # file name stays: it is in the freeze manifest

find_one <- function(base, dirs = c("R", "scripts")) {
  hit <- unlist(lapply(dirs, function(d) { p <- list.files(file.path(root, d), pattern = paste0("^", gsub("\\.", "\\\\.", base), "$"), recursive = TRUE)
  if (length(p)) file.path(d, p) else character(0) }))
  if (length(hit) != 1L) stop("Cannot locate ", base, " uniquely under R/ or scripts/ (found ", length(hit), ").", call. = FALSE)
  hit
}
code <- vapply(c("11_primary_model.R", "12_secondary_diagnostics.R", "14_spatial_cv.R", "16_robustness.R", "17_exploratory.R", "18_figures_tables.R",
                 "clean_sites.R", "diagnostics.R", "exploratory.R", "robustness.R", "spatial_cv.R", "figures.R", "freeze_manifest.R"), find_one, "")
targets <- c(unname(code), "docs/preregistration.md", "output/tables/S11_decision.csv", "output/logs/S17_parameters.csv")

read_txt  <- function(p) rawToChar(readBin(p, "raw", file.size(p)))
write_txt <- function(p, s) writeBin(charToRaw(s), p)
md5 <- function(p) unname(tools::md5sum(p))
rep_all <- function(s) { n <- 0L
for (r in REPL) { n <- n + lengths(regmatches(s, gregexpr(r[1], s, fixed = TRUE))); s <- gsub(r[1], r[2], s, fixed = TRUE, useBytes = TRUE) }
list(text = s, n = n) }
residual <- function(s) { t <- s; for (a in ALLOWED) t <- gsub(a, "", t, fixed = TRUE, useBytes = TRUE)
ln <- strsplit(t, "\r?\n")[[1]]; which(grepl("pre-regist|preregist|OSF", ln, ignore.case = TRUE, useBytes = TRUE)) }

plan <- list(); bad <- FALSE
for (f in targets) {
  p <- file.path(root, f)
  if (!file.exists(p)) { message("  skip (not found): ", f); next }
  s <- read_txt(p); r <- rep_all(s); res <- residual(r$text)
  message(sprintf("  %-46s %2d replacement(s)%s", f, r$n, if (length(res)) paste0("   RESIDUAL at line(s) ", paste(res, collapse = ",")) else ""))
  if (length(res)) bad <- TRUE
  if (r$n > 0L) plan[[f]] <- list(before = s, after = r$text, n = r$n)
}
if (bad) stop("Residual 'pre-registered' / 'OSF' wording remains after the replacements; nothing was written.", call. = FALSE)

# documented-change support for the freeze check (R/repro.R and scripts/19_reproducibility.R)
repro_p <- file.path("R", "repro.R"); s19_p <- file.path("scripts", "19_reproducibility.R")
extra <- list(
  list(file = repro_p, old = 'm$status <- ifelse(!ex, "missing", ifelse(now == m$md5, "unchanged", "changed"))',
       new = 'm$status <- ifelse(!ex, "missing", ifelse(now == m$md5, "unchanged", "changed"))\n  m$status <- repro_mark_documented(root, m, now)'),
  list(file = s19_p, old = 'chg <- fz$table$file[fz$table$status == "changed"]; mis <-',
       new = 'chg <- fz$table$file[fz$table$status == "changed"]; chg_doc <- fz$table$file[fz$table$status == "changed (documented)"]; mis <-'),
  list(file = s19_p, old = "%d unchanged, %d changed, %d missing", new = "%d unchanged, %d changed (documented), %d changed (undocumented), %d missing"),
  list(file = s19_p, old = 'sum(fz$table$status == "unchanged"), length(chg), length(mis)', new = 'sum(fz$table$status == "unchanged"), length(chg_doc), length(chg), length(mis)'))
APPEND_REPRO <- '

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
'
for (e in extra) {
  p <- file.path(root, e$file); if (!file.exists(p)) stop("Missing ", e$file, call. = FALSE)
  s <- if (!is.null(plan[[e$file]])) plan[[e$file]]$after else read_txt(p)
  k <- lengths(regmatches(s, gregexpr(e$old, s, fixed = TRUE)))
  if (k != 1L) stop("Expected exactly one match in ", e$file, " for: ", e$old, " (found ", k, ")", call. = FALSE)
  plan[[e$file]] <- list(before = if (is.null(plan[[e$file]])) s else plan[[e$file]]$before, after = sub(e$old, e$new, s, fixed = TRUE, useBytes = TRUE),
                         n = (if (is.null(plan[[e$file]])) 0L else plan[[e$file]]$n) + 1L)
  message("  ", e$file, ": documented-change support (1 edit)")
}
sr <- read_txt(file.path(root, repro_p)); if (grepl("repro_mark_documented <- function", sr, fixed = TRUE)) stop("R/repro.R already patched.", call. = FALSE)
plan[[repro_p]]$after <- paste0(sub("\\s+$", "", plan[[repro_p]]$after), APPEND_REPRO)

message("files to change: ", length(plan), " | total edits: ", sum(vapply(plan, function(z) z$n, 0L)))
if (!APPLY) {
  message("DRY RUN only. Re-run with PATCH_APPLY=1 to write the files.")
  cat("\nDONE patch (dry run). Paste this console output back.\n")
} else {
  bk <- file.path(root, "output", "logs", "post_freeze_backup"); dir.create(bk, recursive = TRUE, showWarnings = FALSE)
  log <- NULL
  for (f in names(plan)) {
    p <- file.path(root, f); dir.create(file.path(bk, dirname(f)), recursive = TRUE, showWarnings = FALSE)
    file.copy(p, file.path(bk, f), overwrite = FALSE)
    m0 <- md5(p); write_txt(p, plan[[f]]$after); m1 <- md5(p)
    log <- rbind(log, data.frame(date = as.character(Sys.Date()), file = gsub("\\\\", "/", f), md5_before = m0, md5_after = m1, edits = plan[[f]]$n,
                                 change = if (f %in% c(repro_p, s19_p)) "documented-change support for the freeze check" else "wording only: pre-registered -> frozen design (no logic change)",
                                 deviation_id = "A32", stringsAsFactors = FALSE))
  }
  pf <- file.path(root, "docs", "post_freeze_changes.csv")
  old <- if (file.exists(pf)) utils::read.csv(pf, stringsAsFactors = FALSE) else NULL
  utils::write.csv(rbind(old, log), pf, row.names = FALSE)
  message("wrote docs/post_freeze_changes.csv (", nrow(log), " new row(s)); backups in output/logs/post_freeze_backup/")
  cat("\nDONE patch (applied). Paste this console output back.\n")
}