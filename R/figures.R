# R/figures.R  (S18 final; REPLACES the earlier R/figures.R)
# Base graphics only. Rule: every plotted value is a saved-table value.
#  - Each figure has a builder fig*_data(tb) that returns a FigSource table (what is drawn).
#  - Each draw function draws ONLY from that table.
#  - fig_audit() re-reads the WRITTEN FigSource CSVs and the original saved tables and checks
#    every y / lo / hi against the cited table, row and column (0 mismatches required).
# Needs: sf. Optional: ragg. Uses project helpers: P(), CRS_ANALYSIS.

FIG_SRC_COLS <- c("figure", "panel", "series", "x", "x_label", "y", "lo", "hi", "interval_type",
                  "label", "note", "src_table", "src_row", "y_col", "lo_col", "hi_col")
FIG_COL <- c(S1 = "#0072B2", S2 = "#D55E00", S3 = "#009E73", S4 = "#7F7F7F")
FIG_PCH <- c(S1 = 19, S2 = 15, S3 = 17, S4 = 18)
# rd_csv(): reads one saved table from output/tables (searches sub-folders if not at top level); NULL if absent.
# Defined here only if your utils.R does not already define it (it was probably in your old R/figures.R).
if (!exists("rd_csv", mode = "function")) {
  rd_csv <- function(name) {
    f <- P("output", "tables", name)
    if (!file.exists(f)) { h <- list.files(P("output", "tables"), pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"), recursive = TRUE, full.names = TRUE)
    h <- h[!grepl("/paper/", h)]; if (!length(h)) return(NULL); f <- h[1] }
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  }
}
fig_or <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
fig_set_name <- function(s) unname(c(S1 = "S1 open-air", S2 = "S2 all sites", S3 = "S3 riverbed", S4 = "S4 other")[s])

# ------------------------------------------------------------------ formatting
# 2 decimals, but 3 decimals when |x| < 0.01 so that -0.002 never prints as "-0.00" or 0.003 as "0.00".
fig_fmt <- function(x, d = 2) {
  out <- sprintf(paste0("%.", d, "f"), x)
  small <- !is.na(x) & x != 0 & abs(x) < 0.01
  out[small] <- sprintf("%.3f", x[small]); out[is.na(x)] <- "NA"; out
}
fig_est_ci <- function(b, lo, hi, d = 2)
  ifelse(is.na(lo) | is.na(hi), fig_fmt(b, d), paste0(fig_fmt(b, d), " (", fig_fmt(lo, d), " to ", fig_fmt(hi, d), ")"))
fig_md_table <- function(df) {
  df[] <- lapply(df, function(z) { z <- as.character(z); z[is.na(z)] <- ""; gsub("\\|", "/", z) })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    vapply(seq_len(nrow(df)), function(i) paste0("| ", paste(unlist(df[i, ]), collapse = " | "), " |"), ""))
}

# ------------------------------------------------------------------ source tables + audit
fig_finalise <- function(df) {
  for (nm in setdiff(FIG_SRC_COLS, names(df))) df[[nm]] <- if (nm %in% c("x", "y", "lo", "hi", "src_row")) NA_real_ else ""
  df[, FIG_SRC_COLS, drop = FALSE]
}
fig_evalcol <- function(tab, col, i) {
  if (is.na(col) || col == "") return(rep(NA_real_, length(i)))
  if (grepl(" - ", col, fixed = TRUE)) { p <- strsplit(col, " - ", fixed = TRUE)[[1]]; return(tab[[p[1]]][i] - tab[[p[2]]][i]) }
  if (grepl(" + ", col, fixed = TRUE)) { p <- strsplit(col, " + ", fixed = TRUE)[[1]]; return(tab[[p[1]]][i] + tab[[p[2]]][i]) }
  tab[[col]][i]
}
fig_pick <- function(tab, tabname, idx, ycol, locol = NA_character_, hicol = NA_character_, ...) {
  df <- data.frame(y = fig_evalcol(tab, ycol, idx), lo = fig_evalcol(tab, locol, idx), hi = fig_evalcol(tab, hicol, idx),
                   src_table = tabname, src_row = idx, y_col = ycol, lo_col = locol, hi_col = hicol, stringsAsFactors = FALSE)
  ex <- list(...); for (nm in names(ex)) df[[nm]] <- ex[[nm]]
  fig_finalise(df)
}
# reader: function(file_name) -> data.frame (use the project's rd_csv). tol: CSV round trip.
fig_audit <- function(figsource_files, reader, tol = 1e-9) {
  cache <- list(); n <- 0L; bad <- character(0)
  same <- function(a, b) (is.na(a) && is.na(b)) || (!is.na(a) && !is.na(b) && abs(a - b) <= tol)
  for (f in figsource_files) {
    s <- utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    names(s)[names(s) == "y_val"] <- "y"; names(s)[names(s) == "x_pos"] <- "x"   # written names avoid coordinate-like headers
    s <- s[!is.na(s$src_row), , drop = FALSE]
    for (i in seq_len(nrow(s))) {
      tb <- s$src_table[i]; if (is.null(cache[[tb]])) cache[[tb]] <- reader(tb)
      if (is.null(cache[[tb]])) { bad <- c(bad, paste("cannot read", tb)); next }
      for (pr in list(c("y", "y_col"), c("lo", "lo_col"), c("hi", "hi_col"))) {
        n <- n + 1L
        if (!same(s[[pr[1]]][i], fig_evalcol(cache[[tb]], s[[pr[2]]][i], s$src_row[i])))
          bad <- c(bad, sprintf("%s row %d: %s row %d (%s)", basename(f), i, tb, s$src_row[i], pr[1]))
      }
    }
  }
  list(n_checked = n, n_bad = length(bad), bad = bad)
}

# ------------------------------------------------------------------ devices
fig_open_png <- function(file, width, height, res = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file, width = width, height = height, units = "in", res = res)
    return(if (.Platform$OS.type == "windows") "Times New Roman" else "serif")
  }
  if (.Platform$OS.type == "windows") grDevices::png(file, width = width, height = height, units = "in", res = res)
  else grDevices::png(file, width = width, height = height, units = "in", res = res, type = "cairo")
  "serif"
}
fig_save <- function(stem, width, height, draw_fn, outdir) {
  f_png <- file.path(outdir, paste0(stem, ".png")); fam <- fig_open_png(f_png, width, height)
  r1 <- try(draw_fn(fam), silent = TRUE); grDevices::dev.off()
  f_pdf <- file.path(outdir, paste0(stem, ".pdf"))
  if (isTRUE(capabilities("cairo"))) { grDevices::cairo_pdf(f_pdf, width = width, height = height, family = fam); pf <- fam }   # embeds fonts
  else { grDevices::pdf(f_pdf, width = width, height = height, family = "Times", useDingbats = FALSE); pf <- "Times" }
  r2 <- try(draw_fn(pf), silent = TRUE); grDevices::dev.off()
  for (r in list(r1, r2)) if (inherits(r, "try-error")) stop("Drawing failed for ", stem, ": ", as.character(r), call. = FALSE)
  c(f_png, f_pdf)
}
fig_ci <- function(pos, lo, hi, col = "black", lwd = 1.2, cap = 0.03, horiz = FALSE) {
  ok <- is.finite(lo) & is.finite(hi); if (!any(ok)) return(invisible())
  if (horiz) graphics::arrows(lo[ok], pos[ok], hi[ok], pos[ok], angle = 90, code = 3, length = cap, col = col, lwd = lwd)
  else graphics::arrows(pos[ok], lo[ok], pos[ok], hi[ok], angle = 90, code = 3, length = cap, col = col, lwd = lwd)
}
fig_panel_title <- function(lab, txt, cex = 0.95) graphics::mtext(bquote(bold(.(lab)) ~ .(txt)), side = 3, line = 0.5, adj = 0, cex = cex)
fig_stamp <- function(txt) graphics::mtext(txt, side = 3, outer = TRUE, line = 0.2, font = 2, cex = 0.85, col = "firebrick4")

# ================================================================== FIG 2: threshold profile
fig2_data <- function(tb) {
  p <- tb$S11_profile; s <- tb$S12_sample_sets_profile
  lab6 <- ifelse(p$k == 6, sprintf("intensity ratio q25 vs q75 = %s (%s to %s)", fig_fmt(p$ratio_q25_vs_q75), fig_fmt(p$ratio_lo), fig_fmt(p$ratio_hi)), "")
  nt <- function(k) ifelse(k == 6, "confirmatory k", ifelse(k >= 8, "exploratory (k8-k9)", ""))
  s1 <- function(panel) fig_pick(p, "S11_profile.csv", seq_len(nrow(p)), "beta_iqr", "lo", "hi", figure = "Fig2", panel = panel, series = "S1",
                                 x = p$k, x_label = paste0("k", p$k), interval_type = sprintf("95%% CI, full-refit cluster bootstrap, 40 km blocks, %d resamples", p$n_boot_ok),
                                 label = if (panel == "a") lab6 else "", note = nt(p$k))
  other <- function(set, ser) { idx <- which(s$subset == set); idx <- idx[order(s$k[idx])]
  fig_pick(s, "S12_sample_sets_profile.csv", idx, "beta_iqr", "lo", "hi", figure = "Fig2", panel = "b", series = ser, x = s$k[idx],
           x_label = paste0("k", s$k[idx]), interval_type = sprintf("95%% CI, one-step cluster bootstrap, 40 km blocks, %d resamples", s$n_boot_ok[idx]), note = nt(s$k[idx])) }
  rbind(s1("a"), s1("b"), other("S2_all", "S2"), other("S3_riverbed", "S3"))
}
fig2_draw <- function(src, fam) {
  op <- graphics::par(family = fam, cex.lab = 0.85, mfrow = c(1, 2), mar = c(4.2, 4.6, 2.4, 0.8), mgp = c(2.7, 0.7, 0), tcl = -0.3, las = 1, cex.axis = 0.85)
  on.exit(graphics::par(op))
  yl <- range(c(src$lo, src$hi, 0), na.rm = TRUE); yl <- yl + c(-0.06, 0.16) * diff(yl)
  for (pn in c("a", "b")) {
    d <- src[src$panel == pn, ]
    graphics::plot(NA, xlim = c(2.5, 9.5), ylim = yl, xaxt = "n", xlab = "Threshold k (WRIS order >= k)",
                   ylab = if (pn == "a") "Log intensity ratio per IQR (beta)" else "")
    graphics::rect(7.5, yl[1], 9.5, yl[2], col = "grey94", border = NA)
    graphics::text(8.5, yl[2] - 0.03 * diff(yl), "exploratory", cex = 0.7, font = 3, col = "grey35")
    graphics::abline(h = 0, lty = 2, col = "grey35"); graphics::axis(1, at = 3:9); graphics::box()
    ser <- unique(d$series); off <- if (length(ser) == 1) c(S1 = 0) else c(S1 = -0.22, S2 = 0, S3 = 0.22)
    for (sr in ser) { e <- d[d$series == sr, ]; xx <- e$x + off[[sr]]
    fig_ci(xx, e$lo, e$hi, FIG_COL[[sr]])
    graphics::points(xx, e$y, pch = ifelse(e$x >= 8, 1, FIG_PCH[[sr]]), col = FIG_COL[[sr]], cex = ifelse(e$x == 6 & sr == "S1", 1.5, 1.0), lwd = 1.2) }
    if (pn == "a") { k6 <- d[d$x == 6, ]; graphics::text(2.6, yl[2] - 0.06 * diff(yl), paste(strwrap(paste0("k = 6: ", k6$label), 30), collapse = "\n"), cex = 0.62, adj = c(0, 1)) }
    else graphics::legend("bottomright", legend = fig_set_name(ser), pch = FIG_PCH[ser], col = FIG_COL[ser], bty = "n", cex = 0.8)
    fig_panel_title(pn, if (pn == "a") "S1 open-air, full-refit bootstrap" else "S1, S2, S3 (S2, S3: one-step CI)")
  }
  invisible()
}

# ================================================================== FIG 3: distance classes (EXPLORATORY)
fig3_data <- function(tb) {
  d <- tb$S17_distance_shape; out <- list(); j <- 0; pl <- letters[1:9]
  for (set in c("S1_openair", "S2_all", "S3_riverbed")) for (k in c(4, 6, 8)) {
    j <- j + 1; idx <- which(d$subset == set & d$k == k); idx <- idx[order(d$bin[idx])]; ser <- substr(set, 1, 2)
    out[[j]] <- fig_pick(d, "S17_distance_shape.csv", idx, "rel_intensity", "lo", "hi", figure = "Fig3", panel = pl[j], series = ser, x = d$bin[idx],
                         x_label = sprintf("%g-%g", round(d$km_from[idx], 2), round(d$km_to[idx], 2)), interval_type = "95% interval as saved in S17 (reference class has none)",
                         label = sprintf("%s, k = %d", fig_set_name(ser), k), note = ifelse(is.na(d$lo[idx]), "reference class", ""))
  }
  do.call(rbind, out)
}
fig3_draw <- function(src, fam) {
  op <- graphics::par(family = fam, cex.lab = 0.85, mfrow = c(3, 3), mar = c(5.0, 3.8, 2.0, 0.6), oma = c(0, 0, 1.6, 0), mgp = c(2.4, 0.7, 0), tcl = -0.3, las = 1, cex.axis = 0.8)
  on.exit(graphics::par(op))
  yl <- range(c(src$lo, src$hi, src$y), na.rm = TRUE) * c(0.85, 1.15)
  tk <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50); tk <- tk[tk >= yl[1] & tk <= yl[2]]
  for (pn in unique(src$panel)) {
    d <- src[src$panel == pn, ]; sr <- d$series[1]
    graphics::plot(NA, xlim = c(0.5, 5.5), ylim = yl, log = "y", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::abline(h = 1, lty = 3, col = "grey40")
    graphics::axis(1, at = d$x, labels = d$x_label, las = 2, cex.axis = 0.62); graphics::axis(2, at = tk, labels = tk); graphics::box()
    graphics::mtext("Distance class (km)", side = 1, line = 3.7, cex = 0.55)
    if (pn %in% c("a", "d", "g")) graphics::mtext("Relative intensity (class 3 = 1)", side = 2, line = 2.4, cex = 0.62, las = 0)
    fig_ci(d$x, d$lo, d$hi, FIG_COL[[sr]], lwd = 1.1, cap = 0.025)
    ref <- is.na(d$lo)
    graphics::points(d$x[!ref], d$y[!ref], pch = FIG_PCH[[sr]], col = FIG_COL[[sr]], cex = 1)
    graphics::points(d$x[ref], d$y[ref], pch = 1, col = FIG_COL[[sr]], cex = 1.1, lwd = 1.2)
    fig_panel_title(pn, d$label[1], cex = 0.85)
  }
  fig_stamp("EXPLORATORY: five equal-area distance classes per k; open symbol = reference class")
  invisible()
}

# ================================================================== FIG 4: spatial cross-validation
fig4_data <- function(tb) {
  d <- tb$S14_cv_summary; out <- list(); j <- 0; pl <- letters[1:3]; mods <- c(paste0("M2_river_k", 3:9), "M3_lineament")
  for (set in c("S1_openair", "S2_all", "S3_riverbed")) { j <- j + 1
  for (blk in c(40, 20)) {
    idx <- which(d$set == set & d$block_km == blk & d$model %in% mods); idx <- idx[order(match(d$model[idx], mods))]; ser <- substr(set, 1, 2)
    out[[length(out) + 1]] <- fig_pick(d, "S14_cv_summary.csv", idx, "delta_vs_nuisance", "delta_vs_nuisance - se_delta_sum", "delta_vs_nuisance + se_delta_sum",
                                       figure = "Fig4", panel = pl[j], series = paste0(ser, " ", blk, " km"), x = match(d$model[idx], mods),
                                       x_label = ifelse(d$model[idx] == "M3_lineament", "Lin", sub("M2_river_", "", d$model[idx])),
                                       interval_type = "+-1 SE of summed fold differences (not a CI)", label = fig_set_name(ser),
                                       note = sprintf("folds improved %d, worse %d, used %d", d$folds_improved[idx], d$folds_worse[idx], d$n_folds_used[idx]))
  } }
  do.call(rbind, out)
}
fig4_draw <- function(src, fam) {
  op <- graphics::par(family = fam, cex.lab = 0.85, mfrow = c(1, 3), mar = c(4.2, 4.4, 2.4, 0.6), mgp = c(2.6, 0.7, 0), tcl = -0.3, las = 1, cex.axis = 0.8)
  on.exit(graphics::par(op))
  yl <- range(c(src$lo, src$hi, 0), na.rm = TRUE); yl <- yl + c(-0.04, 0.04) * diff(yl)
  for (pn in unique(src$panel)) {
    d <- src[src$panel == pn, ]; sr <- substr(d$series[1], 1, 2)
    graphics::plot(NA, xlim = c(0.5, 8.5), ylim = yl, xaxt = "n", xlab = "River threshold k, or lineament (Lin)",
                   ylab = if (pn == "a") "Held-out log-likelihood gain vs nuisance-only" else "")
    graphics::abline(h = 0, lty = 2, col = "grey35")
    u <- !duplicated(d$x); graphics::axis(1, at = d$x[u], labels = d$x_label[u], cex.axis = 0.8); graphics::box()
    for (s in unique(d$series)) { e <- d[d$series == s, ]; is40 <- grepl(" 40 km", s); xx <- e$x + ifelse(is40, -0.13, 0.13)
    fig_ci(xx, e$lo, e$hi, FIG_COL[[sr]], lwd = 1.1, cap = 0.025)
    graphics::points(xx, e$y, pch = if (is40) FIG_PCH[[sr]] else 1, col = FIG_COL[[sr]], cex = 1.0, lwd = 1.2) }
    if (pn == "a") graphics::legend("topright", legend = c("40 km blocks (+-1 SE)", "20 km blocks (+-1 SE)"), pch = c(19, 1), col = "grey20", bty = "n", cex = 0.7)
    fig_panel_title(pn, d$label[1])
  }
  invisible()
}

# ================================================================== FIG 5: robustness (k = 6)
fig5_data <- function(tb) {
  keys <- character(0); out <- list()
  push <- function(group, label, set, tab, tabname, idx, interval) {
    if (!length(idx)) return(invisible())
    key <- paste(group, label, sep = "||"); if (!key %in% keys) keys <<- c(keys, key)
    out[[length(out) + 1]] <<- fig_pick(tab, tabname, idx, "beta_iqr", "lo", "hi", figure = "Fig5", panel = if (set == "S1") "a" else "b",
                                        series = set, x = match(key, keys), x_label = label, interval_type = interval, note = group)
  }
  g1 <- "Bootstrap block size (full refit, S1 only)"; g2 <- "Model and domain variants (one-step bootstrap)"
  g3 <- "Coordinate tier (S1 only)"; g4 <- "Cluster definition (one-step bootstrap)"; g5 <- "Reference-conditional likelihood"
  pk <- tb$S11_primary_k6
  for (b in c(40, 20, 10)) { i <- which(pk$block_km == b & pk$method == "full-refit bootstrap")
  push(g1, if (b == 40) "40 km blocks (primary)" else paste0(b, " km blocks"), "S1", pk, "S11_primary_k6.csv", i, sprintf("95%% CI, full-refit cluster bootstrap, %d resamples", pk$n_boot_ok[i])) }
  v <- tb$S16_variants_profile
  tidy <- function(z) { z <- sub("^ADDITIONAL: ", "extra (not pre-registered): ", z); z <- sub("^baseline$", "baseline (one-step comparator)", z)
  z <- sub("d_road_any instead of d_road_major", "any road instead of major road", z); sub("HydroRIVERS extension off \\(WRIS network only\\)", "HydroRIVERS extension off", z) }
  for (var in unique(v$variant)) for (set in c("S1", "S2")) {
    i <- which(v$variant == var & v$subset == paste0(set, ifelse(set == "S1", "_openair", "_all")) & v$k == 6)
    push(g2, tidy(var), set, v, "S16_variants_profile.csv", i, sprintf("95%% CI, one-step cluster bootstrap, %d resamples", v$n_boot_ok[i])) }
  t <- tb$S12_tiers_profile
  for (sb in c("S1_tier_A", "S1_tier_AB", "S1_tier_B", "S1_tier_C")) { i <- which(t$subset == sb & t$k == 6)
  push(g3, sprintf("tier %s (n = %d)", sub("S1_tier_", "", sb), t$n_sites[i]), "S1", t, "S12_tiers_profile.csv", i, sprintf("95%% CI, one-step cluster bootstrap, %d resamples", t$n_boot_ok[i])) }
  cdf <- tb$S15_cluster_definitions
  for (cl in setdiff(unique(cdf$clusters), "block 40 km (primary)")) for (set in c("S1", "S2")) {
    i <- which(cdf$clusters == cl & cdf$sample_set == paste0(set, ifelse(set == "S1", "_openair", "_all")) & cdf$k == 6)
    push(g4, paste("clusters:", cl), set, cdf, "S15_cluster_definitions.csv", i, "95% interval as saved in S15 (cluster bootstrap)") }
  cp <- tb$S15_conditional_profile
  for (set in c("S1", "S2")) { i <- which(cp$sample_set == paste0(set, ifelse(set == "S1", "_openair", "_all")) & cp$k == 6)
  push(g5, "reference-conditional likelihood", set, cp, "S15_conditional_profile.csv", i, sprintf("95%% interval as saved in S15 conditional, %d resamples", cp$n_boot_ok[i])) }
  do.call(rbind, out)
}
fig5_draw <- function(src, fam) {
  op <- graphics::par(family = fam, cex.lab = 0.85, tcl = -0.3, las = 1, mgp = c(2.5, 0.7, 0)); on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(2.25, 1))
  rows <- unique(src[, c("x", "x_label", "note")]); rows <- rows[order(rows$x), ]
  rows$ypos <- rows$x + c(0, cumsum(rows$note[-1] != rows$note[-nrow(rows)]) * 0.9); ymax <- max(rows$ypos)
  xl <- range(c(src$lo, src$hi, 0), na.rm = TRUE); xl <- xl + c(-0.04, 0.04) * diff(xl)
  for (pn in c("a", "b")) {
    d <- src[src$panel == pn, ]; sr <- d$series[1]
    graphics::par(mar = c(4.2, if (pn == "a") 17.5 else 0.6, 2.4, 0.8), cex.lab = 0.85)
    graphics::plot(NA, xlim = xl, ylim = c(ymax + 0.6, 0.4), yaxt = "n", xaxt = "n", xlab = "beta per IQR (k = 6)", ylab = "")
    graphics::axis(1, at = seq(-0.5, 2.5, 0.5), labels = sprintf("%.1f", seq(-0.5, 2.5, 0.5)), cex.axis = 0.7)
    graphics::abline(v = 0, lty = 2, col = "grey35")
    if (pn == "a") { gs <- !duplicated(rows$note)
    ux <- graphics::par("usr")[1] - 0.02 * diff(graphics::par("usr")[1:2])
    graphics::text(ux, rows$ypos, rows$x_label, adj = 1, xpd = NA, cex = 0.7, font = ifelse(rows$x == 1, 2, 1))
    graphics::text(ux, rows$ypos[gs] - 0.8, rows$note[gs], adj = 1, xpd = NA, cex = 0.62, font = 3, col = "grey30") }
    yy <- rows$ypos[match(d$x, rows$x)]
    fig_ci(yy, d$lo, d$hi, FIG_COL[[sr]], lwd = 1.1, cap = 0.025, horiz = TRUE)
    graphics::points(d$y, yy, pch = FIG_PCH[[sr]], col = FIG_COL[[sr]], cex = ifelse(d$x == 1 & sr == "S1", 1.4, 0.95))
    graphics::box(); fig_panel_title(pn, fig_set_name(sr))
  }
  graphics::layout(1); invisible()
}

# ================================================================== FIG 6: site type (EXPLORATORY)
fig6_data <- function(tb) {
  s <- tb$S13_site_type_profile; dc <- tb$S13_decomposition_k6; out <- list()
  for (m in c("distance only", "distance + nuisance")) { i <- which(s$model == m & s$fit_ok %in% TRUE); i <- i[order(s$k[i])]
  out[[length(out) + 1]] <- fig_pick(s, "S13_site_type_profile.csv", i, "beta_iqr", "lo", "hi", figure = "Fig6", panel = "a", series = m, x = s$k[i],
                                     x_label = paste0("k", s$k[i]), interval_type = sprintf("95%% CI, cluster bootstrap, %d resamples", s$n_boot_ok[i]),
                                     label = sprintf("P(riverbed) at q25 %s, at q75 %s", fig_fmt(s$p_riverbed_at_q25[i], 3), fig_fmt(s$p_riverbed_at_q75[i], 3)),
                                     note = "EXPLORATORY; riverbed vs open-air among S1+S3 sites") }
  out[[length(out) + 1]] <- fig_pick(dc, "S13_decomposition_k6.csv", seq_len(nrow(dc)), "value", figure = "Fig6", panel = "b", series = "decomposition",
                                     x = seq_len(nrow(dc)), x_label = dc$quantity, interval_type = "none (point value, derived)", note = "EXPLORATORY")
  do.call(rbind, out)
}
fig6_draw <- function(src, fam) {
  op <- graphics::par(family = fam, cex.lab = 0.85, cex.axis = 0.8, tcl = -0.3, las = 1, mgp = c(2.5, 0.7, 0), oma = c(0, 0, 1.6, 0)); on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(1, 1.45))
  a <- src[src$panel == "a", ]; graphics::par(mar = c(4.2, 4.6, 2.4, 0.8))
  yl <- range(c(a$lo, a$hi, 0), na.rm = TRUE); yl <- yl + c(-0.04, 0.04) * diff(yl)
  graphics::plot(NA, xlim = c(2.5, 9.5), ylim = yl, xaxt = "n", xlab = "Threshold k (WRIS order >= k)", ylab = "Log-odds of riverbed per IQR")
  graphics::abline(h = 0, lty = 2, col = "grey35"); graphics::axis(1, at = 3:9); graphics::box()
  for (m in unique(a$series)) { e <- a[a$series == m, ]; o <- ifelse(m == "distance only", -0.12, 0.12)
  fig_ci(e$x + o, e$lo, e$hi, "grey20", lwd = 1.1, cap = 0.025); graphics::points(e$x + o, e$y, pch = ifelse(m == "distance only", 1, 19), col = "grey20") }
  graphics::legend("bottomright", legend = c("distance only", "distance + nuisance"), pch = c(1, 19), col = "grey20", bty = "n", cex = 0.75)
  fig_panel_title("a", "Riverbed vs open-air (negative = nearer rivers)", cex = 0.8)
  b <- src[src$panel == "b", ]; graphics::par(mar = c(4.2, 12.5, 2.4, 0.8))
  lab <- vapply(strwrap(b$x_label, width = 34, simplify = FALSE), paste, "", collapse = "\n")
  xl <- range(c(b$y, 0)); xl <- xl + c(-0.2, 0.12) * diff(xl)
  graphics::plot(NA, xlim = xl, ylim = c(nrow(b) + 0.5, 0.5), yaxt = "n", xlab = "Value at k = 6", ylab = "")
  graphics::abline(v = 0, lty = 2, col = "grey35")
  graphics::text(graphics::par("usr")[1] - 0.03 * diff(graphics::par("usr")[1:2]), b$x, lab, adj = 1, xpd = NA, cex = 0.6, lheight = 0.85)
  graphics::segments(0, b$x, b$y, b$x, col = "grey60"); graphics::points(b$y, b$x, pch = 19, col = "grey20")
  graphics::text(b$y, b$x - 0.3, fig_fmt(b$y, 3), cex = 0.62); graphics::box()
  fig_panel_title("b", "Decomposition, k = 6 (no intervals)", cex = 0.8)
  fig_stamp("EXPLORATORY: site-type composition among open-air and riverbed sites")
  graphics::layout(1); invisible()
}

# ================================================================== FIG 1: study-area map
# Inputs (your layout): data/processed/hydro_prep.gpkg (layers boundary, rivers[ordsh], removed_core),
# data/public/sites_public_1km.csv (x_c, y_c, S1_openair, S3_riverbed), data/external/india_boundary.gpkg (optional).
fig1_load <- function(crs) {
  fh <- P("data", "processed", "hydro_prep.gpkg"); fp <- P("data", "public", "sites_public_1km.csv"); fi <- P("data", "external", "india_boundary.gpkg")
  if (!file.exists(fh)) stop("Missing ", fh, call. = FALSE); if (!file.exists(fp)) stop("Missing ", fp, call. = FALSE)
  pub <- utils::read.csv(fp, stringsAsFactors = FALSE)
  need <- c("x_c", "y_c", "S1_openair", "S3_riverbed"); if (!all(need %in% names(pub))) stop("Public file lacks: ", paste(setdiff(need, names(pub)), collapse = ", "), call. = FALSE)
  if (!(all(abs((pub$x_c %% 1000) - 500) < 1e-6) && all(abs((pub$y_c %% 1000) - 500) < 1e-6)))
    stop("sites_public_1km.csv is not on the 1 km cell-centre grid: refusing to draw it.", call. = FALSE)
  s1 <- as.logical(pub$S1_openair); s3 <- as.logical(pub$S3_riverbed)
  cls <- ifelse(s1, "S1 open-air", ifelse(s3, "S3 riverbed", "other"))
  rd <- function(l) sf::st_transform(sf::st_read(fh, layer = l, quiet = TRUE), crs)
  rv <- rd("rivers"); if (!"ordsh" %in% names(rv)) stop("rivers layer has no ordsh column", call. = FALSE)
  india <- if (file.exists(fi)) sf::st_transform(sf::st_read(fi, quiet = TRUE), 4326) else NULL
  list(pts = data.frame(x = pub$x_c, y = pub$y_c, cls = cls, stringsAsFactors = FALSE), bnd = rd("boundary"), rivers = rv, core = rd("removed_core"),
       india = india, md5 = unname(tools::md5sum(fp)))
}
fig_grat_step <- function(w) { cand <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 5); cand[which(w / cand <= 7)[1]] }
fig_deg_lab <- function(v, step, hemi) {
  dg <- if (step >= 1) 0 else if (abs(step * 10 - round(step * 10)) < 1e-9) 1 else 2
  paste0(formatC(abs(v), format = "f", digits = dg), "\u00b0", hemi)
}
fig_graticule <- function(xlim, ylim, crs) {
  bb <- sf::st_as_sfc(sf::st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]), crs = crs))
  ll <- sf::st_bbox(sf::st_transform(sf::st_segmentize(bb, 1000), 4326))
  sx <- fig_grat_step(ll[["xmax"]] - ll[["xmin"]]); sy <- fig_grat_step(ll[["ymax"]] - ll[["ymin"]])
  lons <- seq(ceiling(ll[["xmin"]] / sx) * sx, floor(ll[["xmax"]] / sx) * sx, by = sx); lats <- seq(ceiling(ll[["ymin"]] / sy) * sy, floor(ll[["ymax"]] / sy) * sy, by = sy)
  tr <- function(lon, lat) sf::st_coordinates(sf::st_transform(sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326), crs))
  lat_seq <- seq(ll[["ymin"]] - 0.3, ll[["ymax"]] + 0.3, length.out = 200); lon_seq <- seq(ll[["xmin"]] - 0.3, ll[["xmax"]] + 0.3, length.out = 200)
  mer <- lapply(lons, function(L) tr(rep(L, 200), lat_seq)); prl <- lapply(lats, function(Y) tr(lon_seq, rep(Y, 200)))
  list(lines = c(mer, prl), xt = vapply(mer, function(m) stats::approx(m[, 2], m[, 1], xout = ylim[1])$y, 0), xlab = fig_deg_lab(lons, sx, "E"),
       yt = vapply(prl, function(m) stats::approx(m[, 1], m[, 2], xout = xlim[1])$y, 0), ylab = fig_deg_lab(lats, sy, "N"))
}
fig_scalebar <- function(usr, pos = "bottomright") {
  w <- (usr[2] - usr[1]) / 1000; ch <- c(8, 20, 40, 80, 100); len <- ch[which.min(abs(ch - 0.22 * w))]
  x0 <- if (grepl("right", pos)) usr[2] - 0.04 * (usr[2] - usr[1]) - len * 1000 else usr[1] + 0.04 * (usr[2] - usr[1]); y0 <- usr[3] + 0.045 * (usr[4] - usr[3]); h <- 0.012 * (usr[4] - usr[3]); sg <- len * 1000 / 4
  for (i in 0:3) graphics::rect(x0 + i * sg, y0, x0 + (i + 1) * sg, y0 + h, col = if (i %% 2 == 0) "black" else "white", border = "black", lwd = 0.5)
  graphics::text(x0 + c(0, 1, 2, 3, 4) * sg, y0 + h + 0.014 * (usr[4] - usr[3]), c("0", len / 4, len / 2, 3 * len / 4, paste(len, "km")), cex = 0.62)
}
fig_northarrow <- function(usr) {
  dy <- usr[4] - usr[3]; dx <- usr[2] - usr[1]; H <- 0.085 * dy; xc <- usr[2] - 0.06 * dx; yb <- usr[4] - 0.07 * dy - H; W <- 0.3 * H
  graphics::polygon(c(xc, xc - W, xc), c(yb + H, yb, yb + 0.3 * H), col = "black", border = "black", lwd = 0.5)
  graphics::polygon(c(xc, xc + W, xc), c(yb + H, yb, yb + 0.3 * H), col = "white", border = "black", lwd = 0.5)
  graphics::text(xc, yb + H + 0.025 * dy, "N", font = 2, cex = 0.9)
}
fig_inset_india <- function(india, ext, crs, pos = "bottomleft") {
  usr <- graphics::par("usr"); w <- 0.24 * (usr[2] - usr[1]); m <- 0.03 * (usr[2] - usr[1])
  x1 <- if (grepl("left", pos)) usr[1] + m else usr[2] - m - w; y2 <- if (grepl("top", pos)) usr[4] - m else usr[3] + m + w
  graphics::rect(x1, y2 - w, x1 + w, y2, col = "white", border = "black", lwd = 0.6)
  fg <- c(graphics::grconvertX(c(x1, x1 + w), "user", "ndc"), graphics::grconvertY(c(y2 - w, y2), "user", "ndc"))
  bb <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = ext$xlim[1], ymin = ext$ylim[1], xmax = ext$xlim[2], ymax = ext$ylim[2]), crs = crs)), 4326))
  op <- graphics::par(fig = fg, new = TRUE, mar = c(0.1, 0.1, 0.1, 0.1)); on.exit(graphics::par(op))
  graphics::plot(sf::st_geometry(india), col = "grey88", border = "grey30", lwd = 0.4, reset = FALSE)
  graphics::rect(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]], border = "red3", lwd = 1.0, col = NA)
}
fig1_data <- function(m) {
  cnt <- table(m$pts$cls)
  rows <- lapply(names(cnt), function(nm) data.frame(figure = "Fig1", panel = "a", series = nm, y = as.numeric(cnt[[nm]]), interval_type = "count of plotted public 1 km points",
                                                     note = paste("sites_public_1km.csv md5", m$md5), stringsAsFactors = FALSE))
  rows[[length(rows) + 1]] <- data.frame(figure = "Fig1", panel = "a", series = "total points drawn", y = nrow(m$pts), note = "1 km cell centres only (grid verified)", stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <- data.frame(figure = "Fig1", panel = "a", series = "river features", y = nrow(m$rivers), note = "hydro_prep.gpkg layer rivers; classes by ordsh", stringsAsFactors = FALSE)
  do.call(rbind, lapply(rows, fig_finalise))
}
fig1_draw <- function(m, fam, ext, crs, inset_pos = "bottomleft", scalebar_pos = "bottomright") {
  op <- graphics::par(family = fam, mar = c(7, 4.8, 1, 1), tcl = -0.3, las = 1, mgp = c(3, 0.7, 0)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = ext$xlim, ylim = ext$ylim, asp = 1, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")   # asp = 1: true shape and scale
  usr <- graphics::par("usr"); gr <- fig_graticule(ext$xlim, ext$ylim, crs)
  for (l in gr$lines) graphics::lines(l[, 1], l[, 2], col = "grey85", lwd = 0.5, lty = 3)
  graphics::plot(sf::st_geometry(m$core), add = TRUE, col = "#9ECAE1", border = NA)
  g <- sf::st_geometry(m$rivers); ord <- suppressWarnings(as.numeric(m$rivers$ordsh))
  drw <- function(ix, col, lwd) if (length(ix)) graphics::plot(g[ix], add = TRUE, col = col, lwd = lwd)
  drw(which(ord < 5), "#9ECAE1", 0.35); drw(which(ord == 5), "#4292C6", 1.0); drw(which(ord >= 6), "#08519C", 1.9)
  graphics::plot(sf::st_geometry(m$bnd), add = TRUE, border = "black", lwd = 0.8)
  sty <- list("S1 open-air" = c(21, "#0072B2"), "S3 riverbed" = c(24, "#009E73"), "other" = c(22, "#E69F00"))
  lt <- c("Other rivers (WRIS order < 5)", "WRIS order 5", "WRIS order >= 6 (confirmatory k = 6 network)", "Removed areas (ponds, reservoirs, quarries)")
  lp <- c(NA, NA, NA, 15); ll <- c(1, 1, 1, NA); lw <- c(0.6, 1.0, 1.9, NA); lc <- c("#9ECAE1", "#4292C6", "#08519C", "#9ECAE1"); lb <- rep(NA, 4)
  for (nm in names(sty)) { s <- m$pts$cls == nm; if (!any(s)) next
  graphics::points(m$pts$x[s], m$pts$y[s], pch = as.numeric(sty[[nm]][1]), bg = sty[[nm]][2], col = "black", lwd = 0.4, cex = 0.85)
  lt <- c(lt, sprintf("%s (n = %d, 1 km grid)", nm, sum(s))); lp <- c(lp, as.numeric(sty[[nm]][1])); ll <- c(ll, NA); lw <- c(lw, NA); lc <- c(lc, "black"); lb <- c(lb, sty[[nm]][2]) }
  graphics::box(lwd = 0.8)
  graphics::axis(1, at = gr$xt, labels = gr$xlab, cex.axis = 0.8); graphics::axis(2, at = gr$yt, labels = gr$ylab, cex.axis = 0.8)
  fig_scalebar(usr, scalebar_pos); fig_northarrow(usr)
  graphics::par(xpd = NA)
  graphics::legend(x = mean(usr[1:2]), y = usr[3] - 0.075 * (usr[4] - usr[3]), xjust = 0.5, yjust = 1, legend = lt, pch = lp, lty = ll, lwd = lw, col = lc, pt.bg = lb, ncol = 2, bty = "n", cex = 0.7, pt.cex = 0.9)
  graphics::par(xpd = FALSE)
  if (!is.null(m$india)) fig_inset_india(m$india, ext, crs, inset_pos)   # last: it switches the plot region
  invisible()
}