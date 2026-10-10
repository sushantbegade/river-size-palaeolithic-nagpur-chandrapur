REGION_BBOX <- c(xmin = 76.5, ymin = 17.5, xmax = 82.5, ymax = 24.0)

EXT <- function(...) P("data", "external", ...)

region_wkt <- function(bbox = REGION_BBOX) {
  sf::st_as_text(sf::st_as_sfc(sf::st_bbox(bbox, crs = sf::st_crs(4326))))
}

is_valid_zip <- function(path) {
  if (!file.exists(path)) return(FALSE)
  ok <- tryCatch({ utils::unzip(path, list = TRUE); TRUE }, error = function(e) FALSE,
                 warning = function(w) FALSE)
  isTRUE(ok)
}

download_if_missing <- function(url, dest, md5 = NA_character_) {
  if (is_valid_zip(dest)) {
    message("already present: ", basename(dest))
    return(invisible(dest))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(7200, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  message("downloading ", basename(dest), " (large file, be patient) ...")
  res <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl"),
                  error = function(e) e, warning = function(w) w)
  if (inherits(res, "condition") || !is_valid_zip(dest)) {
    if (file.exists(dest)) unlink(dest)
    stop("Download failed: ", url,
         "\n  MANUAL OPTION: open the link in a browser, save the file exactly as",
         "\n  ", dest, "\n  then run this script again.", call. = FALSE)
  }
  if (!is.na(md5) && nzchar(md5)) {
    got <- unname(tools::md5sum(dest))
    if (!identical(tolower(got), tolower(md5))) {
      warning("md5 differs from the published value for ", basename(dest),
              " (got ", got, ", expected ", md5, ")", call. = FALSE)
    }
  }
  invisible(dest)
}

figshare_file <- function(article_id, name_pattern) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is needed: install.packages('jsonlite')", call. = FALSE)
  }
  meta <- jsonlite::fromJSON(paste0("https://api.figshare.com/v2/articles/", article_id))
  f <- meta$files
  hit <- f[grepl(name_pattern, f$name, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(hit)) {
    stop("No file matching '", name_pattern, "' in figshare article ", article_id,
         ". Files listed: ", paste(f$name, collapse = ", "), call. = FALSE)
  }
  md5 <- if ("computed_md5" %in% names(hit)) hit$computed_md5[1] else NA_character_
  list(name = hit$name[1], url = hit$download_url[1], size = hit$size[1], md5 = md5)
}

unzip_once <- function(zip, exdir) {
  if (!dir.exists(exdir) || !length(list.files(exdir, recursive = TRUE))) {
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(zip, exdir = exdir)
  }
  invisible(exdir)
}

find_shp <- function(dir, pattern) {
  all <- list.files(dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  hit <- all[grepl(pattern, basename(all), ignore.case = TRUE)]
  if (!length(hit)) {
    stop("No shapefile matching '", pattern, "' under ", dir,
         ". Found: ", paste(basename(all), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

read_region <- function(shp, label) {
  x <- sf::st_read(shp, wkt_filter = region_wkt(), quiet = TRUE)
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  message(label, ": ", nrow(x), " features in region")
  x <- sf::st_transform(x, CRS_ANALYSIS)
  if (any(!sf::st_is_valid(x))) x <- sf::st_make_valid(x)
  x[!sf::st_is_empty(x), ]
}

drop_coord_cols <- function(df) {
  df[, !grepl("^(lat|lon|long|latitude|longitude|x|y)$", names(df), ignore.case = TRUE), drop = FALSE]
}