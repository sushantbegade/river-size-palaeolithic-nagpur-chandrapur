# R/utils.R -- small helpers

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
  message("wrote ", sub(paste0(here::here(), "/"), "", path, fixed = TRUE),
          "  (", nrow(x), " rows)")
  invisible(path)
}

file_checksums <- function(paths) {
  data.frame(file = basename(paths),
             md5 = unname(tools::md5sum(paths)),
             bytes = file.size(paths),
             stringsAsFactors = FALSE)
}
