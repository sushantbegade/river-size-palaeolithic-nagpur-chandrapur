# R/config.R  -- central constants and path helpers. Reads NO data.

CRS_ANALYSIS   <- 32644L   # WGS84 / UTM zone 44N (state this in Methods)
CRS_GEOGRAPHIC <- 4326L    # WGS84 lon/lat
SEED           <- 20261002L

# Coordinate sanity box for the two districts (generous). Order: lon, lat.
STUDY_BBOX <- c(xmin = 78.0, xmax = 80.5, ymin = 19.0, ymax = 22.0)

# Project-relative paths (always via here::here, never setwd())
P <- function(...) here::here(...)

paleo_data_root <- function() {
  root <- Sys.getenv("PALEO_DATA", unset = "")
  if (!nzchar(root)) {
    stop("PALEO_DATA is not set.\n",
         "  1) copy .Renviron.template to .Renviron\n",
         "  2) edit the path inside it\n",
         "  3) restart R (Session > Restart R)", call. = FALSE)
  }
  if (!dir.exists(root)) stop("PALEO_DATA folder not found: ", root, call. = FALSE)
  root
}

site_file_path <- function() {
  f <- Sys.getenv("PALEO_SITE_FILE", unset = "Nagpur Chandrapur Palaeolithic Site Raw Data.xlsx")
  path <- file.path(paleo_data_root(), f)
  if (!file.exists(path)) stop("Site file not found: ", path, call. = FALSE)
  path
}
