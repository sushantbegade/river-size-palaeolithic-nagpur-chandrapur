DOR_THRESHOLDS <- c(2, 5, 10, 25, 50)
DOR_PRIMARY    <- 5
MCM_PER_CMS_YR <- 31.536

dam_rows_for_dor <- function(bd) {
  tl <- as.character(bd$TIMELINE)
  which((bd$INSTREAM %in% "Instream") & !(tl %in% EXCL_TIMELINE) &
          !is.na(bd$CAP_MCM) & bd$CAP_MCM > 0 & !is.na(bd$HYRIV_ID))
}

cumulative_storage <- function(hr, bd, rows) {
  ids <- hr$HYRIV_ID
  nxt <- match(hr$NEXT_DOWN, ids)
  st  <- numeric(nrow(hr))
  at  <- match(bd$HYRIV_ID[rows], ids)
  ok  <- !is.na(at)
  if (any(ok)) {
    add <- tapply(bd$CAP_MCM[rows][ok], at[ok], sum)
    st[as.integer(names(add))] <- as.numeric(add)
  }
  for (i in order(hr$UPLAND_SKM)) {
    j <- nxt[i]
    if (!is.na(j)) st[j] <- st[j] + st[i]
  }
  st
}

dor_index <- function(hr, storage_mcm) {
  100 * storage_mcm / (pmax(hr$DIS_AV_CMS, 0.001) * MCM_PER_CMS_YR)
}

is_estimated_capacity <- function(bd, rows) {
  if (!"DATA_INFO" %in% names(bd)) return(rep(FALSE, length(rows)))
  info <- as.character(bd$DATA_INFO[rows])
  !is.na(info) & grepl("statistics|estimated", info, ignore.case = TRUE)
}