## =============================================================================
## utils.R — Small utility functions for ACH analysis
## =============================================================================

#' Standardise a numeric vector to [0, 1]
#' Returns all-zero vector if range is zero or non-finite.
std01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

#' Benjamini–Hochberg FDR correction
bh_adjust <- function(p) p.adjust(p, method = "BH")

#' Classify ACH result from adjusted rho and adjusted p-value
#' @param rho_adj  Sign-adjusted Spearman rho (margin distances already negated)
#' @param p_adj    BH-adjusted p-value
#' @param alpha    Significance threshold (default 0.05)
classify_ach <- function(rho_adj, p_adj, alpha = 0.05) {
  dplyr::case_when(
    is.na(rho_adj) | is.na(p_adj) ~ "Insufficient data",
    p_adj < alpha & rho_adj < 0   ~ "Supported",
    p_adj < alpha & rho_adj > 0   ~ "Opposite",
    TRUE                           ~ "Not significant"
  )
}

#' Convert wind direction (degrees) + speed to u/v components
calc_wind_uv <- function(dir_deg, spd) {
  rad <- dir_deg * pi / 180
  tibble::tibble(
    Wind_u = -spd * sin(rad),
    Wind_v = -spd * cos(rad)
  )
}

#' Significance label from p-value
get_sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

#' Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
