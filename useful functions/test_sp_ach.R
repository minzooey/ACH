## =============================================================================
## test_sp_ach.R — Spearman ACH test for one species across all settings
## =============================================================================

#' Run Spearman correlation (distance vs. log-abundance) for one species
#'
#' @param sp_obj       One element of species_objs (from make_envelope)
#' @param dist_list    One element of dist_all (from calc_all_dist),
#'                     named list E1…E9
#' @param settings_meta Tibble with columns: setting, centre, dist_type, is_margin
#' @param min_n        Minimum number of valid presence observations
#' @return Tibble with one row per setting
test_sp_ach <- function(sp_obj, dist_list, settings_meta, min_n = 10L) {
  purrr::map_dfr(settings_meta$setting, function(st) {
    d_full <- dist_list[[st]]
    if (is.null(d_full)) return(tibble::tibble(
      species = sp_obj$sp, setting = st, n = 0L,
      rho = NA_real_, rho_adj = NA_real_, pval = NA_real_
    ))

    x  <- d_full[sp_obj$pres_mask]
    y  <- sp_obj$log_abund[sp_obj$pres_mask]
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)

    if (n_ok < min_n) return(tibble::tibble(
      species = sp_obj$sp, setting = st, n = n_ok,
      rho = NA_real_, rho_adj = NA_real_, pval = NA_real_
    ))

    ct  <- suppressWarnings(
      cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    )
    rho <- unname(ct$estimate)
    is_mg <- settings_meta$is_margin[match(st, settings_meta$setting)]

    tibble::tibble(
      species = sp_obj$sp,
      setting = st,
      n       = n_ok,
      rho     = rho,
      rho_adj = if (is_mg) -rho else rho,   # margin: invert sign for ACH
      pval    = ct$p.value
    )
  })
}
