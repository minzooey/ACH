## =============================================================================
## build_model_df.R — Build long-format data frame for predictive models
## =============================================================================

#' Assemble long-format model data for all species × settings
#'
#' @param species_objs Named list of species envelopes (from make_envelope)
#' @param dist_all     Named list (species) of named lists (settings) of
#'                     distance vectors (from calc_all_dist)
#' @param settings_meta Tibble with columns: setting, centre, dist_type, is_margin
#' @param min_occ      Minimum valid observations to include a species×setting
#' @return Long tibble with columns:
#'   species, setting, centre, dist_type, is_margin,
#'   dist_std, dist_signed, log_abund, CountDay
#'   (CountDay: day-index from make_envelope()$species_objs[[sp]]$count_day,
#'    required by the AR(1) Bayesian model's `ar(time = CountDay, gr = species, p = 1)`
#'    term in run_ACH.R.)
build_model_df <- function(species_objs, dist_all, settings_meta,
                           min_occ = 5L) {
  purrr::map_dfr(names(species_objs), function(sp) {
    env_obj <- species_objs[[sp]]
    purrr::map_dfr(settings_meta$setting, function(st) {
      is_mg <- settings_meta$is_margin[match(st, settings_meta$setting)]
      d_full <- dist_all[[sp]][[st]]
      if (is.null(d_full)) return(NULL)
      
      x  <- d_full[env_obj$pres_mask]
      y  <- env_obj$log_abund[env_obj$pres_mask]
      cd_full <- env_obj$count_day
      if (is.null(cd_full)) {
        warning("species_objs[['", sp, "']] has no $count_day (older make_envelope() ",
                "output?) -> falling back to row index for CountDay.")
        cd_full <- seq_along(env_obj$abund)
      }
      cd <- cd_full[env_obj$pres_mask]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < min_occ) return(NULL)
      
      tibble::tibble(
        species     = sp,
        setting     = st,
        centre      = settings_meta$centre[match(st, settings_meta$setting)],
        dist_type   = settings_meta$dist_type[match(st, settings_meta$setting)],
        is_margin   = is_mg,
        dist_std    = std01(x[ok]),
        dist_signed = if (is_mg) -std01(x[ok]) else std01(x[ok]),
        log_abund   = y[ok],
        CountDay    = cd[ok]
      ) %>%
        dplyr::mutate(species = factor(species))
    })
  })
}