## =============================================================================
## fit_quad_lm_sp.R — Species-specific quadratic LM fitted curves
## =============================================================================

#' Fit quadratic LM and return prediction curve for one species × setting
#'
#' @param df1 Subset of model_df for one species × setting
#' @return Tibble with predicted values at 100 equally-spaced dist_std points,
#'         or NULL if fitting fails
fit_quad_lm_sp <- function(df1) {   # FIX B4: corrected function name
  fit <- tryCatch(
    stats::lm(log_abund ~ poly(dist_std, 2, raw = TRUE), data = df1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  newd <- tibble::tibble(dist_std = seq(0, 1, length.out = 100))
  pr   <- stats::predict(fit, newdata = newd)

  tibble::tibble(
    species   = unique(df1$species),
    setting   = unique(df1$setting),
    centre    = unique(df1$centre),
    dist_type = unique(df1$dist_type),
    dist_std  = newd$dist_std,
    pred      = pr,
    adj_r2    = summary(fit)$adj.r.squared
  )
}
