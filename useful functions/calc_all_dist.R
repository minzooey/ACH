## =============================================================================
## calc_all_dist.R — Compute all 9 distance settings for one species
##
## Settings
##   Centre   E1  E2  E3    E4   E5   E6     E7   E8   E9
##   -------  --  --  --    --   --   --     --   --   --
##   CH       ●   ●   ●
##   MVE                    ●    ●    ●
##   DMF                               ●    ●    ●
##   Metric   Eu  Ma  Mg    Eu   Ma   Mg    Eu   Ma   Mg
##
##   Eu = Euclidean, Ma = Mahalanobis, Mg = Margin
##
## Margin definitions
##   E3 (CH-Margin)  : interior distance to CH polygon boundary via p2seg_dist()
##                     — FIX L1: replaces vertex-only approximation
##   E6 (MVE-Margin) : thr − d_maha, where thr = 97.5th pct of presence Maha
##   E9 (DMF-Margin) : interior distance to DMF-point CH boundary (≥3 DMF pts)
##                     returns NA vector when < 3 DMF points available
##
## Note: all distances are for ALL sites (not just presences).
##       Subsetting to presence sites is done in test_sp_ach().
## =============================================================================

#' Compute the 9 distance vectors for one species object
#'
#' @param sp_obj   One element of the `species_objs` list from make_envelope()
#' @param site_xy  Matrix (n_sites × n_axes) — all PCA site coordinates
#' @return Named list E1…E9, each a numeric vector of length n_sites.
#'         DMF settings (E7–E9) are NA vectors when dmf_c is NULL.
calc_all_dist <- function(sp_obj, site_xy) {

  pts  <- as.matrix(site_xy)
  n    <- nrow(pts)

  ## ── Helper: CH boundary margin via p2seg_dist ────────────────────────────
  ch_margin <- function(query_pts, ch_vertices) {
    poly_closed <- rbind(ch_vertices, ch_vertices[1L, , drop = FALSE])
    n_edges <- nrow(poly_closed) - 1L

    inside <- sp::point.in.polygon(
      query_pts[, 1], query_pts[, 2],
      ch_vertices[, 1], ch_vertices[, 2]
    )

    vapply(seq_len(nrow(query_pts)), function(i) {
      if (inside[i] == 0L) return(0)
      p <- query_pts[i, ]
      min(vapply(seq_len(n_edges), function(j) {
        p2seg_dist(p, poly_closed[j, ], poly_closed[j + 1, ])
      }, numeric(1)))
    }, numeric(1))
  }

  ## ── E1–E3: CH centre ─────────────────────────────────────────────────────
  e1 <- eucl_dist(pts, sp_obj$ch_c)
  e2 <- maha_dist(pts, sp_obj$ch_c, sp_obj$sp_cov)
  e3 <- ch_margin(pts, sp_obj$ch_v)

  ## ── E4–E6: MVE centre ────────────────────────────────────────────────────
  e4 <- eucl_dist(pts, sp_obj$mve_c)
  e5 <- maha_dist(pts, sp_obj$mve_c, sp_obj$mve_cov)
  e6 <- margin_mve_dist(pts, sp_obj$mve_c, sp_obj$mve_cov, sp_obj$sp_pts)

  ## ── E7–E9: DMF centre ────────────────────────────────────────────────────
  na_vec <- rep(NA_real_, n)

  if (is.null(sp_obj$dmf_c)) {
    e7 <- e8 <- e9 <- na_vec
  } else {
    e7 <- eucl_dist(pts, sp_obj$dmf_c)
    e8 <- maha_dist(pts, sp_obj$dmf_c, sp_obj$dmf_cov)

    # E9: DMF-CH margin — requires ≥3 DMF points to form a polygon
    e9 <- if (!is.null(sp_obj$dmf_pts) && nrow(sp_obj$dmf_pts) >= 3L) {
      dmf_ch_idx <- chull(sp_obj$dmf_pts)
      dmf_ch_v   <- sp_obj$dmf_pts[dmf_ch_idx, , drop = FALSE]
      ch_margin(pts, dmf_ch_v)
    } else {
      na_vec
    }
  }

  list(E1 = e1, E2 = e2, E3 = e3,
       E4 = e4, E5 = e5, E6 = e6,
       E7 = e7, E8 = e8, E9 = e9)
}
