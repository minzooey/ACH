## =============================================================================
## make_envelope.R — Build species niche envelopes in PCA ordination space
##
## Redesign notes (vs. original):
##   • Single entry-point replaces the fragmented original structure
##   • chull() used for 2D CH (replaces geometry::convhulln which returns edge
##     indices causing vertex duplication and centroid bias)
##   • MASS::cov.mve() for minimum-volume ellipsoid
##   • ade4::niche() called internally for OMI parameters
##   • DMF centre computed from phenology output (new feature)
##   • No library() calls inside the function body
##   • omi_input = "abundance" (default, unchanged) | "presence": lets the
##     ade4::niche() call (marginality/tolerance) run on a binary 0/1
##     occurrence matrix instead of raw abundance, so OMI niche parameters
##     can be estimated independently of abundance weighting. CH/MVE/DMF
##     centroids are unaffected either way (already presence-based).
##   • species_objs[[sp]]$count_day: full-length day-index vector (aligned
##     with abund_vec/pres_mask), derived from site_dates. Needed downstream
##     by build_model_df.R for the `CountDay` column used in the AR(1)
##     Bayesian model's `ar(time = CountDay, gr = species, p = 1)` term.
## =============================================================================

#' Build per-species niche envelopes in PCA ordination space
#'
#' @param dudi       `dudi.pca` object (from ade4)
#' @param abun       Data frame: rows = sites, columns = species (raw abundance)
#' @param min_occ    Minimum number of presences required to include a species
#' @param n_axes     Number of PCA axes to use (default 2)
#' @param pheno_df   Optional: phenology data frame with columns
#'                   `Species`, `DMF_Date` (Date class).
#'                   If NULL, DMF centre is skipped.
#' @param site_dates Optional: Date vector of length nrow(abun) used to match
#'                   DMF dates to site rows. Required when pheno_df is supplied.
#' @param count_day  Optional: numeric/integer vector of length nrow(abun),
#'                   same row order as abun/site_xy. Passed through unchanged
#'                   to species_objs[[sp]]$count_day for use as the AR(1)
#'                   model's `time=` index (see build_model_df.R). Pass the
#'                   same `CountDay` column already used elsewhere in
#'                   run_ACH.R (input$CountDay <- seq_len(n())) so the
#'                   definition of "CountDay" stays consistent across the
#'                   whole pipeline. If NULL, falls back to site_dates-based
#'                   calendar-day arithmetic, and finally to a plain row
#'                   index (1:n) if site_dates is also NULL.
#' @param omi_input  Either "abundance" (default) or "presence". Controls only
#'                   the input passed to `ade4::niche()` for the OMI niche
#'                   parameters (marginality, tolerance). "abundance" keeps
#'                   the original behaviour (abundance-weighted `niche()`
#'                   call). "presence" converts `abun` to a binary 0/1
#'                   occurrence matrix immediately before the `niche()` call,
#'                   so marginality/tolerance are estimated from occurrence
#'                   alone. This does NOT affect CH/MVE/DMF centroid
#'                   calculation below, which already uses presence-only
#'                   points (`pres_mask`) regardless of this option.
#'
#' @return Named list:
#'   \item{site_xy}{Matrix (n × n_axes) of PCA site scores}
#'   \item{global_cov}{Global covariance of site_xy}
#'   \item{species_objs}{Named list, one element per qualifying species}
#'   \item{niche}{ade4 niche object for OMI parameter extraction}
#'   \item{omi_input}{Echoes back which input type ("abundance"/"presence")
#'                     was used for the niche()/niche.param() call above,
#'                     for traceability in downstream outputs.}
make_envelope <- function(dudi,
                          abun,
                          min_occ    = 5L,
                          n_axes     = 2L,
                          pheno_df   = NULL,
                          site_dates = NULL,
                          count_day  = NULL,
                          omi_input  = c("abundance", "presence")) {
  
  omi_input <- match.arg(omi_input)
  
  ## ── 0) Validate inputs ───────────────────────────────────────────────────
  stopifnot(
    inherits(dudi, "dudi"),
    is.data.frame(abun),
    nrow(dudi$li) == nrow(abun)
  )
  if (!is.null(pheno_df) && is.null(site_dates))
    stop("site_dates must be provided when pheno_df is supplied.")
  
  ## ── 1) Site coordinates ──────────────────────────────────────────────────
  ax_names <- colnames(dudi$li)[seq_len(n_axes)]
  site_xy  <- as.matrix(dudi$li[, ax_names, drop = FALSE])
  global_cov <- stats::cov(site_xy)
  
  ## ── 2) OMI niche (for niche.param export) ────────────────────────────────
  # Filter to species meeting min_occ before passing to ade4::niche
  keep_sp <- names(which(
    colSums(abun > 0, na.rm = TRUE) >= min_occ
  ))
  abun <- abun[, keep_sp, drop = FALSE]
  
  # `niche_input`: what actually goes into ade4::niche(). Kept separate from
  # `abun` because `abun` (raw abundance) is still needed below for each
  # species' abund_vec / log_abund, regardless of omi_input.
  niche_input <- if (omi_input == "presence") {
    as.data.frame((abun > 0) * 1)   # binary 0/1 occurrence matrix, same dims/names as abun
  } else {
    abun
  }
  
  niche_obj <- tryCatch(
    ade4::niche(dudi, niche_input, scannf = FALSE, nf = n_axes),
    error = function(e) {
      warning("ade4::niche() failed: ", conditionMessage(e))
      NULL
    }
  )
  
  ## ── 2b) CountDay (day index, needed for AR(1) Bayesian model's `time=` arg) ──
  # Full-length vector aligned with abun's rows (same order as site_xy).
  # Prefer the caller-supplied `count_day` so the definition matches
  # whatever is already used elsewhere in the pipeline (in run_ACH.R,
  # input$CountDay <- seq_len(n()), i.e. a row index, NOT a calendar-day
  # gap-aware count). Only derive from site_dates / row index as fallbacks.
  if (!is.null(count_day)) {
    stopifnot(length(count_day) == nrow(abun))
    count_day_full <- count_day
  } else if (!is.null(site_dates)) {
    warning("make_envelope(): count_day not supplied -> deriving CountDay from ",
            "site_dates (calendar-day arithmetic). This may NOT match the ",
            "CountDay convention used elsewhere in run_ACH.R (seq_len(n())) ",
            "if there are date gaps -- pass count_day explicitly to avoid ambiguity.")
    count_day_full <- as.numeric(site_dates - min(site_dates, na.rm = TRUE)) + 1
  } else {
    warning("make_envelope(): neither count_day nor site_dates supplied -> ",
            "CountDay falls back to row index (1:n).")
    count_day_full <- seq_len(nrow(abun))
  }
  
  ## ── 3) Per-species envelopes ─────────────────────────────────────────────
  species_objs <- lapply(keep_sp, function(sp) {
    abund_vec <- abun[[sp]]
    pres_mask <- !is.na(abund_vec) & abund_vec > 0
    if (sum(pres_mask) < min_occ) return(NULL)
    
    pts <- site_xy[pres_mask, , drop = FALSE]
    
    ## — Convex Hull (2D) using base chull() — FIX L2
    ch_idx <- chull(pts)            # indices into pts (unique vertices)
    ch_v   <- pts[ch_idx, , drop = FALSE]   # CH vertex matrix
    ch_c   <- colMeans(ch_v)        # CH centroid (unbiased mean of vertices)
    
    ## — Minimum Volume Ellipsoid
    n_mve <- max(floor(0.975 * nrow(pts)), ncol(pts) + 1L)
    mve   <- tryCatch(
      MASS::cov.mve(pts, quantile.used = n_mve),
      error = function(e) list(center = colMeans(pts), cov = stats::cov(pts))
    )
    mve_c   <- mve$center
    mve_cov <- mve$cov
    
    ## — DMF centre (new) ────────────────────────────────────────────────────
    dmf_c   <- NULL
    dmf_cov <- NULL
    dmf_pts <- NULL
    
    if (!is.null(pheno_df) && sp %in% pheno_df$Species) {
      ph_sp     <- dplyr::filter(pheno_df, Species == sp)
      dmf_dates <- unique(stats::na.omit(ph_sp$DMF_Date))
      
      if (length(dmf_dates) > 0) {
        dmf_idx <- which(site_dates %in% dmf_dates)
        
        if (length(dmf_idx) > 0) {
          dmf_pts <- site_xy[dmf_idx, , drop = FALSE]
          dmf_c   <- colMeans(dmf_pts)   # centroid of all DMF env. positions
          # Covariance: use species cov when < 3 DMF points (singular otherwise)
          dmf_cov <- if (nrow(dmf_pts) >= 3) {
            stats::cov(dmf_pts)
          } else {
            stats::cov(pts)             # fallback: species-wide covariance
          }
        }
      }
    }
    
    list(
      sp        = sp,
      abund     = abund_vec,
      log_abund = log10(abund_vec + 1),
      count_day = count_day_full,   # full-length, aligned with abund_vec/pres_mask
      pres_mask = pres_mask,
      sp_pts    = pts,          # presence-only PCA coordinates
      ch_v      = ch_v,         # CH vertex matrix  (k × n_axes)
      ch_c      = ch_c,         # CH centroid
      mve_c     = mve_c,        # MVE centre
      mve_cov   = mve_cov,      # MVE covariance
      sp_cov    = stats::cov(pts),  # species-wide covariance (for CH Maha)
      dmf_c     = dmf_c,        # DMF centroid (NULL if unavailable)
      dmf_cov   = dmf_cov,      # Covariance for DMF Mahalanobis
      dmf_pts   = dmf_pts,      # Raw DMF env. points (for DMF margin)
      n_pres    = sum(pres_mask),
      n_dmf     = if (!is.null(dmf_pts)) nrow(dmf_pts) else 0L
    )
  })
  
  names(species_objs) <- keep_sp
  species_objs        <- Filter(Negate(is.null), species_objs)
  
  # niche.param() must be called HERE while abun still exists in this
  # environment. Returning the niche object alone is insufficient: ade4::niche()
  # captures a reference to the local frame, which is destroyed on return.
  niche_param <- if (!is.null(niche_obj)) {
    tryCatch(
      as.data.frame(ade4::niche.param(niche_obj)) %>%
        tibble::rownames_to_column("Species"),
      error = function(e) {
        warning("niche.param() failed: ", conditionMessage(e)); NULL
      }
    )
  } else NULL
  
  list(
    site_xy      = site_xy,
    global_cov   = global_cov,
    species_objs = species_objs,
    niche        = niche_obj,
    niche_param  = niche_param,  # pre-computed — use directly, not niche.param()
    omi_input    = omi_input     # "abundance" or "presence" — which input fed ade4::niche()
  )
}