## =============================================================================
## make_envelope.R — Build species niche envelopes in PCA ordination space
##
## Redesign notes (vs. original):
##   • Single entry-point replaces the fragmented original structure
##   • chull() used for 2D CH (replaces geometry::convhulln which returns edge
##     indices causing vertex duplication and centroid bias — FIX L2)
##   • MASS::cov.mve() for minimum-volume ellipsoid
##   • ade4::niche() called internally for OMI parameters
##   • DMF centre computed from phenology output (new feature)
##   • No library() calls inside the function body (FIX L3)
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
#'
#' @return Named list:
#'   \item{site_xy}{Matrix (n × n_axes) of PCA site scores}
#'   \item{global_cov}{Global covariance of site_xy}
#'   \item{species_objs}{Named list, one element per qualifying species}
#'   \item{niche}{ade4 niche object for OMI parameter extraction}
make_envelope <- function(dudi,
                          abun,
                          min_occ    = 10L,
                          n_axes     = 2L,
                          pheno_df   = NULL,
                          site_dates = NULL) {
  
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
  abun_filt <- abun[, keep_sp, drop = FALSE]
  
  niche_obj <- tryCatch(
    ade4::niche(dudi, abun_filt, scannf = FALSE, nf = n_axes),
    error = function(e) {
      warning("ade4::niche() failed: ", conditionMessage(e))
      NULL
    }
  )
  
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
  
  # niche.param() must be called HERE while abun_filt still exists in this
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
    niche_param  = niche_param   # pre-computed — use directly, not niche.param()
  )
}