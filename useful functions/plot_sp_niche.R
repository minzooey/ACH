## =============================================================================
## plot_sp_niche.R — Per-species niche panel plot
## =============================================================================

#' Plot realized niche for one species (CH polygon + MVE ellipse + DMF point)
#'
#' @param sp       Species name (character)
#' @param env_obj  One element of species_objs (from make_envelope)
#' @param site_xy  Matrix (n_sites × 2) of all PCA site scores
#' @return ggplot object
plot_sp_niche <- function(sp, env_obj, site_xy) {
  pts_all <- tibble::as_tibble(site_xy, .name_repair = "minimal") %>%
    stats::setNames(c("Axis1", "Axis2"))

  pts_sp <- pts_all[env_obj$pres_mask, , drop = FALSE] %>%
    dplyr::mutate(
      abund     = env_obj$abund[env_obj$pres_mask],
      log_abund = log10(1 + abund)
    )

  hull_df <- tibble::as_tibble(env_obj$ch_v, .name_repair = "minimal") %>%
    stats::setNames(c("Axis1", "Axis2"))

  ell_df <- ellipse_df(env_obj$mve_c, env_obj$mve_cov, level = 0.975)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey90") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey90") +
    ggplot2::geom_point(data = pts_all,
                        ggplot2::aes(Axis1, Axis2),
                        colour = "grey70", size = 0.5, alpha = 0.8) +
    ggplot2::geom_polygon(data = hull_df,
                          ggplot2::aes(Axis1, Axis2),
                          fill    = scales::alpha("#FFD230", 0.20),
                          colour  = "#FFD230", linewidth = 0.8) +
    ggplot2::geom_path(data = ell_df,
                       ggplot2::aes(Axis1, Axis2),
                       colour = "#1A7595", linewidth = 0.9) +
    ggplot2::geom_point(data = pts_sp,
                        ggplot2::aes(Axis1, Axis2, size = log_abund),
                        colour = "black", alpha = 0.9) +
    # CH centroid
    ggplot2::geom_point(ggplot2::aes(x = env_obj$ch_c[1], y = env_obj$ch_c[2]),
                        shape = 4, size = 4, stroke = 1.2, colour = "#FFD230") +
    # MVE centroid
    ggplot2::geom_point(ggplot2::aes(x = env_obj$mve_c[1], y = env_obj$mve_c[2]),
                        shape = 4, size = 4, stroke = 1.2, colour = "#1A7595") +
    # DMF centroid (only if available)
    {
      if (!is.null(env_obj$dmf_c)) {
        ggplot2::geom_point(
          ggplot2::aes(x = env_obj$dmf_c[1], y = env_obj$dmf_c[2]),
          shape = 8, size = 4, stroke = 1.2, colour = "#E76F51"
        )
      }
    } +
    ggplot2::scale_size_continuous(
      name = expression(log[10](1 + abundance))
    ) +
    ggplot2::labs(title = sp, x = "Axis1", y = "Axis2") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "none"
    )

  p
}
