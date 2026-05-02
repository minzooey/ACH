## =============================================================================
## geometry.R — Distance and geometric primitives for niche envelope analysis
## =============================================================================

#' Euclidean distance from each row of pts to centre vector
#' @param pts   Numeric matrix (n × d)
#' @param ctr   Numeric vector of length d
eucl_dist <- function(pts, ctr) {
  sqrt(rowSums(sweep(as.matrix(pts), 2, ctr, "-")^2))
}

#' Mahalanobis distance from each row of pts to centre
#' Falls back to Moore–Penrose pseudo-inverse if cov_mat is singular.
#' @param pts     Numeric matrix (n × d)
#' @param ctr     Numeric vector (length d)
#' @param cov_mat Covariance matrix (d × d)
maha_dist <- function(pts, ctr, cov_mat) {
  inv_cov <- tryCatch(
    solve(cov_mat),
    error = function(e) MASS::ginv(cov_mat)
  )
  dx <- sweep(as.matrix(pts), 2, ctr, "-")
  sqrt(pmax(0, rowSums((dx %*% inv_cov) * dx)))
}

#' Shortest distance from point p to line segment [a, b]
#' @param p  Numeric vector (length 2)
#' @param a  Segment start (length 2)
#' @param b  Segment end   (length 2)
p2seg_dist <- function(p, a, b) {
  ab <- b - a
  t0 <- max(0, min(1, sum((p - a) * ab) / (sum(ab^2) + 1e-12)))
  sqrt(sum((p - (a + t0 * ab))^2))
}

#' Shortest distance from each row of pts to the convex-hull polygon boundary
#'
#' Points outside the polygon receive distance 0 (not penalised).
#' Points inside receive the minimum perpendicular distance to any edge.
#' Uses p2seg_dist() so the margin is exact (fixes the vertex-only approximation).
#'
#' @param pts        Numeric matrix (n × 2) of query points
#' @param polygon_xy Numeric matrix (k × 2) of CH vertices (open polygon,
#'                   does NOT repeat the first vertex at the end)
margin_nr_dist <- function(pts, polygon_xy) {
  poly_closed <- rbind(polygon_xy, polygon_xy[1, , drop = FALSE])
  n_edges <- nrow(poly_closed) - 1L

  inside_flag <- sp::point.in.polygon(
    point.x = pts[, 1], point.y = pts[, 2],
    pol.x   = polygon_xy[, 1], pol.y = polygon_xy[, 2]
  )

  vapply(seq_len(nrow(pts)), function(i) {
    if (inside_flag[i] == 0L) return(0)
    p <- pts[i, ]
    min(vapply(seq_len(n_edges), function(j) {
      p2seg_dist(p, poly_closed[j, ], poly_closed[j + 1, ])
    }, numeric(1)))
  }, numeric(1))
}

#' Margin distance inside an MVE ellipse
#'
#' Value is 0 outside the ellipse; positive (and largest at centre) inside.
#' Defined as: thr − d_maha, where thr is the 97.5th percentile Mahalanobis
#' distance of the reference (presence) points.
#'
#' @param pts     Query points (n × 2)
#' @param center  MVE centre
#' @param cov_mat MVE covariance
#' @param ref_pts Presence points used to calibrate the ellipse radius
margin_mve_dist <- function(pts, center, cov_mat, ref_pts) {
  d_all <- maha_dist(pts,     center, cov_mat)
  d_ref <- maha_dist(ref_pts, center, cov_mat)
  thr   <- stats::quantile(d_ref, probs = 0.975, na.rm = TRUE)
  pmax(0, thr - d_all)
}

#' Build a data frame of ellipse coordinates for plotting
#'
#' @param center  Centre vector (length 2)
#' @param cov_mat Covariance matrix (2 × 2)
#' @param level   Probability level passed to ellipse::ellipse()
#' @param n       Number of points on the ellipse perimeter
ellipse_df <- function(center, cov_mat, level = 0.975, n = 200) {
  # FIX L4: removed the erroneous `if (as.numeric(level))` branch —
  # ellipse::ellipse() already accepts numeric level directly.
  ee <- ellipse::ellipse(cov_mat, centre = center, level = level, npoints = n)
  tibble::as_tibble(ee, .name_repair = "minimal") %>%
    stats::setNames(c("Axis1", "Axis2"))
}
