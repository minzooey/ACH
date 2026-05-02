## =============================================================================
## run_ACH.R — Abundant-Centre Hypothesis validation
##             Jangcheon Harbor dinoflagellate long-term abundance (411 days)
## =============================================================================
## Description
##   Tests the ACH for 32 dinoflagellate taxa using PCA-reduced environmental
##   space. Compares three niche centre definitions (CH, MVE, DMF) × three
##   distance metrics (Euclidean, Mahalanobis, Margin) = 9 settings (E1–E9).
##
## Inputs
##   JC_envs_daily.xlsx          — daily environmental variables
##   JC_abundance.xlsx           — daily raw abundance (sheet: Dino_daily)
##   dino_5day.xlsx              — 5-day smoothed abundance
##   5_phenology_variables.xlsx  — phenology output (DMF_Date per species/event)
##
## Outputs (written to OUT_DIR/)
##   dudi_loadings.csv
##   OMI_niche_parameters.csv
##   ACH_sp_results.csv
##   ACH_centre_comparison.csv         ← new: 3-centre comparison summary
##   ACH_model_data_long.csv
##   Fig1_species_niches_all_panels.pdf
##   Fig2_*.pdf
##   Fig3_*.pdf
##   Fig4_*/
##   Fig5_centre_comparison/           ← new
##
## Author      : [Author Name]
## Affiliation : [Institution]
## Last updated: 2025
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ── 0) Packages ──────────────────────────────────────────────────────────────
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr", "stringr", "lubridate",
  "ggplot2", "patchwork", "scales",
  "ade4", "MASS", "sp", "ellipse", "mgcv"
)

new_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

bayes_pkgs <- c("brms", "posterior")
bayes_ok   <- all(vapply(bayes_pkgs, requireNamespace, logical(1), quietly = TRUE))
if (bayes_ok) invisible(lapply(bayes_pkgs, library, character.only = TRUE))

## ── 1) Paths — edit INPUT_DIR to your project root ───────────────────────────
INPUT_DIR <- "."   # FIX W1: replaced hard-coded absolute path

ENVS_FILE  <- "data/JC_envs_daily.xlsx"
ABUN_FILE  <- "data/JC_abundance.xlsx"
SMOO_FILE  <- "data/dino_5day.xlsx"
PHENO_FILE <- "data/5_phenology_variables.xlsx"  # phenology output

OUT_DIR <- file.path(INPUT_DIR,
                     paste0("Output_ACH_", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

FUNC_DIR <- file.path(INPUT_DIR, "R")  # FIX W1: relative path

## ── 2) Load helper functions ──────────────────────────────────────────────────
source(file.path(FUNC_DIR, "utils.R"))
source(file.path(FUNC_DIR, "geometry.R"))
source(file.path(FUNC_DIR, "make_envelope.R"))
source(file.path(FUNC_DIR, "calc_all_dist.R"))
source(file.path(FUNC_DIR, "test_sp_ach.R"))
source(file.path(FUNC_DIR, "build_model_df.R"))
source(file.path(FUNC_DIR, "plot_sp_niche.R"))
source(file.path(FUNC_DIR, "fit_quad_lm_sp.R"))

## ── 3) Global settings ────────────────────────────────────────────────────────
min_occ <- 10L
alpha   <- 0.05
n_axes  <- 2L
set.seed(123)

# 9 settings: 3 centres × 3 distance metrics
# FIX B2 (original): settings_meta now matches the redesigned make_envelope /
#   calc_all_dist interfaces.
settings_meta <- tibble(
  setting   = paste0("E", 1:9),
  centre    = c(rep("CH", 3),  rep("MVE", 3),  rep("DMF", 3)),
  dist_type = rep(c("Euclidean", "Mahalanobis", "Margin"), 3),
  is_margin = rep(c(FALSE, FALSE, TRUE), 3)
)

ACH_cols <- c(
  "Supported"         = "#38D5BE",
  "Opposite"          = "#FF6F61",
  "Not significant"   = "#BDBDBD",
  "Insufficient data" = "#F0F0F0"
)

centre_cols <- c("CH" = "#FFD230", "MVE" = "#1A7595", "DMF" = "#E76F51")

## ── 4) Input data ─────────────────────────────────────────────────────────────
# Environmental variables
envs_raw <- openxlsx::read.xlsx(
  file.path(INPUT_DIR, ENVS_FILE), detectDates = TRUE
) %>%
  rename_with(~ sub("_.*", "", .x)) %>%
  mutate(
    Date     = lubridate::make_date(Year, Month, Day),
    DOY      = as.integer(format(Date, "%j")),   # FIX L8: integer, not character
    CountDay = seq_len(n()),
    NOX      = NO3 + NO2
  ) %>%
  bind_cols(calc_wind_uv(.$AverWD, .$AverWS)) %>%    # FIX: vectorised wind UV
  relocate(c(Date, DOY, CountDay, NOX), .after = Day) %>%
  dplyr::select(-c(MinTide, MaxTide, MaxWS, AverWS, AverWD)) %>%
  na.omit() %>%
  arrange(Date)                                        # FIX B5: no `by=`

# Smoothed dinoflagellate abundance
dino_sm5 <- openxlsx::read.xlsx(
  file.path(INPUT_DIR, SMOO_FILE), detectDates = TRUE
) %>%
  mutate(
    Date     = lubridate::make_date(Year, Month, Day),
    DOY      = as.integer(format(Date, "%j")),         # FIX L8
    CountDay = seq_len(n())
  ) %>%
  dplyr::filter(Date %in% envs_raw$Date) %>%
  dplyr::select(-any_of("Undefined")) %>%
  relocate(c(Date, DOY, CountDay), .after = Day) %>%
  na.omit()

# FIX B1: `dino_sm7` replaced everywhere with `dino_sm5`
# FIX L7: name-based column selection instead of hardcoded indices
SP_ALL  <- names(dino_sm5)[sapply(dino_sm5, is.numeric) &
                              !names(dino_sm5) %in% c("Year","Month","Day",
                                                       "DOY","CountDay")]
keep_sp <- names(which(
  colSums(dino_sm5[, SP_ALL] > 0, na.rm = TRUE) >= min_occ
))
dino_sm5 <- dino_sm5 %>%
  dplyr::select(Year, Month, Day, Date, DOY, CountDay, all_of(keep_sp))

# Environmental variable selection
env_vars <- c("Temperature", "Salinity", "NOX", "NH4", "PO4", "SiO2",
              "DLI", "SumPrec", "AverTide", "Wind_u", "Wind_v")
# Note: Tchla (from pigments) removed here since JC_pigments_daily.xlsx is
# not included in the repository. Add it back if the file is available:
#   pigs_raw <- openxlsx::read.xlsx(PIGS_FILE, detectDates = TRUE) %>%
#     mutate(Date = make_date(Year, Month, Day), Tchla = `Chl-a` + ...) %>%
#     dplyr::select(Date, Tchla)
#   envs_raw <- left_join(envs_raw, pigs_raw, by = "Date")
#   env_vars  <- c("Temperature", "Salinity", "Tchla", env_vars)

# Align sites
input <- envs_raw %>%
  dplyr::filter(Date %in% dino_sm5$Date) %>%
  left_join(dino_sm5, by = c("Year", "Month", "Day", "Date", "DOY", "CountDay")) %>%
  dplyr::filter(if_all(all_of(env_vars), is.finite)) %>%
  arrange(Date)

env   <- input %>% dplyr::select(all_of(env_vars))
abun  <- input %>% dplyr::select(all_of(keep_sp))

# Phenology output (for DMF centre)
pheno_raw <- tryCatch(
  openxlsx::read.xlsx(file.path(INPUT_DIR, PHENO_FILE), detectDates = TRUE),
  error = function(e) { warning("Phenology file not found — DMF centre skipped."); NULL }
)

# Standardise species name column and convert DMF_Date
pheno_df <- if (!is.null(pheno_raw)) {
  pheno_raw %>%
    dplyr::rename_with(~ stringr::str_to_title(.x), any_of(c("species","SPECIES"))) %>%
    dplyr::rename_with(~ "Species", any_of(c("species","Species"))) %>%
    dplyr::mutate(
      DMF_Date = as.Date(DMF_Date,
                         origin = if (is.numeric(DMF_Date)) "1899-12-30" else NULL)
    ) %>%
    dplyr::filter(!is.na(DMF_Date))
} else NULL

## ── 5) PCA ordination ────────────────────────────────────────────────────────
dudi <- ade4::dudi.pca(
  env[, env_vars],
  center = TRUE, scale = TRUE, scannf = FALSE, nf = n_axes
)

cat("PCA eigenvalues:\n"); print(round(dudi$eig, 4))

load_tab <- as.data.frame(dudi$c1[, seq_len(n_axes), drop = FALSE]) %>%
  tibble::rownames_to_column("Variable") %>%
  setNames(c("Variable", paste0("PC", seq_len(n_axes), "_loading")))
write.csv(load_tab, file.path(OUT_DIR, "dudi_loadings.csv"), row.names = FALSE)

## ── 6) Niche envelopes ────────────────────────────────────────────────────────
niche_res <- make_envelope(
  dudi       = dudi,
  abun       = abun,
  min_occ    = min_occ,
  n_axes     = n_axes,
  pheno_df   = pheno_df,
  site_dates = input$Date
)

site_xy      <- niche_res$site_xy
species_objs <- niche_res$species_objs
nic          <- niche_res$niche

# OMI niche parameters
if (!is.null(nic)) {
  niche_param_tab <- as.data.frame(ade4::niche.param(nic)) %>%
    tibble::rownames_to_column("Species")
  write.csv(niche_param_tab,
            file.path(OUT_DIR, "OMI_niche_parameters.csv"), row.names = FALSE)
}

# DMF coverage summary
dmf_summary <- purrr::map_dfr(names(species_objs), ~ tibble(
  Species  = .x,
  n_pres   = species_objs[[.x]]$n_pres,
  n_dmf    = species_objs[[.x]]$n_dmf,
  has_dmf  = !is.null(species_objs[[.x]]$dmf_c)
))
write.csv(dmf_summary, file.path(OUT_DIR, "DMF_centre_coverage.csv"), row.names = FALSE)
cat("\nDMF centre available for",
    sum(dmf_summary$has_dmf), "/", nrow(dmf_summary), "species.\n")

## ── 7) Distance calculation ───────────────────────────────────────────────────
dist_all <- purrr::map(
  species_objs,
  ~ calc_all_dist(.x, site_xy)
)

## ── 8) ACH test: Spearman + BH-FDR ──────────────────────────────────────────
ach_res <- purrr::map_dfr(
  names(species_objs),
  ~ test_sp_ach(species_objs[[.x]], dist_all[[.x]], settings_meta, min_n = min_occ)
) %>%
  group_by(setting) %>%
  mutate(
    p_adj = bh_adjust(pval),
    ACH   = classify_ach(rho_adj, p_adj, alpha = alpha)
  ) %>%
  ungroup() %>%
  left_join(settings_meta, by = "setting")

write.csv(ach_res, file.path(OUT_DIR, "ACH_sp_results.csv"), row.names = FALSE)

## ── 9) Centre comparison (new) ────────────────────────────────────────────────
# Compare ACH support rate and median rho across 3 centre definitions
centre_comparison <- ach_res %>%
  group_by(centre, dist_type) %>%
  summarise(
    n_species        = dplyr::n_distinct(species),
    n_valid          = sum(!is.na(rho_adj)),
    n_supported      = sum(ACH == "Supported", na.rm = TRUE),
    n_opposite       = sum(ACH == "Opposite",  na.rm = TRUE),
    pct_supported    = 100 * n_supported / n_valid,
    median_rho_adj   = median(rho_adj, na.rm = TRUE),
    mean_rho_adj     = mean(rho_adj,   na.rm = TRUE),
    p_wilcox         = tryCatch(
      wilcox.test(rho_adj[!is.na(rho_adj)], mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(wilcox_label = get_sig_label(p_wilcox)) %>%
  arrange(centre, dist_type)

write.csv(centre_comparison,
          file.path(OUT_DIR, "ACH_centre_comparison.csv"), row.names = FALSE)
print(centre_comparison)

## ── 10) Model-ready data ──────────────────────────────────────────────────────
model_df <- build_model_df(species_objs, dist_all, settings_meta,
                            min_occ = min_occ)
write.csv(model_df, file.path(OUT_DIR, "ACH_model_data_long.csv"), row.names = FALSE)

## ── 11) Visualisation ─────────────────────────────────────────────────────────
dist_lab_order <- c("Euclidean", "Mahalanobis", "Margin")

# ── Fig 1: Species niches ────────────────────────────────────────────────────
FIG1_DIR <- file.path(OUT_DIR, "Fig1_sp_niches")
dir.create(FIG1_DIR, showWarnings = FALSE)

fig1_plots <- purrr::map(
  names(species_objs),
  ~ plot_sp_niche(.x, species_objs[[.x]], site_xy)
)

p_fig1 <- patchwork::wrap_plots(fig1_plots,
                                 ncol   = 4,
                                 guides = "collect")
ggplot2::ggsave(
  file.path(OUT_DIR, "Fig1_species_niches_all_panels.pdf"),
  p_fig1,
  width  = 4.8 * 4,
  height = 4.5 * ceiling(length(fig1_plots) / 4)
)

# ── Fig 2-1: Violin plot of Spearman rho ────────────────────────────────────
p_2_1 <- ach_res %>%
  dplyr::mutate(
    setting = factor(setting, levels = settings_meta$setting),
    centre  = factor(centre,  levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = setting, y = rho_adj)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_violin(ggplot2::aes(fill = centre),
                       colour = "grey85", linewidth = 0.4, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    alpha = 0.7, size = 1.4
  ) +
  ggplot2::scale_fill_manual(values = centre_cols) +
  ggplot2::labs(
    x = NULL,
    y = "Spearman rho (margin sign reversed)",
    fill = "Centre",
    title = "Figure 2-1. Spearman correlation by setting"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    panel.grid = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(OUT_DIR, "Fig2_1_Spearman_violin.pdf"),
                p_2_1, width = 10, height = 5, dpi = 300)

# ── Fig 2-2: ACH result heatmap ──────────────────────────────────────────────
p_2_2 <- ach_res %>%
  dplyr::mutate(
    species = factor(species, levels = rev(sort(unique(species)))),
    setting = factor(setting, levels = settings_meta$setting),
    rho_lab = ifelse(!is.na(p_adj) & p_adj < alpha,
                     sprintf("%.2f", rho_adj), "")
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = setting, y = species, fill = ACH)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = rho_lab), size = 2.8) +
  ggplot2::scale_fill_manual(values = ACH_cols, drop = FALSE) +
  ggplot2::labs(
    x = NULL, y = NULL, fill = NULL,
    title = "Figure 2-2. ACH result heatmap (species × setting)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    panel.grid = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(OUT_DIR, "Fig2_2_ACH_heatmap.pdf"),
                p_2_2, width = 11, height = 10, dpi = 300)

# ── Fig 3: ACH support proportion by setting ─────────────────────────────────
prop_df <- ach_res %>%
  dplyr::count(setting, centre, ACH) %>%
  dplyr::group_by(setting) %>%
  dplyr::mutate(prop = n / sum(n) * 100) %>%
  dplyr::ungroup()

p_3 <- ggplot2::ggplot(
  prop_df,
  ggplot2::aes(x = factor(setting, levels = settings_meta$setting),
               y = prop, fill = ACH)
) +
  ggplot2::geom_col(width = 0.65, colour = "white") +
  ggplot2::facet_wrap(~ centre, scales = "free_x", nrow = 1) +
  ggplot2::scale_fill_manual(values = ACH_cols, drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 100),
                               expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of species (%)",
    fill = NULL,
    title = "Figure 3. ACH support proportion by centre × distance metric"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5),
    panel.grid  = ggplot2::element_blank(),
    strip.text  = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(file.path(OUT_DIR, "Fig3_ACH_support_proportion.pdf"),
                p_3, width = 10, height = 5, dpi = 300)

# ── Fig 5: Centre comparison (new) ───────────────────────────────────────────
FIG5_DIR <- file.path(OUT_DIR, "Fig5_centre_comparison")
dir.create(FIG5_DIR, showWarnings = FALSE)

# 5-1: Median rho by centre and distance metric
p_5_1 <- centre_comparison %>%
  dplyr::mutate(
    centre    = factor(centre,    levels = c("CH", "MVE", "DMF")),
    dist_type = factor(dist_type, levels = dist_lab_order)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = dist_type, y = median_rho_adj,
                                fill = centre, colour = centre)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.7),
                    width = 0.6, alpha = 0.85) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::geom_text(
    ggplot2::aes(label = wilcox_label,
                 y = median_rho_adj + ifelse(median_rho_adj >= 0, 0.02, -0.04)),
    position = ggplot2::position_dodge(0.7),
    size = 3.5, fontface = "bold", colour = "black"
  ) +
  ggplot2::scale_fill_manual(values   = centre_cols) +
  ggplot2::scale_colour_manual(values = centre_cols) +
  ggplot2::labs(
    x = "Distance metric", y = "Median Spearman rho",
    fill = "Centre", colour = "Centre",
    title = "Figure 5-1. ACH support by centre definition"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(FIG5_DIR, "Fig5_1_centre_median_rho.pdf"),
                p_5_1, width = 7, height = 5, dpi = 300)

# 5-2: Pairwise rho comparison across centres (Euclidean only for clarity)
rho_wide <- ach_res %>%
  dplyr::filter(dist_type == "Euclidean") %>%
  dplyr::select(species, centre, rho_adj) %>%
  tidyr::pivot_wider(names_from = centre, values_from = rho_adj)

if (all(c("CH","MVE","DMF") %in% names(rho_wide))) {
  p_5_2 <- rho_wide %>%
    ggplot2::ggplot(ggplot2::aes(x = CH, y = DMF)) +
    ggplot2::geom_point(ggplot2::aes(colour = MVE), size = 2.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = 2, colour = "grey40") +
    ggplot2::scale_colour_gradient2(
      low = "#d73027", mid = "grey80", high = "#1a9850", midpoint = 0,
      name = "MVE rho"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "CH centre rho",
      y = "DMF centre rho",
      title = "Figure 5-2. CH vs. DMF centre rho (Euclidean; colour = MVE rho)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

  ggplot2::ggsave(file.path(FIG5_DIR, "Fig5_2_centre_rho_scatter.pdf"),
                  p_5_2, width = 6, height = 5.5, dpi = 300)
}

write.csv(rho_wide, file.path(FIG5_DIR, "Fig5_rho_wide.csv"), row.names = FALSE)

# ── Fig 4: Predictive models ──────────────────────────────────────────────────
FIG4_DIR <- file.path(OUT_DIR, "Fig4_predictive_models")
dir.create(FIG4_DIR, showWarnings = FALSE)

# 4-1: GAM (distance–abundance, by setting)
fit_gam_setting <- function(df1) {
  mgcv::gam(
    log_abund ~ s(dist_signed, k = 4) + s(species, bs = "re"),
    data = df1, method = "REML"
  )
}

gam_fits    <- purrr::map(split(model_df, model_df$setting), fit_gam_setting)
gam_pred_df <- purrr::imap_dfr(gam_fits, function(fit, st) {
  dfi       <- dplyr::filter(model_df, setting == st)
  is_mg     <- unique(dfi$is_margin)
  nd <- tibble::tibble(
    species      = unique(dfi$species)[1],
    dist_std     = seq(0, 1, length.out = 100),
    dist_signed  = if (is_mg) -seq(0, 1, length.out = 100) else seq(0, 1, length.out = 100)
  )
  pr <- stats::predict(fit, newdata = nd, se.fit = TRUE)
  tibble::tibble(
    setting   = st,
    centre    = unique(dfi$centre),
    dist_type = unique(dfi$dist_type),
    dist_std  = nd$dist_std,
    y_med     = pr$fit,
    y_low     = pr$fit - 1.96 * pr$se.fit,
    y_high    = pr$fit + 1.96 * pr$se.fit
  )
})

p_4_1 <- gam_pred_df %>%
  dplyr::mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = dist_std, y = y_med,
                                colour = centre, fill = centre)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = y_low, ymax = y_high),
                       alpha = 0.18, linewidth = 0, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(~ dist_type, ncol = 3, scales = "free_y") +
  ggplot2::scale_colour_manual(values = centre_cols) +
  ggplot2::scale_fill_manual(values   = centre_cols) +
  ggplot2::labs(
    x = "Standardised distance (0–1)",
    y = "log(predicted abundance)",
    colour = "Centre", fill = "Centre",
    title = "Figure 4-1. GAM distance–abundance curves by centre"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(file.path(FIG4_DIR, "Fig4_1_GAM.pdf"),
                p_4_1, width = 11, height = 4.6, dpi = 300)

# 4-2: Species-specific quadratic LM curves
species_curve_df <- model_df %>%
  dplyr::group_by(setting, species) %>%
  dplyr::group_modify(~ fit_quad_lm_sp(.x)) %>%   # FIX B4: correct function name
  dplyr::ungroup()

p_4_2 <- species_curve_df %>%
  dplyr::mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = dist_std, y = pred,
                                group = species, colour = adj_r2)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
  ggplot2::facet_grid(centre ~ dist_type, scales = "free_y") +
  ggplot2::scale_colour_gradient2(
    low = "#d73027", mid = "grey70", high = "#1a9850", midpoint = 0.2,
    name = expression(Adj.~R^2)
  ) +
  ggplot2::labs(
    x = "Standardised distance (0–1)",
    y = "log(predicted abundance)",
    title = "Figure 4-2. Species-specific quadratic fitted curves"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(file.path(FIG4_DIR, "Fig4_2_species_curves.pdf"),
                p_4_2, width = 9, height = 8, dpi = 300)

# 4-3: Bayesian hierarchical model (optional — requires brms)
if (bayes_ok) {
  fit_brms_quad <- function(df1) {
    brms::brm(
      brms::bf(log_abund ~ 1 + dist_signed + I(dist_signed^2) +
                 (1 + dist_signed + I(dist_signed^2) | species)),
      data    = df1,
      family  = brms::gaussian(),
      chains  = 3, iter = 3000, warmup = 1500, cores = 3, seed = 123,
      refresh = 0,
      prior   = c(
        brms::prior(normal(0, 5),          class = "b"),
        brms::prior(normal(0, 5),          class = "Intercept"),
        brms::prior(student_t(3, 0, 2.5),  class = "sd"),
        brms::prior(lkj(2),                class = "cor"),
        brms::prior(student_t(3, 0, 2.5),  class = "sigma")
      )
    )
  }

  bayes_fits <- purrr::map(
    stats::setNames(settings_meta$setting, settings_meta$setting),
    function(st) {
      message("Bayesian model: ", st)
      fit_brms_quad(dplyr::filter(model_df, setting == st))
    }
  )

  bayes_pred_df <- purrr::imap_dfr(bayes_fits, function(fit, st) {
    dfi   <- dplyr::filter(model_df, setting == st)
    is_mg <- unique(dfi$is_margin)
    nd <- tibble::tibble(
      species     = unique(dfi$species)[1],
      dist_std    = seq(0, 1, length.out = 100),
      dist_signed = if (is_mg) -seq(0, 1, length.out = 100) else seq(0, 1, length.out = 100)
    )
    dr <- posterior::posterior_epred(fit, newdata = nd, re_formula = NA)
    q  <- apply(dr, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    tibble::tibble(
      setting   = st,
      centre    = unique(dfi$centre),
      dist_type = unique(dfi$dist_type),
      dist_std  = nd$dist_std,
      y_low     = q[1, ], y_med = q[2, ], y_high = q[3, ]
    )
  })

  p_4_3 <- bayes_pred_df %>%
    dplyr::mutate(
      dist_type = factor(dist_type, levels = dist_lab_order),
      centre    = factor(centre,    levels = c("CH","MVE","DMF"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = dist_std, y = y_med,
                                  colour = centre, fill = centre)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = y_low, ymax = y_high),
                         alpha = 0.18, linewidth = 0, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~ dist_type, ncol = 3, scales = "free_y") +
    ggplot2::scale_colour_manual(values = centre_cols) +
    ggplot2::scale_fill_manual(values   = centre_cols) +
    ggplot2::labs(
      x      = "Standardised distance (0–1)",
      y      = "log(predicted abundance)",
      colour = "Centre", fill = "Centre",
      title  = "Figure 4-3. Bayesian hierarchical quadratic model"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(file.path(FIG4_DIR, "Fig4_3_Bayesian.pdf"),
                  p_4_3, width = 11, height = 4.6, dpi = 300)
}

message("\n✓ All outputs written to: ", OUT_DIR)
