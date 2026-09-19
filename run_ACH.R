## =============================================================================
## run_ACH.R — Abundant-Centre Hypothesis validation
##             Jangcheon Harbor dinoflagellate long-term abundance (411 days)
## =============================================================================
## Description
##   Tests the ACH for 32 dinoflagellate taxa using PCA-reduced environmental
##   space. Compares three niche centre definitions (CH, MVE, DMF) × two
##   distance metrics (Euclidean, Mahalanobis) = 6 settings (E1–E6).
##
##   Temporal autocorrelation in the 411-day time series is mitigated with
##   three complementary approaches (nominal results are always retained
##   alongside the corrected ones, never overwritten):
##     1. Moving block bootstrap CI for Spearman's rho          -> ACH_boot
##     2. Effective-N correction (Pyper & Peterman 1998)        -> ACH_effN
##     3. AR(1) residual structure in the Bayesian hierarchical
##        quadratic model (brms)                                -> ACH_AR1
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
##   ACH_sp_results.csv               — original / nominal results (unchanged)
##   ACH_sp_results_effN.csv          — + effective-N corrected p-values
##   ACH_sp_results_corrected.csv     — ACH, ACH_boot, ACH_effN, ACH_AR1 side by side
##   ACH_correction_comparison.csv    — agreement summary across the 3 methods
##   ACH_centre_comparison.csv
##   ACH_model_data_long.csv
##   Fig1_species_niches_all_panels.pdf
##   Fig2_*.pdf / Fig2_*_corrected.pdf
##   Fig3_*.pdf / Fig3_*_corrected.pdf
##   Fig4_*/
##   Fig5_centre_comparison/
##
## Author      : Juhee Min
## Affiliation : Department of Oceanography, Chonnam National University
## Last updated: 30 July 2026
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
INPUT_DIR <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology"

ENVS_FILE  <- "JC_envs_daily.xlsx"
PIGS_FILE  <- "JC_pigments_daily.xlsx"
ABUN_FILE  <- "JC_abundance.xlsx"
SMOO_FILE  <- "dino_5day.xlsx"
PHENO_FILE <- "ONE/Output_260801/5_phenology_variables.xlsx"  # phenology output

OUT_DIR <- file.path(INPUT_DIR,
                     paste0("ACH/Output_", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

FUNC_DIR <- file.path(INPUT_DIR, "/ACH/useful functions")

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
min_occ <- 5L
alpha   <- 0.05
n_axes  <- 2L
set.seed(123)

settings_meta <- tibble(
  setting   = paste0("E", 1:9),
  centre    = c(rep("CH", 3),  rep("MVE", 3),  rep("DMF", 3)),
  dist_type = rep(c("Euclidean", "Mahalanobis", "Margin"), 3),
  is_margin = rep(c(rep(FALSE, 2), TRUE),3)
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
tchla <- openxlsx::read.xlsx(file.path(INPUT_DIR, PIGS_FILE)) %>%
  mutate(TChla = `Chlide-a` + `Chl-a` + `DVChl-a`,
         Date  = lubridate::make_date(Year, Month, Day),) %>%
  dplyr::select(Date, TChla)

envs_raw <- openxlsx::read.xlsx(
  file.path(INPUT_DIR, ENVS_FILE), detectDates = TRUE
) %>%
  rename_with(~ sub("_.*", "", .x)) %>%
  mutate(
    Date     = lubridate::make_date(Year, Month, Day),
    DOY      = as.integer(format(Date, "%j")),
    CountDay = seq_len(n()),
    NOX      = NO3 + NO2
  ) %>%
  left_join(tchla, by = "Date") %>%
  bind_cols(calc_wind_uv(.$AverWD, .$AverWS)) %>%
  relocate(c(Date, DOY, CountDay), .after = Day) %>%
  dplyr::select(-c(Fluorescence, MinTide, MaxTide, MaxWS, AverWS, AverWD)) %>%
  na.omit() %>%
  arrange(Date)

# Smoothed dinoflagellate abundance
dino_abun <- openxlsx::read.xlsx(
  file.path(INPUT_DIR, ABUN_FILE), sheet = "Dino_daily", detectDates = TRUE
) %>%
  mutate(
    Date     = lubridate::make_date(Year, Month, Day),
    DOY      = as.integer(format(Date, "%j")),
    CountDay = seq_len(n())
  ) %>%
  dplyr::filter(Date %in% envs_raw$Date) %>%
  dplyr::select(-c("Undefined", "Meso")) %>%
  relocate(c(Date, DOY, CountDay), .after = Day) %>%
  na.omit()

SP_ALL  <- names(dino_abun)[sapply(dino_abun, is.numeric) &
                              !names(dino_abun) %in% c("Year","Month","Day", "DOY","CountDay")]
keep_sp <- names(which(
  colSums(dino_abun[, SP_ALL] > 0, na.rm = TRUE) >= min_occ
))

dino_abun <- dino_abun %>%
  dplyr::select(Year, Month, Day, Date, DOY, CountDay, all_of(keep_sp))

# Environmental variable selection
env_vars <- envs_raw[, 7:20]

vif_res <- usdm::vifstep(env_vars, th = 5)
print(vif_res)

cor_table <- cor(env_vars, use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation") %>%
  filter(Var1 != Var2)

# final select
env_vars <- c("Temperature", "Salinity", "TChla", "NH4", "PO4", "SiO2",
              "DLI", "SumPrec", "AverTide", "Wind_u", "Wind_v")

# Align sites
input <- envs_raw %>%
  dplyr::filter(Date %in% dino_abun$Date) %>%
  left_join(dino_abun, by = c("Year", "Month", "Day", "Date", "DOY", "CountDay")) %>% # dino_occ or dino_sm5
  dplyr::filter(if_all(all_of(env_vars), is.finite)) %>%
  arrange(Date)

env  <- input %>% dplyr::select(all_of(env_vars))
abun <- input %>% dplyr::select(all_of(keep_sp))

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
  site_dates = input$Date,
  count_day  = input$CountDay,
  omi_input  = "abundance"   # unchanged default; drives all downstream distance/ACH analysis
)

site_xy      <- niche_res$site_xy
species_objs <- niche_res$species_objs
nic          <- niche_res$niche

# OMI niche parameters
if (!is.null(nic)) {
  niche_param_tab <- niche_res$niche_param
  write.csv(niche_param_tab,
            file.path(OUT_DIR, "OMI_niche_parameters.csv"), row.names = FALSE)
}

## ── Add: occurrence-based OMI niche parameters (sensitivity analysis) ───────
## Re-runs ade4::niche() on a binary 0/1 occurrence matrix instead of raw
## abundance, so marginality/tolerance no longer depend on abundance
## weighting. CH/MVE/DMF centroids and all distance/ACH results above are
## completely unaffected -- this only produces a companion OMI parameter
## table + comparison, addressing the abundance-weighting non-independence
## concern for the tolerance/niche-breadth trait predictor.
niche_res_presence <- make_envelope(
  dudi       = dudi,
  abun       = abun,
  min_occ    = min_occ,
  n_axes     = n_axes,
  pheno_df   = pheno_df,
  site_dates = input$Date,
  count_day  = input$CountDay,
  omi_input  = "presence"
)

if (!is.null(niche_res_presence$niche_param)) {
  niche_param_tab_presence <- niche_res_presence$niche_param
  write.csv(niche_param_tab_presence,
            file.path(OUT_DIR, "OMI_niche_parameters_presence.csv"), row.names = FALSE)
  
  omi_comparison <- niche_param_tab %>%
    dplyr::select(Species, OMI_abund = OMI, Tol_abund = Tol) %>%
    dplyr::inner_join(
      niche_param_tab_presence %>% dplyr::select(Species, OMI_presence = OMI, Tol_presence = Tol),
      by = "Species"
    ) %>%
    dplyr::mutate(
      OMI_spearman_rho = suppressWarnings(cor(OMI_abund, OMI_presence, method = "spearman")),
      Tol_spearman_rho  = suppressWarnings(cor(Tol_abund,  Tol_presence,  method = "spearman"))
    )
  write.csv(omi_comparison,
            file.path(OUT_DIR, "OMI_niche_parameters_comparison.csv"), row.names = FALSE)
  message("\nOMI marginality/tolerance, abundance- vs presence-weighted Spearman rho: ",
          "OMI = ", round(unique(omi_comparison$OMI_spearman_rho), 3),
          ", Tolerance = ", round(unique(omi_comparison$Tol_spearman_rho), 3))
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
    ACH   = classify_ach(rho_adj, pval, alpha = alpha)
  ) %>%
  ungroup() %>%
  left_join(settings_meta, by = "setting")

ach_res_wide <- ach_res %>% 
  dplyr::select(species, setting, rho_adj) %>% 
  pivot_wider(, names_from = "setting", values_from = "rho_adj")

write.csv(ach_res, file.path(OUT_DIR, "ACH_sp_results.csv"), row.names = FALSE)
write.csv(ach_res_wide, file.path(OUT_DIR, "ACH_sp_results_wide.csv"), row.names = FALSE)

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

## ── 10a) Model-ready data ──────────────────────────────────────────────────────
model_df <- build_model_df(species_objs, dist_all, settings_meta,
                            min_occ = min_occ)
write.csv(model_df, file.path(OUT_DIR, "ACH_model_data_long.csv"), row.names = FALSE)

## ── 10b) ICC: inter-setting agreement ────────────────────────────────────────
if (!requireNamespace("irr", quietly = TRUE)) install.packages("irr")
library(irr)

# 28 species × 6 settings (Margin settings already excluded from settings_meta)
rho_wide_icc <- ach_res %>%
  dplyr::select(species, setting, rho_adj) %>%
  tidyr::pivot_wider(names_from = setting, values_from = rho_adj) %>%
  tibble::column_to_rownames("species")

# ICC(2,1): two-way random, absolute agreement, single measures
icc_res <- irr::icc(rho_wide_icc,
                    model   = "twoway",
                    type    = "consistency",
                    unit    = "single",
                    conf.level = 0.95)

icc_summary <- tibble::tibble(
  ICC   = icc_res$value,
  CI_lo = icc_res$lbound,
  CI_hi = icc_res$ubound,
  F_val = icc_res$Fvalue,
  df1   = icc_res$df1,
  df2   = icc_res$df2,
  p     = icc_res$p.value
)

print(icc_summary)
write.csv(icc_summary, file.path(OUT_DIR, "ICC_inter_setting.csv"), row.names = FALSE)

## ── 10c) Centre comparison: Friedman + post-hoc Wilcoxon ─────────────────────
# representative rho (Euclidean + Mahalanobis average) in each centroid
centre_rho <- ach_res %>%
  dplyr::filter(dist_type != "Margin") %>%
  dplyr::group_by(species, centre) %>%
  dplyr::summarise(rho_mean = mean(rho_adj, na.rm = TRUE), .groups = "drop")

# Friedman test
friedman_res <- stats::friedman.test(rho_mean ~ centre | species, data = centre_rho)

friedman_summary <- tibble::tibble(
  statistic = friedman_res$statistic,
  df        = friedman_res$parameter,
  p         = friedman_res$p.value
)

# Post-hoc pairwise Wilcoxon (BH correction, paired by species)
posthoc_res <- pairwise.wilcox.test(
  centre_rho$rho_mean,
  centre_rho$centre,
  paired          = TRUE,
  p.adjust.method = "BH"
) 

posthoc_summary <- as.data.frame(posthoc_res$p.value) %>%
  tibble::rownames_to_column("centre_row") %>%
  tidyr::pivot_longer(-centre_row,
                      names_to  = "centre_col",
                      values_to = "p_BH") %>%
  dplyr::filter(!is.na(p_BH)) %>%
  dplyr::rename(centre_A = centre_row, centre_B = centre_col)

print(friedman_summary)
print(posthoc_summary)

# Save
write.csv(friedman_summary, file.path(OUT_DIR, "centre_friedman.csv"),    row.names = FALSE)
write.csv(posthoc_summary,  file.path(OUT_DIR, "centre_posthoc.csv"),     row.names = FALSE)

## ── Add: Effective sample size correction (Pyper & Peterman, 1998) ────────────
# precondition: species in model_df ordered in time sequence (build_model_df.R odered by Date/CountDay)
# species x setting Spearman's rho (p value) -> resampled subset recorrects
effective_n <- function(x, y, max_lag = NULL) {
  n <- length(x)
  if (is.null(max_lag)) max_lag <- floor(n / 5)     # 관행적 컷오프
  rx <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE)$acf)[-1]
  ry <- as.numeric(acf(y, lag.max = max_lag, plot = FALSE)$acf)[-1]
  tau <- 1 + 2 * sum((1 - (seq_len(max_lag) / n)) * rx * ry)
  n_eff <- n / tau
  pmax(pmin(n_eff, n), 3)                            # 3 <= N_eff <= N
}

eff_n_tab <- model_df %>%
  dplyr::filter(!setting %in% c("E3","E6","E9")) %>%
  dplyr::group_by(setting, species) %>%
  dplyr::filter(dplyr::n() >= min_occ) %>%
  dplyr::group_modify(~ tibble::tibble(
    n_eff = effective_n(rank(.x$dist_signed), rank(.x$log_abund))
  )) %>%
  dplyr::ungroup()

ach_res <- ach_res %>%
  dplyr::filter(dist_type != "Margin") %>%
  dplyr::left_join(eff_n_tab, by = c("setting", "species")) %>%
  dplyr::mutate(
    t_eff = rho * sqrt(pmax(n_eff - 2, 1) / pmax(1 - rho^2, 1e-6)),
    p_eff = 2 * pt(-abs(t_eff), df = pmax(n_eff - 2, 1))
  ) %>%
  dplyr::group_by(setting) %>%
  dplyr::mutate(
    p_eff_adj = bh_adjust(p_eff),
    ACH_effN  = classify_ach(rho_adj, p_eff, alpha = alpha)
  ) %>%
  dplyr::ungroup()

write.csv(ach_res, file.path(OUT_DIR, "ACH_sp_results_effN.csv"), row.names = FALSE)

## ── Add: Moving block bootstrap for Spearman rho ───────────────────────

block_bootstrap_rho <- function(x, y, block_len = 10, n_boot = 9999) {
  n <- length(x)
  obs_rho  <- suppressWarnings(cor(x, y, method = "spearman"))
  n_blocks <- ceiling(n / block_len)
  
  boot_rho <- vapply(seq_len(n_boot), function(b) {
    starts <- sample(seq_len(max(n - block_len + 1, 1)), n_blocks, replace = TRUE)
    idx <- unlist(lapply(starts, function(s) s:min(s + block_len - 1, n)))
    idx <- idx[seq_len(min(length(idx), n))]
    if (length(idx) < 5) return(NA_real_)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }, numeric(1))
  
  boot_rho <- boot_rho[!is.na(boot_rho)]
  
  p_boot <- 2 * min(mean(boot_rho <= 0), mean(boot_rho >= 0))
  p_boot <- min(p_boot, 1)
  
  tibble::tibble(
    rho_obs = obs_rho,
    ci_low  = quantile(boot_rho, 0.025, na.rm = TRUE),
    ci_high = quantile(boot_rho, 0.975, na.rm = TRUE),
    ci_excludes_zero = (quantile(boot_rho, 0.025, na.rm = TRUE) > 0) ||
      (quantile(boot_rho, 0.975, na.rm = TRUE) < 0),
    p_boot = p_boot
  )
}

model_df_neff <- model_df %>%
  dplyr::filter(! setting %in% c("E3","E6","E9")) %>%
  dplyr::left_join(eff_n_tab, by = c("setting","species"))

block_boot_res <- model_df_neff %>%
  dplyr::group_by(setting, species) %>%
  dplyr::filter(dplyr::n() >= min_occ) %>%
  dplyr::group_modify(~ {
    n_grp   <- nrow(.x)
    n_eff_g <- unique(.x$n_eff)
    tau     <- n_grp / n_eff_g
    bl      <- min(max(round(tau), 3), floor(n_grp / 4))
    block_bootstrap_rho(.x$dist_signed, .x$log_abund,
                        block_len = bl, n_boot = 9999) %>%
      dplyr::mutate(block_len_used = bl, .before = 1)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    ACH_boot = dplyr::case_when(
      is.na(rho_obs)                    ~ "Insufficient data",
      ci_excludes_zero & rho_obs < 0    ~ "Supported",
      ci_excludes_zero & rho_obs > 0    ~ "Opposite",
      TRUE                              ~ "Not significant"
    )
  )

write.csv(block_boot_res, file.path(OUT_DIR, "ACH_block_bootstrap.csv"), row.names = FALSE)

## ── Add: merge Method-1 (block bootstrap) classification into ach_res ────────
## (Method-3, AR(1) Bayesian classification, is appended later once
##  bayes_fits_ar1 has been fitted — see "Add-on B" below — and the final
##  ACH_sp_results_corrected.csv / ACH_correction_comparison.csv are written
##  at the very end of the script.)
ach_res <- ach_res %>%
  dplyr::left_join(
    block_boot_res %>% dplyr::select(setting, species, block_len_used,
                                      boot_rho = rho_obs,
                                      boot_ci_low = ci_low, boot_ci_high = ci_high,
                                      ACH_boot = ACH_boot),
    by = c("setting", "species")
  )

## ── 11) Visualisation ─────────────────────────────────────────────────────────
dist_lab_order <- c("Euclidean", "Mahalanobis", "Margin")  # Margin removed

# ── Fig 1: Species niches ────────────────────────────────────────────────────
FIG1_DIR <- file.path(OUT_DIR, "Fig1_sp_niches")
dir.create(FIG1_DIR, showWarnings = FALSE)

fig1_plots <- purrr::map(
  names(species_objs),
  ~ plot_sp_niche(.x, species_objs[[.x]], site_xy)
)

p_fig1 <- patchwork::wrap_plots(fig1_plots,
                                 ncol   = 5,
                                 guides = "collect")
ggplot2::ggsave(file.path(OUT_DIR, "Fig1_species_niches_all_panels.pdf"), p_fig1, width  = 14, height = 14, dpi = 600)
ggplot2::ggsave(file.path(OUT_DIR, "Fig1_species_niches_all_panels.tiff"), p_fig1, width  = 14, height = 14, dpi = 600)

# ── Fig 2-1: Violin plot of Spearman rho ────────────────────────────────────
p_2_1 <- ach_res %>%
  dplyr::filter(setting %in% settings_meta$setting) %>%  # all 6 settings (Margin already excluded)
  dplyr::mutate(
    setting = factor(setting, levels = settings_meta$setting),
    centre  = factor(centre,  levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot(aes(x = setting, y = rho_adj)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  #geom_violin(aes(fill = centre),
  #            colour = "grey85", linewidth = 0.4, trim = FALSE) +
  geom_boxplot(aes(fill = centre), width = 0.5, outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.7, size = 1.4
  ) +
  ylim(c(-0.7,0.7)) +
  scale_fill_manual(values = centre_cols) +
  labs(
    x = NULL,
    y = "Spearman rho",    # If including "margin" setting, "(margin sign reversed)" include in x axis title
    fill = "Centroid",
    title = NULL  
  ) +
  coord_fixed(ratio = 3) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig2_1_Spearman_violin.pdf"),
       p_2_1, width = 8, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_1_Spearman_violin.tiff"),
       p_2_1, width = 8, height = 5, dpi = 300)

# ── Fig 2-2: ACH result heatmap ──────────────────────────────────────────────
p_2_2 <- ach_res %>%
  dplyr::mutate(
    species = factor(species, levels = rev(sort(unique(species)))),
    setting = factor(setting, levels = settings_meta$setting),
    rho_lab = ifelse(!is.na(pval) & pval < alpha,
                     sprintf("%.2f", rho_adj), "")
  ) %>%
  dplyr::filter(is_margin != TRUE) %>%  # all 6 settings (Margin already excluded)
  ggplot(aes(x = setting, y = species, fill = ACH)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = rho_lab), size = 2.8) +
  scale_fill_manual(values = ACH_cols, drop = FALSE) +
  labs(
    x = NULL, y = NULL, fill = NULL,
    title = NULL
    ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig2_2_ACH_heatmap.pdf"),
       p_2_2, width = 7, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_2_ACH_heatmap.tiff"),
       p_2_2, width = 7, height = 7, dpi = 300)

# ── Fig 2-3: species-level ACH result barplot ──────────────────────────────────────
ord_species <- ach_res %>%
  group_by(species) %>%
  summarise(mean_rho = mean(rho_adj, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_rho) %>%
  pull(species)

p_2_3 <- ach_res %>%
  dplyr::mutate(
    species = factor(species, levels = ord_species)
    ) %>%
  dplyr::filter(is_margin != TRUE) %>%  # all 6 settings (Margin already excluded)
  ggplot(aes(x = rho_adj, y = species, fill = ACH)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_col() +
  facet_wrap(~ setting, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = ACH_cols) +
  labs(x = "Spearman rho", y = NULL, fill = NULL,   # If including "margin" setting, "(margin sign reversed)" include in x axis title
       title = NULL
       ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig2_3_species-level_ACH_result.pdf"),
       p_2_3, width = 8, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_3_species-level_ACH_result.tiff"),
       p_2_3, width = 8, height = 10, dpi = 300)

# ── Fig 3: ACH support proportion by setting ─────────────────────────────────
prop_df <- ach_res %>%
  dplyr::count(setting, centre, ACH) %>%
  dplyr::filter(!setting %in% c("E3","E6","E9")) %>%  # all 6 settings (Margin already excluded)
  dplyr::group_by(setting) %>%
  dplyr::mutate(prop = n / sum(n) * 100) %>%
  dplyr::ungroup()

p_3 <- ggplot(
  prop_df,
  aes(x = factor(setting, levels = settings_meta$setting),
      y = prop, fill = ACH)
  ) +
  geom_col(width = 0.65, colour = "white") +
  facet_wrap(~ factor(centre, levels = c("CH", "MVE", "DMF")), scales = "free_x", nrow = 1) +
  scale_fill_manual(values = ACH_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100),
                     expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
  labs(
    x = NULL,
    y = "Proportion of species (%)",
    fill = NULL,
    title = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "Fig3_ACH_support_proportion.pdf"),
       p_3, width = 8, height = 3, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig3_ACH_support_proportion.tiff"),
       p_3, width = 8, height = 3, dpi = 300)

# ── Fig 5: Centre comparison (new) ───────────────────────────────────────────
FIG5_DIR <- file.path(OUT_DIR, "Fig5_centre_comparison")
dir.create(FIG5_DIR, showWarnings = FALSE)

# 5-1: Median rho by centre and distance metric
p_5_1 <- centre_comparison %>%
  dplyr::mutate(
    centre    = factor(centre,    levels = c("CH", "MVE", "DMF")),
    dist_type = factor(dist_type, levels = dist_lab_order)
  ) %>%
  ggplot(aes(x = dist_type, y = median_rho_adj,
                      fill = centre, colour = centre)) +
  geom_col(position = position_dodge(0.7),
           width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_text(aes(label = wilcox_label,
                y = median_rho_adj + ifelse(median_rho_adj >= 0, 0.02, -0.04)),
            position = position_dodge(0.7),
            size = 3.5, fontface = "bold", colour = "black"
            ) +
  scale_fill_manual(values   = centre_cols) +
  scale_colour_manual(values = centre_cols) +
  labs(
    x = "Distance metric", y = "Median Spearman rho",
    fill = "Centre", colour = "Centre",
    title = NULL
  ) +
  coord_fixed(ratio = 7) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(FIG5_DIR, "Fig5_1_centre_median_rho.pdf"),
                p_5_1, width = 7, height = 5, dpi = 300)

# 5-2: Pairwise rho comparison across centres (Euclidean only for clarity)
rho_wide <- ach_res %>%
  dplyr::filter(dist_type != "Margin") %>%
  dplyr::select(species, centre, dist_type, rho_adj) %>%
  tidyr::pivot_wider(names_from = centre, values_from = rho_adj)

if (all(c("CH","MVE","DMF") %in% names(rho_wide))) {
  p_5_2_Eucl <- rho_wide %>%
    dplyr::filter(dist_type == "Euclidean") %>%
    ggplot(aes(x = CH, y = DMF)) +
    geom_point(ggplot2::aes(colour = MVE), size = 2.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 2, colour = "grey40") +
    scale_colour_gradient2(
      low = "#d73027", mid = "grey80", high = "#1a9850", midpoint = 0,
      name = "MVE rho"
    ) +
    coord_fixed(ratio = 0.65) +
    labs(
      x = "CH centre rho",
      y = "DMF centre rho",
      title = NULL  
    ) +
    scale_x_continuous(limits = c(-0.4, 0.45), breaks = seq(-0.4, 0.45, 0.2)) +
    scale_y_continuous(limits = c(-0.6, 0.65), breaks = seq(-0.6, 0.6, 0.3)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p_5_2_Maha <- rho_wide %>%
    dplyr::filter(dist_type == "Mahalanobis") %>%
    ggplot(aes(x = CH, y = DMF)) +
    geom_point(ggplot2::aes(colour = MVE), size = 2.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 2, colour = "grey40") +
    scale_colour_gradient2(
      low = "#d73027", mid = "grey80", high = "#1a9850", midpoint = 0,
      name = "MVE rho"
    ) +
    coord_fixed(ratio = 0.8) +
    labs(
      x = "CH centre rho",
      y = "DMF centre rho",
      title = NULL  
    ) +
    scale_x_continuous(limits = c(-0.4, 0.45), breaks = seq(-0.4, 0.45, 0.2)) +
    scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p_5_2 <- (p_5_2_Eucl | p_5_2_Maha) + plot_layout(widths = c(1,1))

  ggsave(file.path(FIG5_DIR, "Fig5_2_centre_rho_scatter.pdf"),
         p_5_2, width = 9, height = 4, dpi = 300)
  ggsave(file.path(FIG5_DIR, "Fig5_2_centre_rho_scatter.tiff"),
         p_5_2, width = 9, height = 4, dpi = 300)
}

write.csv(rho_wide, file.path(FIG5_DIR, "Fig5_rho_wide.csv"), row.names = FALSE)

comp_cent_panel <- (p_2_1 / p_5_2) + plot_layout(heights = c(1.5, 1),)
ggsave(file.path(FIG5_DIR, "Fig5_comp_centre.pdf"),comp_cent_panel, width = 7, height = 7, dpi = 300)
ggsave(file.path(FIG5_DIR, "Fig5_comp_centre.tiff"),comp_cent_panel, width = 7, height = 7, dpi = 300)

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
  ggplot(aes(x = dist_std, y = y_med, colour = centre, fill = centre)) +
  geom_ribbon(aes(ymin = y_low, ymax = y_high),
              alpha = 0.18, linewidth = 0, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ dist_type, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = centre_cols) +
  scale_fill_manual(values   = centre_cols) +
  labs(
    x = "Standardised distance (0–1)",
    y = "log(predicted abundance)",
    colour = "Centre", fill = "Centre",
    title = "Figure 4-1. GAM distance–abundance curves by centre"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(FIG4_DIR, "Fig4_1_GAM.pdf"),
       p_4_1, width = 11, height = 4.6, dpi = 300)

# 4-2: Species-specific quadratic LM curves
species_curve_df <- model_df %>%
  dplyr::group_by(setting, species) %>%
  dplyr::group_modify(~ fit_quad_lm_sp(.x)) %>%   # FIX B4: correct function name
  dplyr::ungroup()

p_4_2 <- species_curve_df %>%
  mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot2::ggplot(aes(x = dist_std, y = pred,
                               group = species, colour = adj_r2)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  facet_grid(centre ~ dist_type, scales = "free_y") +
  scale_colour_gradient2(
    low = "red", 
    mid = "grey70", 
    high = "green"
    , midpoint = 0.2,
    name = expression(Adj.~R^2)
  ) +
  labs(
    x = "Standardised distance (0–1)",
    y = "log(predicted abundance)",
    title = "Figure 4-2. Species-specific quadratic fitted curves"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(FIG4_DIR, "Fig4_2_species_curves.pdf"),
       p_4_2, width = 9, height = 8, dpi = 300)

# 4-3: Bayesian hierarchical model (optional — requires brms)
if (bayes_ok) {
  fit_brms_quad <- function(df1) {
    brms::brm(
      brms::bf(log_abund ~ 1 + dist_signed + I(dist_signed^2) +
                 (1 + dist_signed + I(dist_signed^2) | species)),
      data    = df1,
      family  = stats::gaussian(),
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
    dr <- brms::posterior_epred(fit, newdata = nd, re_formula = NA)
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
    ggplot(aes(x = dist_std, y = y_med,
                                  colour = centre, fill = centre)) +
    geom_ribbon(aes(ymin = y_low, ymax = y_high),
                alpha = 0.18, linewidth = 0, colour = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ dist_type, ncol = 3, scales = "free_y") +
    scale_colour_manual(values = centre_cols) +
    scale_fill_manual(values   = centre_cols) +
    labs(
      x      = "Standardised distance (0–1)",
      y      = "log(predicted abundance)",
      colour = "Centre", fill = "Centre",
      title  = "Figure 4-3. Bayesian hierarchical quadratic model"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold")
    )
  ggsave(file.path(FIG4_DIR, "Fig4_3_Bayesian.pdf"),
         p_4_3, width = 11, height = 4.6, dpi = 300)
}

## ── Add: AR(1) residual correlation per species ────────────────────────
## Precondition: present the column of CountDay in model_df
fit_brms_quad_ar1 <- function(df1) {
  df1 <- df1 %>% dplyr::arrange(species, CountDay)
  brms::brm(
    brms::bf(log_abund ~ 1 + dist_signed + I(dist_signed^2) +
               (1 + dist_signed + I(dist_signed^2) | species) + 
               ar(time = CountDay, gr = species, p = 1)),
    data    = df1,
    family  = stats::gaussian(),
    chains  = 3, iter = 4000, warmup = 2000, cores = 3, seed = 123,
    refresh = 0,
    control = list(adapt_delta = 0.95),
    prior   = c(
      brms::prior(normal(0, 5),          class = "b"),
      brms::prior(normal(0, 5),          class = "Intercept"),
      brms::prior(student_t(3, 0, 2.5),  class = "sd"),
      brms::prior(lkj(2),                class = "cor"),
      brms::prior(student_t(3, 0, 2.5),  class = "sigma")
    )
  )
}

bayes_fits_ar1 <- purrr::map(
  stats::setNames(settings_meta$setting, settings_meta$setting),
  function(st) {
    message("Bayesian model (AR1): ", st)
    fit_brms_quad_ar1(dplyr::filter(model_df, setting == st))
  }
)

# 4-3: prediction curve
bayes_pred_df_ar1 <- purrr::imap_dfr(bayes_fits_ar1, function(fit, st) {
  dfi   <- dplyr::filter(model_df, setting == st)
  is_mg <- unique(dfi$is_margin)
  
  nd <- tibble::tibble(
    species     = unique(dfi$species)[1],
    dist_std    = seq(0, 1, length.out = 100),
    dist_signed = if (is_mg) -seq(0, 1, length.out = 100)
    else        seq(0, 1, length.out = 100),
    CountDay    = median(dfi$CountDay)   # AR항 제외 시 값은 무의미하나 컬럼은 필요
  )
  
  dr <- brms::posterior_epred(
    fit, newdata = nd,
    re_formula   = NA,       # 커뮤니티 수준: 종별 랜덤효과 제외
    incl_autocor = FALSE     # AR(1) 잔차 구조 예측에서 제외
  )
  q <- apply(dr, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  tibble::tibble(
    setting   = st,
    centre    = unique(dfi$centre),
    dist_type = unique(dfi$dist_type),
    dist_std  = nd$dist_std,
    y_low     = q[1, ], y_med = q[2, ], y_high = q[3, ]
  )
})

p_4_3_ar1 <- bayes_pred_df_ar1 %>%
  dplyr::mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH","MVE","DMF"))
  ) %>%
  ggplot(aes(x = dist_std, y = y_med,
             colour = centre, fill = centre)) +
  geom_ribbon(aes(ymin = y_low, ymax = y_high),
              alpha = 0.18, linewidth = 0, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ dist_type, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = centre_cols) +
  scale_fill_manual(values   = centre_cols) +
  labs(
    x      = "Standardised distance (0–1)",
    y      = "log(predicted abundance)",
    colour = "Centre", fill = "Centre",
    title  = "Figure 4-3. Bayesian hierarchical quadratic model with AR(1) residuals"
  ) +
  ylim(c(1.5,4)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(FIG4_DIR, "Fig4_3_Bayesian_AR1_corrected.pdf"),
       p_4_3_ar1, width = 11, height = 4.6, dpi = 300)

## LOO 비교: AR(1) 추가가 모델 적합도를 개선하는지 확인
loo_compare_tab <- purrr::map_dfr(names(bayes_fits), function(st) {
  loo0 <- loo::loo(bayes_fits[[st]])
  loo1 <- loo::loo(bayes_fits_ar1[[st]])
  cmp  <- loo::loo_compare(loo0, loo1)
  tibble::tibble(setting = st,
                 elpd_diff = cmp[2, "elpd_diff"],
                 se_diff   = cmp[2, "se_diff"])
})
write.csv(loo_compare_tab, file.path(OUT_DIR, "LOO_AR1_comparison.csv"), row.names = FALSE)

## ── Add: Method 3 classification — AR(1) Bayesian P(beta1 < 0) ──────────────
## For each setting, extract the population slope (b_dist_signed) plus the
## species-level random-slope deviation from the AR(1) model, and classify
## each species from the posterior probability that its linear slope is
## negative (ACH support), mirroring the one used for Table S6 in Fig2_final.R.
extract_ar1_beta1 <- function(fit, sp_list) {
  draws_all <- posterior::as_draws_df(fit)
  b1_pop    <- draws_all[["b_dist_signed"]]

  purrr::map_dfr(sp_list, function(sp) {
    re_col <- paste0("r_species[", sp, ",dist_signed]")
    if (!re_col %in% names(draws_all)) return(NULL)
    b1_sp <- b1_pop + draws_all[[re_col]]
    tibble::tibble(species = sp, P_beta1_neg_AR1 = mean(b1_sp < 0))
  })
}

ar1_beta_tab <- purrr::imap_dfr(bayes_fits_ar1, function(fit, st) {
  sp_list <- unique(dplyr::filter(model_df, setting == st)$species)
  tryCatch(
    extract_ar1_beta1(fit, sp_list) %>% dplyr::mutate(setting = st),
    error = function(e) {
      message("AR1 beta extraction failed for ", st, ": ", e$message)
      NULL
    }
  )
})

ar1_beta_tab <- ar1_beta_tab %>%
  dplyr::mutate(
    ACH_AR1 = dplyr::case_when(
      P_beta1_neg_AR1 >= 0.95 ~ "Supported",
      P_beta1_neg_AR1 <= 0.05 ~ "Opposite",
      TRUE                    ~ "Not significant"
    )
  )

write.csv(ar1_beta_tab, file.path(OUT_DIR, "ACH_AR1_beta_classification.csv"), row.names = FALSE)

## ── Add: consolidated corrected results + 3-method agreement summary ────────
ach_res <- ach_res %>%
  dplyr::left_join(
    ar1_beta_tab %>% dplyr::select(setting, species, P_beta1_neg_AR1, ACH_AR1),
    by = c("setting", "species")
  )

write.csv(ach_res, file.path(OUT_DIR, "ACH_sp_results_corrected.csv"), row.names = FALSE)

ACH_correction_comparison <- ach_res %>%
  dplyr::mutate(
    agree_boot_effN = ACH_boot == ACH_effN,
    agree_boot_AR1  = ACH_boot == ACH_AR1,
    agree_effN_AR1  = ACH_effN == ACH_AR1,
    agree_all_3     = agree_boot_effN & agree_boot_AR1 & agree_effN_AR1
  ) %>%
  dplyr::group_by(setting) %>%
  dplyr::summarise(
    n_species                = dplyr::n(),
    pct_supported_nominal    = 100 * mean(ACH      == "Supported", na.rm = TRUE),
    pct_supported_boot       = 100 * mean(ACH_boot == "Supported", na.rm = TRUE),
    pct_supported_effN       = 100 * mean(ACH_effN == "Supported", na.rm = TRUE),
    pct_supported_AR1        = 100 * mean(ACH_AR1  == "Supported", na.rm = TRUE),
    pct_agree_all_3          = 100 * mean(agree_all_3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(settings_meta, by = "setting")

write.csv(ACH_correction_comparison,
          file.path(OUT_DIR, "ACH_correction_comparison.csv"), row.names = FALSE)
print(ACH_correction_comparison)

# (optional) 4-4: Species-specific posterior predicted curves
# Bayesian hierarchical quadratic model
bayes_sp_pred_df <- purrr::imap_dfr(bayes_fits, function(fit, st) {
  dfi    <- dplyr::filter(model_df, setting == st)
  is_mg  <- unique(dfi$is_margin)
  sp_list <- unique(dfi$species)
  
  purrr::map_dfr(sp_list, function(sp) {
    nd <- tibble::tibble(
      species     = sp,
      dist_std    = seq(0, 1, length.out = 50),
      dist_signed = if (is_mg) -seq(0, 1, length.out = 50)
      else        seq(0, 1, length.out = 50)
    )
    
    dr <- tryCatch(
      brms::posterior_epred(fit, newdata = nd, re_formula = NULL),
      error = function(e) NULL
    )
    if (is.null(dr)) return(NULL)
    
    q <- apply(dr, 2, quantile,
               probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    
    tibble::tibble(
      setting   = st,
      centre    = unique(dfi$centre),
      dist_type = unique(dfi$dist_type),
      species   = sp,
      dist_std  = nd$dist_std,
      y_low     = q[1, ],
      y_med     = q[2, ],
      y_high    = q[3, ]
    )
  })
})

# Adj. R² 계산 (OLS quadratic, 색상 매핑용)
adjr2_df <- model_df %>%
  dplyr::group_by(setting, species) %>%
  dplyr::summarise(
    adj_r2 = tryCatch({
      m <- lm(log_abund ~ dist_signed + I(dist_signed^2), data = dplyr::cur_data())
      summary(m)$adj.r.squared
    }, error = function(e) NA_real_),
    .groups = "drop"
  )

bayes_sp_pred_df <- bayes_sp_pred_df %>%
  dplyr::left_join(adjr2_df, by = c("setting", "species"))

p_4_4 <- bayes_sp_pred_df %>%
  dplyr::mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH", "MVE", "DMF"))
  ) %>%
  ggplot(aes(x = dist_std, y = y_med,
             group = species, colour = adj_r2)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  facet_grid(centre ~ dist_type, scales = "free_y") +
  scale_colour_gradient2(
    low = "red", 
    mid = "grey70", 
    high = "green",
    midpoint = 0.2,
    name     = expression(Adj.~R^2)
  ) +
  labs(
    x     = "Standardised distance (0–1)",
    y     = "log(predicted abundance)",
    title = "Figure 4-4. Species-specific Bayesian posterior predicted curves"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(FIG4_DIR, "Fig4_4_Bayesian_sp_curves.pdf"),
                p_4_4, width = 9, height = 8, dpi = 300)

## ── Add: AR(1)-corrected counterpart of Fig 4-4 (species posterior curves) ──
bayes_sp_pred_df_ar1 <- purrr::imap_dfr(bayes_fits_ar1, function(fit, st) {
  dfi     <- dplyr::filter(model_df, setting == st)
  is_mg   <- unique(dfi$is_margin)
  sp_list <- unique(dfi$species)

  purrr::map_dfr(sp_list, function(sp) {
    nd <- tibble::tibble(
      species     = sp,
      dist_std    = seq(0, 1, length.out = 50),
      dist_signed = if (is_mg) -seq(0, 1, length.out = 50)
      else        seq(0, 1, length.out = 50)
    )
    dr <- tryCatch(
      brms::posterior_epred(fit, newdata = nd, re_formula = NULL, 
                            incl_autocor = FALSE),
      error = function(e) NULL
    )
    if (is.null(dr)) return(NULL)
    q <- apply(dr, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    tibble::tibble(
      setting = st, centre = unique(dfi$centre), dist_type = unique(dfi$dist_type),
      species = sp, dist_std = nd$dist_std,
      y_low = q[1, ], y_med = q[2, ], y_high = q[3, ]
    )
  })
}) %>%
  dplyr::left_join(adjr2_df, by = c("setting", "species"))

p_4_4_ar1 <- bayes_sp_pred_df_ar1 %>%
  dplyr::mutate(
    dist_type = factor(dist_type, levels = dist_lab_order),
    centre    = factor(centre,    levels = c("CH", "MVE", "DMF"))
  ) %>%
  ggplot(aes(x = dist_std, y = y_med, group = species, colour = adj_r2)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  facet_grid(centre ~ dist_type, scales = "free_y") +
  scale_colour_gradient2(
    low = "red", mid = "grey70", high = "green", midpoint = 0.2,
    name = expression(Adj.~R^2)
  ) +
  labs(
    x = "Standardised distance (0–1)", y = "log(predicted abundance)",
    title = "Figure 4-4 (AR(1)-corrected). Species-specific Bayesian posterior predicted curves"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(FIG4_DIR, "Fig4_4_Bayesian_sp_curves_AR1_corrected.pdf"),
       p_4_4_ar1, width = 9, height = 8, dpi = 300)

## ── Add: temporal-AC-corrected versions of Fig2-2 / Fig2-3 / Fig3 ───────────
## Uses ACH_effN (effective-N / Pyper & Peterman corrected significance) as
## the primary "corrected" classification for the main manuscript figures,
## since it shares the same rho point-estimate as the nominal analysis and
## only adjusts the p-value for temporal autocorrelation — making it the
## most directly comparable corrected counterpart to the original figures.
## ACH_boot and ACH_AR1 are reported as sensitivity checks in
## ACH_correction_comparison.csv / ACH_sp_results_corrected.csv.
## Original (nominal) Fig2_2 / Fig2_3 / Fig3 files above are left untouched.

p_2_2_corrected <- ach_res %>%
  dplyr::mutate(
    species = factor(species, levels = rev(sort(unique(species)))),
    setting = factor(setting, levels = settings_meta$setting),
    rho_lab = ifelse(!is.na(p_eff) & p_eff < alpha,
                     sprintf("%.2f", rho_adj), "")
  ) %>%
  dplyr::filter(setting %in% settings_meta$setting) %>%
  ggplot(aes(x = setting, y = species, fill = ACH_effN)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = rho_lab), size = 2.8) +
  scale_fill_manual(values = ACH_cols, drop = FALSE, name = "ACH (effective-N corrected)") +
  labs(x = NULL, y = NULL, fill = NULL, title = NULL) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig2_2_ACH_heatmap_corrected.pdf"),
       p_2_2_corrected, width = 7, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_2_ACH_heatmap_corrected.tiff"),
       p_2_2_corrected, width = 7, height = 7, dpi = 300)

p_2_3_corrected <- ach_res %>%
  dplyr::mutate(species = factor(species, levels = ord_species)) %>%
  dplyr::filter(setting %in% settings_meta$setting) %>%
  ggplot(aes(x = rho_adj, y = species, fill = ACH_effN)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_col() +
  facet_wrap(~ setting, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = ACH_cols, name = "ACH (effective-N corrected)") +
  labs(x = "Spearman rho", y = NULL, fill = NULL, title = NULL) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig2_3_species-level_ACH_result_corrected.pdf"),
       p_2_3_corrected, width = 8, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_3_species-level_ACH_result_corrected.tiff"),
       p_2_3_corrected, width = 8, height = 10, dpi = 300)

prop_df_corrected <- ach_res %>%
  dplyr::count(setting, centre, ACH_effN) %>%
  dplyr::filter(setting %in% settings_meta$setting) %>%
  dplyr::group_by(setting) %>%
  dplyr::mutate(prop = n / sum(n) * 100) %>%
  dplyr::ungroup()

p_3_corrected <- ggplot(
  prop_df_corrected,
  aes(x = factor(setting, levels = settings_meta$setting),
      y = prop, fill = ACH_effN)
  ) +
  geom_col(width = 0.65, colour = "white") +
  facet_wrap(~ factor(centre, levels = c("CH", "MVE", "DMF")), scales = "free_x", nrow = 1) +
  scale_fill_manual(values = ACH_cols, drop = FALSE, name = "ACH (effective-N corrected)") +
  scale_y_continuous(limits = c(0, 100),
                     expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
  labs(x = NULL, y = "Proportion of species (%)", fill = NULL, title = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "Fig3_ACH_support_proportion_corrected.pdf"),
       p_3_corrected, width = 8, height = 3, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig3_ACH_support_proportion_corrected.tiff"),
       p_3_corrected, width = 8, height = 3, dpi = 300)

message("\n✓ All outputs written to: ", OUT_DIR)
message("  - Nominal figures   : Fig2_2_ACH_heatmap.pdf / Fig2_3_*.pdf / Fig3_*.pdf")
message("  - Corrected figures : *_corrected.pdf (ACH_effN classification)")
message("  - 3-method comparison: ACH_correction_comparison.csv / ACH_sp_results_corrected.csv")


