# ============================================================
#  Functional Traits vs. ACH Support
#  Fisher's Z Regression Analysis
#
#  Method : Effective-N weighted regression (PRIMARY)
#           + naive (unweighted) OLS kept as a reported comparison
#  Response: Fisher's Z of Spearman's rho (abundance–distance)
#  Input   : ACH_sp_results_corrected.csv / dino_functraits.xlsx
#            / OMI_results.xlsx (occurrence-based, post circularity fix)
#
#  ── Design note (why effective-N weighting is PRIMARY here) ──────────────
#  Each species' z_rho carries a different precision depending on how much
#  temporal autocorrelation was present in its own 411-day series (captured
#  by n_eff, already computed upstream). Treating all species' z_rho as
#  equally precise (plain OLS) ignores this. Weighting by 1/z_se_eff^2 is a
#  standard inverse-variance meta-regression specification — this is the
#  statistically appropriate model, not a "looser" one. Naive OLS is kept
#  only as a secondary, reported cross-check (mirrors how ACH_sp_results
#  already reports both raw rho/p and effective-N-adjusted rho/p_adj).
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ── 0. Packages ───────────────────────────────────────────────
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr", "stringr",
  "ggplot2", "patchwork", "emmeans"
)

if (length(pkgs)) install.packages(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. Paths ──────────────────────────────────────────────────
INPUT_DIR  <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology"
RHO_FILE   <- "ACH/Output_260804/ACH_sp_results_corrected.csv"
FT_FILE    <- "dino_functraits.xlsx"
# NOTE: must point to the OCCURRENCE-based OMI parameter table (post
# circularity fix — ade4::niche() re-run on presence/absence), NOT the
# original abundance-weighted version. Otherwise Marginality/Niche_breadth
# as predictors of an abundance-derived rho would be circular. Confirm the
# file below is the corrected output before running.
NICHE_FILE <- "OMI/Output_260804/OMI_results_260804.xlsx"

WORK_DIR <- file.path(INPUT_DIR, "ACH")
OUT_DIR  <- file.path(WORK_DIR, paste0("Output_FT", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 2. Global settings ────────────────────────────────────────
METHOD      <- "E8"
TROPHIC_REF <- "CM"
CONT_TRAITS <- c("log_Biovolume", "Speed_max", "Marginality", "Niche_breadth")
CAT_TRAITS  <- c("Spicule", "Colony", "Cyst", "Trophic_type")

troph_cols <- c(CM = "#EFC050", pSNCM = "#DD4124", HET = "#0F4C81", OPA = "#378661")
bin_cols   <- c("0"= "#B2BEB5", "1" = "#36454F") 

# ── 3. Load & merge input data ────────────────────────────────
rho_df <- read.csv(file.path(INPUT_DIR, RHO_FILE)) %>%
  dplyr::filter(setting == METHOD) %>%
  dplyr::rename(Species = species)

nic_df <- read.xlsx(file.path(INPUT_DIR, NICHE_FILE), sheet = "OMI_params") %>%
  dplyr::select(Species, Marginality = OMI, Niche_breadth = Tol)

ft_df <- read.xlsx(file.path(INPUT_DIR, FT_FILE), sheet = "FuncTrait32") %>%
  dplyr::select(TaxID, Abbrevration, Biovolume, Spicule, Colony,
                Speed_max, Trophic_type, Cyst) %>%
  dplyr::rename(Species = Abbrevration) %>%
  left_join(nic_df, by = "Species") %>%
  mutate(
    log_Biovolume = log10(Biovolume)
  ) %>%
  dplyr::select(-Biovolume) %>%
  dplyr::filter(
    Species %in% rho_df$Species,
    !(Species == "Nsci" & Trophic_type == "eSNCM"),
    !(Species == "Poly" & Trophic_type == "CM")
  )

# ── 4. Fisher's Z transformation ──────────────────────────────
fisher_z    <- function(r) atanh(pmin(pmax(r, -0.9999), 0.9999))
fisher_z_se <- function(n) 1 / sqrt(n - 3)

rho_df <- rho_df %>%
  mutate(
    z_rho    = fisher_z(rho),
    z_se     = fisher_z_se(n),
    z_se_eff = fisher_z_se(n_eff),
    w_eff    = 1 / (z_se_eff^2)   # inverse-variance weight, computed once here
  )

# final input data
mydf <- rho_df %>%
  dplyr::select(Species, setting, ACH, centre, dist_type, n, z_rho, z_se, z_se_eff, w_eff) %>%
  left_join(ft_df, by = "Species") %>%
  relocate(c(log_Biovolume), .after = "w_eff")

# ── 5. Helper: significance label ─────────────────────────────
sig_label <- function(p) {
  as.character(cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   labels = c("***", "**", "*", "\u2020", "ns")))
}

# ── 6. descriptive statistics of binary and catogorical data ───
lapply(c("Spicule", "Colony", "Cyst", "Trophic_type"), function(x){
  
  mydf %>%
    group_by(.data[[x]]) %>%
    summarise(
      N = n(),
      Median = median(z_rho, na.rm = TRUE),
      IQR = IQR(z_rho, na.rm = TRUE),
      .groups = "drop"
    )
  
})

# ── 7. Trait analysis: naive OLS (secondary) + effective-N WLS (PRIMARY) ──

## 7a. Continuous traits
cont_results <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, w_eff, all_of(trait)) %>% drop_na()
  n   <- nrow(sub)
  sub$x_std <- scale(sub[[trait]])[, 1]
  
  # naive (unweighted) — secondary / reported for comparison only
  ols   <- lm(z_rho ~ x_std, data = sub)
  beta_naive <- coef(ols)["x_std"]
  ci_naive   <- confint(ols)["x_std", ]
  p_naive    <- summary(ols)$coefficients["x_std", "Pr(>|t|)"]
  
  # effective-N weighted (PRIMARY) — inverse-variance meta-regression
  wls   <- lm(z_rho ~ x_std, data = sub, weights = w_eff)
  beta_effN <- coef(wls)["x_std"]
  ci_effN   <- confint(wls)["x_std", ]
  p_effN    <- summary(wls)$coefficients["x_std", "Pr(>|t|)"]
  
  # non-parametric cross-check (unweighted; cross-sectional across species,
  # no temporal autocorrelation at this level, so no effN analogue needed)
  sp <- cor.test(sub[[trait]], sub$z_rho, method = "spearman", exact = FALSE)
  
  data.frame(
    Trait = trait, Type = "continuous", n = n,
    Beta_naive = beta_naive, CI_low_naive = ci_naive[1], CI_high_naive = ci_naive[2], p_naive = p_naive,
    Beta_effN  = beta_effN,  CI_low_effN  = ci_effN[1],  CI_high_effN  = ci_effN[2],  p_effN  = p_effN,
    Spearman_rho = sp$estimate, p_Spearman = sp$p.value,
    row.names = NULL
  )
})

## 7b. Binary categorical: Spicule / Colony / Cyst
bin_results <- lapply(c("Spicule", "Colony", "Cyst"), function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, w_eff, all_of(trait)) %>% drop_na()
  sub[[trait]] <- factor(sub[[trait]])
  n   <- nrow(sub)
  key <- paste0(trait, "1")
  
  # naive — secondary
  ols  <- lm(as.formula(paste("z_rho ~", trait)), data = sub)
  beta_naive <- coef(ols)[key]
  ci_naive   <- confint(ols)[key, ]
  p_naive    <- summary(ols)$coefficients[key, "Pr(>|t|)"]
  
  # effective-N weighted — PRIMARY
  wls  <- lm(as.formula(paste("z_rho ~", trait)), data = sub, weights = w_eff)
  beta_effN <- coef(wls)[key]
  ci_effN   <- confint(wls)[key, ]
  p_effN    <- summary(wls)$coefficients[key, "Pr(>|t|)"]
  
  # Mann-Whitney U test (non-parametric cross-check; unweighted for the
  # same reason as the Spearman check above)
  g    <- split(sub$z_rho, sub[[trait]])
  mw_p <- wilcox.test(g[["0"]], g[["1"]], exact = FALSE)$p.value
  
  data.frame(
    Trait = trait, Type = "binary", n = n,
    Beta_naive = beta_naive, CI_low_naive = ci_naive[1], CI_high_naive = ci_naive[2], p_naive = p_naive,
    Beta_effN  = beta_effN,  CI_low_effN  = ci_effN[1],  CI_high_effN  = ci_effN[2],  p_effN  = p_effN,
    Spearman_rho = NA, mw_p = mw_p,
    row.names = NULL
  )
})

## 7c. Trophic_type — all pairwise comparisons (emmeans)
##     Raw (unadjusted) p-values are pulled here; BH-FDR is applied ONCE,
##     across the full combined set of tests in step 7d — pulling
##     adjust = "BH" here AND again in 7d would double-correct.
sub_t <- mydf %>% dplyr::select(z_rho, w_eff, Trophic_type) %>% drop_na()
sub_t$Trophic_type <- factor(sub_t$Trophic_type)
sub_t$Trophic_type <- relevel(sub_t$Trophic_type, ref = "OPA")

kw_res <- kruskal.test(z_rho ~ Trophic_type, data = sub_t)
kw_p   <- kw_res$p.value

# naive (unweighted) — secondary
ols_t   <- lm(z_rho ~ Trophic_type, data = sub_t)
emm_naive <- emmeans(ols_t, ~ Trophic_type)
pairs_naive <- contrast(emm_naive, method = "pairwise", adjust = "none") %>%
  as.data.frame() %>%
  dplyr::transmute(Trait = contrast, Beta_naive = estimate, p_naive = p.value)
ci_naive_t <- confint(contrast(emm_naive, method = "pairwise", adjust = "none")) %>%
  as.data.frame() %>%
  dplyr::select(CI_low_naive = lower.CL, CI_high_naive = upper.CL)

# effective-N weighted — PRIMARY
wls_t   <- lm(z_rho ~ Trophic_type, data = sub_t, weights = w_eff)
emm_effN <- emmeans(wls_t, ~ Trophic_type)
pairs_effN <- contrast(emm_effN, method = "pairwise", adjust = "none") %>%
  as.data.frame() %>%
  dplyr::transmute(Trait = contrast, Beta_effN = estimate, p_effN = p.value)
ci_effN_t <- confint(contrast(emm_effN, method = "pairwise", adjust = "none")) %>%
  as.data.frame() %>%
  dplyr::select(CI_low_effN = lower.CL, CI_high_effN = upper.CL)

trophic_results <- pairs_naive %>%
  dplyr::bind_cols(ci_naive_t) %>%
  dplyr::left_join(dplyr::bind_cols(pairs_effN, ci_effN_t), by = "Trait") %>%
  dplyr::mutate(
    Trait = paste0("Trophic: ", Trait),
    Type  = "categorical",
    n     = nrow(sub_t),
    Spearman_rho = NA,
    p_Spearman   = kw_p
  ) %>%
  dplyr::select(Trait, Type, n,
                Beta_naive, CI_low_naive, CI_high_naive, p_naive,
                Beta_effN,  CI_low_effN,  CI_high_effN,  p_effN,
                Spearman_rho, p_Spearman)

## 7d. Combine & apply BH-FDR ONCE, to the PRIMARY (effective-N) p-values
##     across the full family of continuous + binary + trophic tests.
##     Naive p-values are retained for reported comparison only and do NOT
##     drive significance calls (mirrors species-level rho_adj/p_adj usage).
res_df <- bind_rows(cont_results, bin_results, trophic_results) %>%
  mutate(
    p_FDR = p.adjust(p_effN, method = "BH"),
    sig   = sig_label(p_FDR)
  )

# ── 8. Trophic-controlled regression (fixed effect) ────────────
#  Replaces lmer — Trophic_type treated as fixed covariate
#  Rationale: only 4 trophic levels → insufficient for random effect
#             variance estimation; Trophic_type is itself a trait of interest
#  Both naive and effective-N weighted versions reported; effN is PRIMARY.
ctrl_results <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, w_eff, all_of(trait), Trophic_type, Colony) %>% drop_na()
  if (nrow(sub) < 15) return(NULL)
  
  sub$x_std        <- scale(sub[[trait]])[, 1]
  sub$Trophic_type <- relevel(factor(sub$Trophic_type), ref = TROPHIC_REF)
  
  ols  <- lm(z_rho ~ x_std + Trophic_type, data = sub)
  beta_naive <- coef(ols)["x_std"]
  ci_naive   <- confint(ols)["x_std", ]
  p_naive    <- summary(ols)$coefficients["x_std", "Pr(>|t|)"]
  
  wls  <- lm(z_rho ~ x_std + Trophic_type, data = sub, weights = w_eff)
  beta_effN <- coef(wls)["x_std"]
  ci_effN   <- confint(wls)["x_std", ]
  p_effN    <- summary(wls)$coefficients["x_std", "Pr(>|t|)"]
  
  data.frame(
    Trait = trait,
    Beta_ctrl_naive = beta_naive, CI_low_ctrl_naive = ci_naive[1], CI_high_ctrl_naive = ci_naive[2], p_ctrl_naive = p_naive,
    Beta_ctrl_effN  = beta_effN,  CI_low_ctrl_effN  = ci_effN[1],  CI_high_ctrl_effN  = ci_effN[2],  p_ctrl_effN  = p_effN,
    row.names = NULL
  )
})

ctrl_df <- bind_rows(ctrl_results) %>%
  mutate(
    p_FDR_ctrl = p.adjust(p_ctrl_effN, method = "BH"),
    sig_ctrl   = sig_label(p_FDR_ctrl)
  )

# ── 9. Forest plot (PRIMARY = effective-N weighted) ────────────
make_forest <- function(data, beta_col, ci_low, ci_high, sig_col,
                        xlim = c(-1, 1)) {
  data <- data %>%
    mutate(color = case_when(
      .data[[ci_low]]  > 0 ~ "pos",
      .data[[ci_high]] < 0 ~ "neg",
      TRUE                 ~ "ns"
    ),
    sig_fill = case_when(
      color == "pos" & .data[[sig_col]] != "ns" ~ "pos",
      color == "neg" & .data[[sig_col]] != "ns" ~ "neg",
      TRUE ~ "ns"
    ))
  
  ggplot(data, aes(x = .data[[beta_col]],
                   y = reorder(Trait, .data[[beta_col]]))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbarh(aes(xmin = .data[[ci_low]], xmax = .data[[ci_high]],
                       color = color), height = 0, linewidth = 0.8) +
    geom_point(aes(color = color, fill = sig_fill), pch = 21, size = 4) +
    scale_color_manual(values = c(pos = "#d62728", neg = "#1f77b4", ns = "#aaaaaa"), guide  = "none") +
    scale_fill_manual(values = c(pos = "#d62728", neg = "#1f77b4", ns = "white"), guide  = "none") +
    coord_cartesian(xlim = xlim) +
    labs(x = "Standardized \u03b2 (95% CI)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "bottom", legend.key = element_rect(color = NA, fill = NA),
          axis.text = element_text(color = "black"))
}

# Primary (effective-N weighted) univariate forest plot
p_uni  <- make_forest(res_df, "Beta_effN", "CI_low_effN", "CI_high_effN", "sig",
                      xlim = c(-1, 0.6))

# Trophic-controlled forest plot (primary = effN)
p_ctrl <- make_forest(ctrl_df %>% dplyr::rename(sig = sig_ctrl),
                      "Beta_ctrl_effN", "CI_low_ctrl_effN", "CI_high_ctrl_effN", "sig",
                      xlim = c(-0.4, 0.4))

forest_combined <- p_uni + p_ctrl +
  plot_annotation(
    title = sprintf("Functional Traits vs. ACH Support (%s, Fisher's Z, effective-N weighted)", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_forest.pdf", METHOD)),
       forest_combined, width = 10, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_forest.tiff", METHOD)),
       forest_combined, width = 10, height = 6, dpi = 300)

# ── 11. Scatter plots (fit line reflects the effective-N weighted model) ──
scatter_plots <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>%
    dplyr::select(Species, z_rho, w_eff, all_of(trait), Trophic_type) %>%
    drop_na()
  r <- res_df %>% filter(Trait == trait)
  
  x_label <- dplyr::case_when(
    trait == "log_ESD"       ~ "log\u2081\u2080 ESD (\u03bcm)",
    trait == "log_Biovolume" ~ "log\u2081\u2080 Biovolume (\u03bcm\u00b3)",
    trait == "Speed_max"       ~ "Max. Speed (\u03bcm s\u207b\u00b9)",
    TRUE                     ~ trait
  )
  
  ggplot(sub, aes(x = .data[[trait]], y = z_rho, color = Trophic_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_smooth(aes(weight = w_eff), method = "lm", se = TRUE,
                color = "black", linewidth = 1, alpha = 0.15) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = troph_cols, name = "Trophic type") +
    labs(x      = x_label,
         y      = sprintf("Fisher's Z (%s \u03c1)", METHOD),
         title  = trait       # If you want to include the statistic results in plot, change to sprintf("%s\n\u03b2=%.3f, p_FDR=%.3f %s", trait, r$Beta_effN, r$p_FDR, as.character(r$sig))
    ) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "right", legend.key = element_rect(color = NA, fill = NA),
          axis.text = element_text(color = "black"))
})

scatter_combined <- wrap_plots(scatter_plots, ncol = 2) +
  plot_annotation(
    title = sprintf("Continuous Traits vs. Fisher's Z (%s \u03c1)", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_scatter.pdf", METHOD)),
       scatter_combined, width = 13, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_scatter.tiff", METHOD)),
       scatter_combined, width = 13, height = 7, dpi = 300)

# ── 12. Categorical box plots ─────────────────────────────────
cat_plots <- lapply(CAT_TRAITS, function(trait) {
  sub <- mydf %>%
    dplyr::select(z_rho, all_of(trait)) %>%
    drop_na() %>%
    mutate(across(all_of(trait), as.character))
  
  fill_cols <- if (trait == "Trophic_type") {troph_cols} else {bin_cols}
  ord <- if (trait == "Trophic_type") {c("OPA", "CM", "pSNCM", "HET")} else {NULL}
  
  ggplot(sub, aes(x = .data[[trait]], y = z_rho, fill = .data[[trait]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
    scale_x_discrete(limits = ord) +
    scale_fill_manual(values = fill_cols, guide = "none") +
    ylim(c(-0.6, 0.6)) +
    labs(x = NULL,
         y = sprintf("Fisher's Z (%s \u03c1)", METHOD),
         title = trait) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "right", legend.key = element_rect(color = NA, fill = NA),
          axis.text = element_text(color = "black"))
})

cat_combined <- wrap_plots(cat_plots, ncol = 2) +
  plot_annotation(
    title = sprintf("Categorical Traits vs. Fisher's Z (%s \u03c1)", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_categorical.pdf", METHOD)),
       cat_combined, width = 8, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_categorical.tiff", METHOD)),
       cat_combined, width = 8, height = 8, dpi = 300)

# Figure 4
scatter_plots_nolegend <- lapply(
  scatter_plots,
  function(p) p + theme(legend.position = "none")
)

combined_plot <- wrap_plots(
  wrap_plots(scatter_plots_nolegend, ncol = 4),
  wrap_plots(cat_plots, ncol = 4),
  ncol = 1
)

ggsave(file.path(OUT_DIR, sprintf("Fig4_FT_ACH_%s.pdf", METHOD)),
       combined_plot, width = 12, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("Fig4_FT_ACH_%s.tiff", METHOD)),
       combined_plot, width = 12, height = 6, dpi = 300)

# ── 13. Save results ──────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Univariate")
addWorksheet(wb, "TrophicControlled")
writeData(wb, "Univariate", res_df)
writeData(wb, "TrophicControlled", ctrl_df)
saveWorkbook(wb, file.path(OUT_DIR, sprintf("FT_ACH_%s_results.xlsx", METHOD)), overwrite = TRUE)

cat("\nDone. All outputs saved to:\n", OUT_DIR, "\n")
cat("\nPrimary results use effective-N weighted regression (Beta_effN / p_FDR).\n",
    "Naive (unweighted) columns are retained for reported comparison only.\n")
