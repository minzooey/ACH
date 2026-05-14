# ============================================================
#  Functional Traits vs. ACH Support
#  Fisher's Z Regression Analysis
#
#  Method : OLS univariate + Trophic_type-controlled OLS
#  Response: Fisher's Z of Spearman's rho (abundance–distance)
#  Input   : ACH_sp_results.csv / dino_functraits_origin.xlsx
#            / OMI_results.xlsx
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
RHO_FILE   <- "ACH/Output_260502/ACH_sp_results.csv"
FT_FILE    <- "dino_functraits_origin.xlsx"
NICHE_FILE <- "OMI/Output_260506/OMI_results_260506.xlsx"

WORK_DIR <- file.path(INPUT_DIR, "ACH")
OUT_DIR  <- file.path(WORK_DIR, paste0("Output_FT", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 2. Global settings ────────────────────────────────────────
METHOD      <- "E2"
TROPHIC_REF <- "CM"
CONT_TRAITS <- c("log_ESD", "log_Biovolume", "Speed_max", "Marginality", "Niche_breadth")
CAT_TRAITS  <- c("Spincule", "Colony", "Toxin", "Trophic_type")

troph_cols <- c(CM = "#EFC050", pSNCM = "#DD4124", HET = "#0F4C81", OPA = "#378661")
bin_cols   <- c("0"= "#B2BEB5", "1" = "#36454F") 

# ── 3. Load & merge input data ────────────────────────────────
rho_df <- read.csv(file.path(INPUT_DIR, RHO_FILE)) %>%
  dplyr::filter(setting == METHOD) %>%
  rename(Species = species)

nic_df <- read.xlsx(file.path(INPUT_DIR, NICHE_FILE), sheet = "OMI_params") %>%
  dplyr::select(Species, Marginality = OMI, Niche_breadth = Tol)

ft_df <- read.xlsx(file.path(INPUT_DIR, FT_FILE), sheet = "FuncTrait32") %>%
  dplyr::select(TaxID, Abbrevration, ESD, Biovolume, Spincule, Colony,
                Speed_max, Trophic_type, Cyst, Toxin) %>%
  rename(Species = Abbrevration) %>%
  left_join(nic_df, by = "Species") %>%
  mutate(
    Toxin         = if_else(Toxin == "None", 0L, 1L),
    log_ESD       = log10(ESD),
    log_Biovolume = log10(Biovolume)
  ) %>%
  dplyr::select(-ESD, -Biovolume) %>%
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
    z_rho = fisher_z(rho),
    z_se  = fisher_z_se(n)
  )

# final input data
mydf <- rho_df %>%
  dplyr::select(Species, setting, ACH, centre, dist_type, n, z_rho, z_se) %>%
  left_join(ft_df, by = "Species") %>%
  relocate(c(log_ESD, log_Biovolume), .after = "z_se")

# ── 5. Helper: significance label ─────────────────────────────
sig_label <- function(p) {
  as.character(cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   labels = c("***", "**", "*", "\u2020", "ns")))
}

# ── 7. Univariate OLS analysis ────────────────────────────────

## 7a. Continuous traits: OLS + Spearman cross-check
cont_results <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, all_of(trait)) %>% drop_na()
  n   <- nrow(sub)
  
  sub$x_std <- scale(sub[[trait]])[, 1]
  
  ols  <- lm(z_rho ~ x_std, data = sub)
  beta <- coef(ols)["x_std"]
  ci   <- confint(ols)["x_std", ]
  p    <- summary(ols)$coefficients["x_std", "Pr(>|t|)"]
  
  sp   <- cor.test(sub[[trait]], sub$z_rho, method = "spearman", exact = FALSE)
  
  data.frame(Trait = trait, Type = "continuous", n = n,
             Beta = beta, CI_low = ci[1], CI_high = ci[2],
             p_OLS = p, Spearman_rho = sp$estimate, p_Spearman = sp$p.value,
             row.names = NULL)
})

## 7b. Binary categorical: Spincule / Colony / Toxin
bin_results <- lapply(c("Spincule", "Colony", "Toxin"), function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, all_of(trait)) %>% drop_na()
  sub[[trait]] <- factor(sub[[trait]])
  n   <- nrow(sub)
  
  ols  <- lm(as.formula(paste("z_rho ~", trait)), data = sub)
  key  <- paste0(trait, "1")
  beta <- coef(ols)[key]
  ci   <- confint(ols)[key, ]
  p    <- summary(ols)$coefficients[key, "Pr(>|t|)"]
  
  # Mann-Whitney U test (non-parametric)
  g    <- split(sub$z_rho, sub[[trait]])
  mw_p <- wilcox.test(g[["0"]], g[["1"]], exact = FALSE)$p.value
  
  data.frame(Trait = trait, Type = "binary", n = n,
             Beta = beta, CI_low = ci[1], CI_high = ci[2],
             p_OLS = p, Spearman_rho = NA, p_Spearman = mw_p,
             row.names = NULL)
})

## 7c. Trophic_type — all pairwise comparisons (emmeans + BH)
sub_t <- mydf %>% dplyr::select(z_rho, Trophic_type) %>% drop_na()
sub_t$Trophic_type <- factor(sub_t$Trophic_type)
sub_t$Trophic_type <- relevel(sub_t$Trophic_type, ref = "OPA")
ols_t  <- lm(z_rho ~ Trophic_type, data = sub_t)
kw_p   <- kruskal.test(z_rho ~ Trophic_type, data = sub_t)$p.value

# Estimated Marginal Mean - EMM (Tukey contrast, BH-FDR correct)
emm    <- emmeans(ols_t, ~ Trophic_type)
pairs_t <- contrast(emm, method = "pairwise", adjust = "BH") %>%
  as.data.frame() %>%
  rename(Trait    = contrast,
         Beta     = estimate,
         p_OLS    = p.value) %>%
  mutate(
    Trait      = paste0("Trophic: ", Trait),
    Type       = "categorical",
    n          = nrow(sub_t),
    p_FDR      = p_OLS,          # adjust="BH" already applied above
    Spearman_rho = NA,
    p_Spearman = kw_p
  )

# confint 별도 추출
ci_t <- confint(contrast(emm, method = "pairwise", adjust = "BH")) %>%
  as.data.frame() %>%
  dplyr::select(lower.CL, upper.CL)

pairs_t$CI_low  <- ci_t$lower.CL
pairs_t$CI_high <- ci_t$upper.CL

trophic_results <- pairs_t %>%
  dplyr::select(Trait, Type, n, Beta, CI_low, CI_high,
                p_OLS, p_FDR, Spearman_rho, p_Spearman)

## 7d. Combine & BH-FDR correction
res_df <- bind_rows(cont_results, bin_results, trophic_results) %>%
  mutate(
    p_FDR = p.adjust(p_OLS, method = "BH"),
    sig   = sig_label(p_FDR)
  )

# ── 8. Trophic-controlled OLS (fixed effect) ──────────────────
#  Replaces lmer — Trophic_type treated as fixed covariate
#  Rationale: only 4 trophic levels → insufficient for random effect
#             variance estimation; Trophic_type is itself a trait of interest
ctrl_results <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>% dplyr::select(z_rho, all_of(trait), Trophic_type, Colony) %>% drop_na()
  if (nrow(sub) < 15) return(NULL)
  
  sub$x_std      <- scale(sub[[trait]])[, 1]
  sub$Trophic_type <- relevel(factor(sub$Trophic_type), ref = TROPHIC_REF)
  
  ols  <- lm(z_rho ~ x_std + Trophic_type, data = sub)    # change control variable (= covariate)
  beta <- coef(ols)["x_std"]
  ci   <- confint(ols)["x_std", ]
  p    <- summary(ols)$coefficients["x_std", "Pr(>|t|)"]
  
  data.frame(Trait = trait, Beta_ctrl = beta,
             CI_low_ctrl = ci[1], CI_high_ctrl = ci[2], p_ctrl = p,
             row.names = NULL)
})

ctrl_df <- bind_rows(ctrl_results) %>%
  mutate(
    p_FDR_ctrl = p.adjust(p_ctrl, method = "BH"),
    sig_ctrl   = sig_label(p_FDR_ctrl)
  )

# ── 9. Forest plot ───────────────────────────────────────────
make_forest <- function(data, beta_col, ci_low, ci_high, sig_col,
                        xlim = c(-1, 1)) {
  data <- data %>%
    mutate(color = case_when(
      .data[[ci_low]]  > 0 ~ "pos",
      .data[[ci_high]] < 0 ~ "neg",
      TRUE                 ~ "ns"
    ),
    sig_fill = case_when(
      color == "pos" & .data[[sig_col]] == "*" ~ "pos",
      color == "neg" & .data[[sig_col]] == "*" ~ "neg",
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

p_uni  <- make_forest(res_df, "Beta", "CI_low", "CI_high", "sig", xlim = c(-0.75, 0.5))

p_ctrl <- make_forest(ctrl_df %>% rename(sig = sig_ctrl),
                      "Beta_ctrl", "CI_low_ctrl", "CI_high_ctrl", "sig", xlim = c(-0.2, 0.2))

forest_combined <- p_uni + p_ctrl +
  plot_annotation(
    title = sprintf("Functional Traits vs. ACH Support (%s, Fisher's Z)", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_forest.png", METHOD)),
       forest_combined, width = 10, height = 6, dpi = 300)

# ── 11. Scatter plots ───────────────────────────────────────────
scatter_plots <- lapply(CONT_TRAITS, function(trait) {
  sub <- mydf %>%
    dplyr::select(Species, z_rho, all_of(trait), Trophic_type) %>%
    drop_na()
  r <- res_df %>% filter(Trait == trait)
  
  x_label <- dplyr::case_when(
    trait == "log_ESD"       ~ "log\u2081\u2080 ESD (\u03bcm)",
    trait == "log_Biovolume" ~ "log\u2081\u2080 Biovolume (\u03bcm\u00b3)",
    TRUE                     ~ trait
  )
  
  ggplot(sub, aes(x = .data[[trait]], y = z_rho, color = Trophic_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_smooth(method = "lm", se = TRUE,
                color = "black", linewidth = 1, alpha = 0.15) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = troph_cols, name = "Trophic type") +
    labs(x      = x_label,
         y      = sprintf("Fisher's Z (%s \u03c1)", METHOD),
         title  = sprintf("%s\n\u03b2=%.3f, p_FDR=%.3f %s",
                          trait, r$Beta, r$p_FDR, as.character(r$sig))) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          legend.position = "right", legend.key = element_rect(color = NA, fill = NA),
          axis.text = element_text(color = "black"))
})

scatter_combined <- wrap_plots(scatter_plots, ncol = 3) +
  plot_annotation(
    title = sprintf("Continuous Traits vs. Fisher's Z (%s \u03c1)", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_scatter.png", METHOD)),
       scatter_combined, width = 13, height = 7, dpi = 300)

# ── 12. Categorical box plots ─────────────────────────────────
cat_plots <- lapply(CAT_TRAITS, function(trait) {
  sub <- mydf %>%
    dplyr::select(z_rho, all_of(trait)) %>%
    drop_na() %>%
    mutate(across(all_of(trait), as.character))
  
  fill_cols <- if (trait == "Trophic_type") {troph_cols} else {bin_cols}
  
  ggplot(sub, aes(x = .data[[trait]], y = z_rho, fill = .data[[trait]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = fill_cols, guide = "none") +
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

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_categorical.png", METHOD)),
       cat_combined, width = 8, height = 8, dpi = 300)

# ── 13. Save results ──────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Univariate")
addWorksheet(wb, "TrophicControlled")
writeData(wb, "Univariate", res_df)
writeData(wb, "TrophicControlled", ctrl_df)
saveWorkbook(wb, file.path(OUT_DIR, sprintf("FT_ACH_%s_results.xlsx", METHOD)), overwrite = TRUE)

cat("\nDone. All outputs saved to:\n", OUT_DIR, "\n")
