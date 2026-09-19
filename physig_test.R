# ============================================================
#  PHYLOGENETIC SIGNAL & OLS vs PGLS COMPARISON
#
#  1. NCBI taxonomy 기반 수동 계통수 구성 (phylostratr::ncbi_tree)
#  2. Blomberg's K 계산 (phytools::phylosig)
#  3. OLS vs PGLS 비교 (nlme::gls + ape::corBrownian / corPagel)
#  4. 시각화 및 결과 저장
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ── 0. Packages ───────────────────────────────────────────────
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr", "stringr",
  "ggplot2", "patchwork", "emmeans",
  "ape", "phytools", "phylostratr", "nlme"
)

if (length(pkgs)) install.packages(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. Paths ──────────────────────────────────────────────────
INPUT_DIR  <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology"
RHO_FILE   <- "ACH/Output_260804/ACH_sp_results_corrected.csv"
FT_FILE    <- "dino_functraits.xlsx"
NICHE_FILE <- "OMI/Output_260804/OMI_results_260804.xlsx"

WORK_DIR <- file.path(INPUT_DIR, "ACH")
OUT_DIR  <- file.path(WORK_DIR, paste0("Output_FT", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 2. Global settings ────────────────────────────────────────
METHOD      <- "E5"
TROPHIC_REF <- "CM"
CONT_TRAITS <- c("log_Biovolume", "Speed_max", "Marginality", "Niche_breadth")
CAT_TRAITS  <- c("Spicule", "Colony", "Trophic_type", "Cyst")

troph_cols <- c(CM = "#ff7f0e", pSNCM = "#d62728", HET = "#1f77b4", OPA = "#2ca02c")

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

taxid <- ft_df %>% dplyr::select(TaxID, Species)

# ── 4. Fisher's Z transformation ──────────────────────────────
fisher_z    <- function(r) atanh(pmin(pmax(r, -0.9999), 0.9999))
fisher_z_se <- function(n) 1 / sqrt(n - 3)

rho_df <- rho_df %>%
  mutate(
    z_rho = fisher_z(rho),
    z_se  = fisher_z_se(n),
    # ── temporal-AC correction (Method 2: effective-N, Pyper & Peterman 1998) ──
    # rho itself is unchanged; only its precision (SE) is corrected for
    # autocorrelation in the 411-day series. z_se uses the nominal n and is
    # kept as-is for comparability with the original analysis.
    z_se_eff = fisher_z_se(n_eff)
  )

# final input data
mydf <- rho_df %>%
  dplyr::select(Species, setting, ACH, ACH_effN, ACH_boot, ACH_AR1,
                centre, dist_type, n, n_eff, z_rho, z_se, z_se_eff) %>%
  left_join(ft_df, by = "Species") %>%
  relocate(c(log_Biovolume), .after = "z_se")

# ── 5. Helper: significance label ─────────────────────────────
sig_label <- function(p) {
  as.character(cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   labels = c("***", "**", "*", "\u2020", "ns")))
}

# ── 6. Phylogenetic signal (Blomberg's K)────────────────────────────────────
#  phytools::phylosig(tree, x, method="K", test=TRUE, nsim=999)
#  연속형 형질 및 응답변수(z_rho) 모두 계산
#  해석 기준: K≈1 (BM), K>1 (계통 보존적), K<1 (약한 신호), p<0.05 (유의)
tree   <- ncbi_tree(taxid$TaxID)
tree   <- compute.brlen(tree, method = "Grafen", power = 1)

tip_name <- setNames(
  as.character(taxid$Species),   # 값: 약어
  as.character(taxid$TaxID)      # 키: TaxID 문자열
)

tree$tip.label <- dplyr::recode(tree$tip.label, !!!tip_name)

k_results <- lapply(c(CONT_TRAITS, "z_rho"), function(trait) {
  
  # drop_na()로 해당 형질의 NA 종 자동 제외
  sub <- mydf %>%
    dplyr::select(Species, all_of(trait)) %>%
    drop_na() %>%
    dplyr::filter(Species %in% tree$tip.label)
  # ↑ arrange(tree$tip.label) 제거 — 길이 불일치 원인
  
  if (nrow(sub) < 8) return(NULL)
  
  # NA 제거 후 남은 종만으로 계통수 pruning
  tree_sub <- ape::keep.tip(tree, sub$Species)
  
  # pruning된 tree 순서에 맞게 x_vec 정렬
  x_vec <- setNames(
    sub[[trait]][match(tree_sub$tip.label, sub$Species)],
    tree_sub$tip.label
  )
  
  k_res <- phytools::phylosig(
    tree_sub, x_vec,
    method = "K", test = TRUE, nsim = 999
  )
  
  data.frame(Trait = trait, K = k_res$K, p_K = k_res$P,
             n = length(x_vec), row.names = NULL)
})

k_df <- bind_rows(k_results) %>%
  mutate(
    sig_K = sig_label(p_K),
    interpretation = dplyr::case_when(
      K >  1 & p_K < 0.05 ~ "Phylogenetically conserved (K > 1)",
      K <= 1 & p_K < 0.05 ~ "Weaker than BM (K < 1, sig.)",
      TRUE                 ~ "No significant phylogenetic signal"
    )
  )

# ── 7. OLS vs PGLS 비교 ─────────────────────────────────────
#  두 가지 공분산 구조를 비교:
#  - corBrownian: λ = 1 고정 (완전한 BM 계통 구조)
#  - corPagel   : λ 자동 추정 (0 = OLS, 1 = BM; 데이터 기반 최적값)

pgls_results <- lapply(CONT_TRAITS, function(trait) {
  
  sub <- mydf %>%
    dplyr::select(Species, z_rho, all_of(trait)) %>%
    drop_na() %>%
    dplyr::filter(Species %in% tree$tip.label) %>%
    as.data.frame()
  
  if (nrow(sub) < 10) return(NULL)
  
  sub$x_std <- scale(sub[[trait]])[, 1]
  tree_sub  <- ape::keep.tip(tree, sub$Species)
  
  # OLS
  ols    <- lm(z_rho ~ x_std, data = sub)
  b_ols  <- coef(ols)["x_std"]
  ci_ols <- confint(ols)["x_std", ]
  p_ols  <- summary(ols)$coefficients["x_std", "Pr(>|t|)"]
  
  # ── temporal-AC correction (Method 2: effective-N weighted OLS) ──────────
  # Complements the phylogenetic (PGLS) correction below: this instead
  # down-weights species whose rho is estimated from a time series with
  # strong autocorrelation (small n_eff relative to n).
  sub_effN <- mydf %>%
    dplyr::select(Species, z_rho, z_se_eff, all_of(trait)) %>%
    drop_na() %>%
    dplyr::filter(Species %in% tree$tip.label) %>%
    as.data.frame()
  sub_effN$x_std <- scale(sub_effN[[trait]])[, 1]
  sub_effN$w     <- 1 / (sub_effN$z_se_eff^2)
  ols_effN   <- lm(z_rho ~ x_std, data = sub_effN, weights = w)
  b_ols_effN <- coef(ols_effN)["x_std"]
  ci_ols_effN<- confint(ols_effN)["x_std", ]
  p_ols_effN <- summary(ols_effN)$coefficients["x_std", "Pr(>|t|)"]
  
  # Useful function - extract CI from GLS
  extract_gls <- function(mod) {
    if (is.null(mod)) {
      return(list(beta = NA_real_, se = NA_real_,
                  ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_))
    }
    tbl    <- summary(mod)$tTable
    b      <- tbl["x_std", "Value"]
    se     <- tbl["x_std", "Std.Error"]
    p      <- tbl["x_std", "p-value"]
    df_res <- mod$dims$N - mod$dims$p   # residual df
    tc     <- qt(0.975, df = df_res)
    list(beta  = b,
         se    = se,
         ci_lo = b - tc * se,
         ci_hi = b + tc * se,
         p     = p)
  }
  
  # PGLS: Brownian motion (λ = 1 고정)
  pgls_bm <- tryCatch(
    nlme::gls(
      z_rho ~ x_std,
      data        = sub,
      correlation = ape::corBrownian(phy = tree_sub, form = ~Species),
      method      = "ML"
    ),
    error = function(e) {
      message(sprintf("  BM-PGLS failed for %s: %s", trait, e$message))
      NULL
    }
  )
  
  # PGLS: Pagel's λ (λ 자동 추정, 0~1)
  pgls_pagel <- tryCatch(
    nlme::gls(
      z_rho ~ x_std,
      data        = sub,
      correlation = ape::corPagel(
        value = 1, phy = tree_sub,
        fixed = FALSE, form = ~Species
      ),
      method = "ML"
    ),
    error = function(e) {
      message(sprintf("  Pagel-PGLS failed for %s: %s", trait, e$message))
      NULL
    }
  )
  
  lambda_est <- if (!is.null(pgls_pagel)) {
    as.numeric(coef(pgls_pagel$modelStruct$corStruct, unconstrained = FALSE))
  } else NA_real_
  
  bm    <- extract_gls(pgls_bm)
  pagel <- extract_gls(pgls_pagel)
  
  data.frame(
    Trait            = trait,   n = nrow(sub),
    Beta_OLS         = b_ols,   CI_lo_OLS    = ci_ols[1], CI_hi_OLS    = ci_ols[2], p_OLS        = p_ols,
    Beta_OLS_effN    = b_ols_effN, CI_lo_OLS_effN = ci_ols_effN[1], CI_hi_OLS_effN = ci_ols_effN[2], p_OLS_effN = p_ols_effN,
    Beta_PGLS_BM     = bm$beta, CI_lo_PGLS_BM     = bm$ci_lo,  CI_hi_PGLS_BM     = bm$ci_hi,  p_PGLS_BM    = bm$p,
    Beta_PGLS_Pagel  = pagel$beta, CI_lo_PGLS_Pagel = pagel$ci_lo, CI_hi_PGLS_Pagel = pagel$ci_hi, p_PGLS_Pagel = pagel$p,
    Lambda           = lambda_est,
    Delta_BM         = bm$beta    - b_ols,
    Delta_Pagel      = pagel$beta - b_ols,
    Delta_effN       = b_ols_effN - b_ols,
    row.names = NULL
  )
})

pgls_df <- bind_rows(pgls_results) %>%
  mutate(
    p_FDR_OLS        = p.adjust(p_OLS, method = "BH"),
    p_FDR_OLS_effN   = p.adjust(p_OLS_effN, method = "BH"),
    p_FDR_PGLS_BM    = p.adjust(p_PGLS_BM, method = "BH"),
    p_FDR_PGLS_Pagel = p.adjust(p_PGLS_Pagel, method = "BH"),
    sig_OLS          = sig_label(p_FDR_OLS),
    sig_OLS_effN     = sig_label(p_FDR_OLS_effN),
    sig_BM           = sig_label(p_FDR_PGLS_BM),
    sig_Pagel        = sig_label(p_FDR_PGLS_Pagel)
  )

# ── 8. 시각화 ────────────────────────────────────────────────
## 1) Blomberg's K bar chart
k_plot <- k_df %>% dplyr::filter(Trait %in% CONT_TRAITS) %>%
  mutate(p_group = cut(p_K,
                       breaks = c(-Inf, 0.05, 0.10, Inf),
                       labels = c("p < 0.05", "p < 0.10", "ns")),
         Trait = factor(Trait, levels = CONT_TRAITS))

p_k_bar <- ggplot(k_plot, aes(x = Trait, y = K, fill = p_group)) +
  geom_col(width = 0.6, color = "white") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "gray60") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray60") +
  geom_text(aes(label  = sprintf("K=%.3f\np=%.3f", K, p_K),
                y      = pmax(K, 0) + 0.015),
            vjust = 0, size = 3.2, color = "gray20") +
  scale_fill_manual(
    values = c("p < 0.05" = "#d62728", "p < 0.10" = "#ff7f0e", "ns" = "#bbbbbb"),
    name   = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(
    x = NULL, y = "Blomberg's K"
  ) +
  coord_fixed(ratio = 2) +
  theme_classic(base_size = 14) +
  theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom", legend.key = element_rect(color = NA, fill = NA),
        axis.text = element_text(color = "black"))

## 2) OLS vs PGLS forest plot
pgls_long <- pgls_df %>%
  mutate(across(starts_with("sig_"), as.character)) %>%
  dplyr::select(Trait,
                Beta_OLS, CI_lo_OLS, CI_hi_OLS, p_OLS, p_FDR_OLS,
                Beta_PGLS_BM, CI_lo_PGLS_BM, CI_hi_PGLS_BM, p_PGLS_BM, p_FDR_PGLS_BM,
                Beta_PGLS_Pagel, CI_lo_PGLS_Pagel, CI_hi_PGLS_Pagel, p_PGLS_Pagel, p_FDR_PGLS_Pagel) %>%
  tidyr::pivot_longer(
    cols          = -Trait,
    names_to      = c("Stat", "Model"),
    names_pattern = "(Beta|CI_lo|CI_hi|p|p_FDR)_(OLS|PGLS_BM|PGLS_Pagel)"
  ) %>%
  tidyr::pivot_wider(names_from = Stat, values_from = value) %>%
  mutate(
    Trait = factor(Trait, levels = CONT_TRAITS),
    Sig = sig_label(p_FDR),
    Model = factor(Model,
                   levels = c("OLS", "PGLS_BM", "PGLS_Pagel"),
                   labels = c("OLS", "PGLS (BM)", "PGLS (Pagel \u03bb)")),
    across(c(Beta, CI_lo, CI_hi), as.numeric),
    fill_group = dplyr::if_else(p_FDR < 0.05, as.character(Model), "ns")
  )

p_forest <- ggplot(pgls_long,
                   aes(x = Beta, y = Trait,
                       color = Model, fill = fill_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 linewidth = 0.5, height = 0,
                 position = position_dodge(width = 0.6)) +
  geom_point(pch = 21, size = 4, position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = c("OLS" = "#0F4C81", "PGLS (BM)" = "#9B2335", "PGLS (Pagel \u03bb)" = "#EFC050"),
    name  = "Model"
  ) +
  scale_fill_manual(values = c("OLS" = "#0F4C81", "PGLS (BM)" = "#9B2335", 
                               "PGLS (Pagel \u03bb)" = "#EFC050", "ns" = "white"), guide = "none") +
  xlim(c(-0.4, 0.4)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Standardised \u03b2 (95% CI)", y = NULL) +
  coord_fixed(ratio = 0.5) +
  theme_classic(base_size = 14) +
  theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom", legend.key = element_rect(color = NA, fill = NA),
        axis.text = element_text(color = "black"))

## 3) Delta Beta bar chart
delta_long <- pgls_df %>%
  dplyr::select(Trait, Delta_BM, Delta_Pagel) %>%
  tidyr::pivot_longer(-Trait, names_to = "Model", values_to = "Delta_Beta") %>%
  mutate(Model = dplyr::recode(Model,
                               Delta_BM    = "PGLS (BM) \u2212 OLS",
                               Delta_Pagel = "PGLS (Pagel \u03bb) \u2212 OLS"),
         Trait = factor(Trait, levels = CONT_TRAITS))

p_delta <- ggplot(delta_long,
                  aes(x = Trait, y = Delta_Beta, fill = Model)) +
  geom_col(position = position_dodge(width = 0.6),
           width = 0.5, color = "white") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("PGLS (BM) \u2212 OLS" = "#9B2335",
                               "PGLS (Pagel \u03bb) \u2212 OLS" = "#EFC050"),
                    name = NULL) +
  labs(x = NULL, y = "\u0394\u03b2") +
  coord_fixed(ratio = 10) +
  theme_classic(base_size = 14) +
  theme(panel.background = element_rect(linewidth = 1, color = "black", fill = NA),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom", legend.key = element_rect(color = NA, fill = NA),
        axis.text = element_text(color = "black")
        )

## 4) 패널 결합 및 저장
phylo_panel <- (p_k_bar / p_delta | p_forest) +
  plot_annotation(title = sprintf("Phylogenetic Signal in Functional traits [%s]", METHOD),
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  ) +
  plot_layout(heights = c(1, 1), widths = c(2, 2)) &
  theme(plot.margin = margin(1, 1, 1, 1))

ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_phylosig.pdf", METHOD)),
       phylo_panel, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("FT_ACH_%s_phylosig.tiff", METHOD)),
       phylo_panel, width = 12, height = 8, dpi = 300)

# ── 9. 결과 저장 ─────────────────────────────────────────────
wb <- createWorkbook()

addWorksheet(wb, "Blomberg_K")
addWorksheet(wb, "OLS_vs_PGLS")
writeData(wb, "Blomberg_K", k_df)
writeData(wb, "OLS_vs_PGLS", pgls_df)
saveWorkbook(wb, file.path(OUT_DIR, sprintf("phylsig_%s_results.xlsx", METHOD)), overwrite = TRUE)
