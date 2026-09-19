## ── Fig 2: Final publication figure (single setting) ───────────────────────────
##
##  Panel A : Species-level Spearman ρ bar chart (sorted)                — model-independent
##  Panel B : Community-level Bayesian posterior prediction curve
##            + x-axis rug for data density                              — per Bayesian model
##  Panel C : Species-level posterior prediction spaghetti plot
##            (lines coloured by OLS-adjusted R²)                        — per Bayesian model
##
##  Table S6: Species-level Bayesian β₁ / β₂ + P(β₁<0) + OLS adj.R²
##            written to Excel (openxlsx)                                — per Bayesian model
##
##  This version runs the SAME panel-B/C/Table-S6 pipeline twice, once for the
##  base quadratic hierarchical model (bayes_fits) and once for the AR(1)
##  residual-correlation model (bayes_fits_ar1), WITHOUT merging the two model
##  fitting procedures. Two full sets of outputs are written
##  (suffix "_Base" and "_AR1"); the user assembles the final 3-panel figure
##  manually by choosing which B/C panels to pair with Panel A.
##
##  Prerequisites (must exist in environment before sourcing this file):
##    ach_res        — species-level ACH results (rho_adj, ACH, p_adj columns)
##    model_df       — long-format data (setting, species, dist_signed, log_abund)
##    bayes_fits     — named list of brms fit objects (base model, names = setting codes)
##    bayes_fits_ar1 — named list of brms fit objects (AR(1) model, names = setting codes)
##    OUT_DIR        — output directory path
##    FIG2_SETTING   — e.g. "E4"  (set once before sourcing)
##
##  Required packages: ggplot2, patchwork, dplyr, tidyr, tibble, purrr,
##                     brms, posterior, openxlsx
## ─────────────────────────────────────────────────────────────────────────────

# ── 1. environment setting ───────────────────────────────────────────────────────
# Choosing setting
FIG2_SETTING <- "E5"

# loading packages
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr",
  "ggplot2", "patchwork",
  "brms", "posterior"
)

if (length(pkgs)) install.packages(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# check presence of input data (both Bayesian model lists required now)
stopifnot(
  exists("ach_res"),
  exists("model_df"),
  exists("bayes_fits"),
  exists("bayes_fits_ar1"),
  exists("OUT_DIR"),
  exists("FIG2_SETTING")
)

# color palette
ACH_cols <- c(
  "Supported"         = "#38D5BE",
  "Opposite"          = "#FF6F61",
  "Not significant"   = "#BDBDBD",
  "Insufficient data" = "#F0F0F0"
)


# ── 2. Subset to chosen setting (model-independent) ────────────────────────
ach_s <- ach_res %>%
  dplyr::filter(setting == FIG2_SETTING) %>%
  dplyr::arrange(rho_adj) %>%
  dplyr::mutate(
    species = factor(species, levels = species),
    ACH     = factor(ACH_effN, levels = names(ACH_cols))
  )

model_s  <- model_df %>% dplyr::filter(setting == FIG2_SETTING)
is_mg    <- unique(model_s$is_margin)
sp_list  <- levels(ach_s$species)          # ordered by rho_adj (used in panel A)
dist_seq <- seq(0, 1, length.out = 120)

# one-sample Wilcoxon for panel A annotation
rho_vec   <- ach_s$rho_adj[!is.na(ach_s$rho_adj)]
wilcox_p  <- wilcox.test(rho_vec, mu = 0)$p.value
wilcox_lab <- dplyr::case_when(
  wilcox_p < 0.001 ~ "***", wilcox_p < 0.01 ~ "**",
  wilcox_p < 0.05  ~ "*",   wilcox_p < 0.10 ~ "†",
  TRUE             ~ "ns"
)
n_sup <- sum(ach_s$ACH == "Supported",  na.rm = TRUE)
n_opp <- sum(ach_s$ACH == "Opposite",   na.rm = TRUE)
n_tot <- sum(!is.na(ach_s$rho_adj))

# Helper: signed distance (margin settings use negative direction)
make_dist_signed <- function(d) if (is_mg) -d else d

# OLS adj.R² per species (for spaghetti colour encoding) — model-independent
ols_r2 <- model_s %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    adj_r2 = {
      fit_lm <- lm(log_abund ~ dist_signed, data = dplyr::cur_data())
      summary(fit_lm)$adj.r.squared
    },
    .groups = "drop"
  )

# rug data: observed dist_std values (all species pooled) — model-independent
rug_df <- model_s %>%
  dplyr::mutate(
    dist_std_plot = abs(dist_signed) / max(abs(dist_signed), na.rm = TRUE)) %>%
  dplyr::select(dist_std_plot, log_abund)


# ── 3. Panel A: species Spearman rho bar chart (shared across models) ───────
p_A <- ggplot(ach_s, aes(x = rho_adj, y = species, fill = ACH_effN)) +
  geom_col() +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = ACH_cols, drop = FALSE, name = "ACH support") +
  scale_x_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, 0.4)) +
  labs(x = expression("Spearman "*rho), y = NULL, title = "A") +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "inside",
    legend.justification = c(0.01, 0.98),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line        = element_blank(),
    plot.title       = element_text(face = "bold", size = 14)
  )


# ── 4. Reusable pipeline: Bayesian panels B/C + combined figure + Table S6 ──
##  model_fit_list : named list of brms fits (names = setting codes)
##  model_label    : short label used in filenames / messages, e.g. "Base", "AR1"
##  has_ar1        : TRUE if the fit includes ar(time = CountDay, gr = species, p = 1).
##                    When TRUE, CountDay is supplied to newdata as a dummy
##                    placeholder (required by validate_data()) and the
##                    autocorrelation term is excluded from prediction via
##                    incl_autocor = FALSE — the AR(1) structure only concerns
##                    residual correlation between *observed* time points and
##                    is not meaningful for the synthetic dist_seq prediction grid.
build_bayes_outputs <- function(model_fit_list, model_label, has_ar1 = FALSE) {
  
  if (!FIG2_SETTING %in% names(model_fit_list))
    stop("No Bayesian fit found for setting: ", FIG2_SETTING,
         " (model: ", model_label, ")")
  
  fit <- model_fit_list[[FIG2_SETTING]]
  
  # dummy CountDay: only needed to satisfy validate_data() when has_ar1 = TRUE;
  # value is irrelevant because incl_autocor = FALSE drops the AR(1) term
  dummy_countday <- if (has_ar1) rep(1L, length(dist_seq)) else NULL
  
  # 4a. Community-level fixed-effect prediction (re_formula = NA)
  nd_comm <- tibble::tibble(
    species     = sp_list[1],          # reference species (marginalised below)
    dist_std    = dist_seq,
    dist_signed = make_dist_signed(dist_seq)
  )
  if (has_ar1) nd_comm$CountDay <- dummy_countday
  
  # incl_autocor = FALSE is only passed for AR1 fits (base model has no autocor term)
  draws_comm <- if (has_ar1) {
    brms::posterior_epred(fit, newdata = nd_comm, re_formula = NA,
                          incl_autocor = FALSE)
  } else {
    brms::posterior_epred(fit, newdata = nd_comm, re_formula = NA)
  }
  q_comm     <- apply(draws_comm, 2,
                      quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  comm_curve <- tibble::tibble(
    dist_std = dist_seq,
    y_med    = q_comm[2, ],
    y_low    = q_comm[1, ],
    y_high   = q_comm[3, ]
  )
  
  # 4b. Species-level posterior predictions (re_formula = full)
  sp_pred_list <- purrr::map(sp_list, function(sp) {
    nd_sp <- tibble::tibble(
      species     = sp,
      dist_std    = dist_seq,
      dist_signed = make_dist_signed(dist_seq)
    )
    if (has_ar1) nd_sp$CountDay <- dummy_countday
    
    dr <- if (has_ar1) {
      brms::posterior_epred(
        fit, newdata = nd_sp,
        re_formula = ~ (1 + dist_signed + I(dist_signed^2) | species),
        allow_new_levels = FALSE,
        incl_autocor = FALSE
      )
    } else {
      brms::posterior_epred(
        fit, newdata = nd_sp,
        re_formula = ~ (1 + dist_signed + I(dist_signed^2) | species),
        allow_new_levels = FALSE
      )
    }
    q <- apply(dr, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    tibble::tibble(
      species  = sp,
      dist_std = dist_seq,
      y_med    = q[2, ],
      y_low    = q[1, ],
      y_high   = q[3, ]
    )
  })
  sp_pred_df <- dplyr::bind_rows(sp_pred_list) %>%
    dplyr::left_join(ols_r2, by = "species")
  
  # 4c. Species-level β₁, β₂ extraction for Table S6
  #     Column names in brms draws: b_dist_signed, b_Idist_signedE2
  #     (verify with posterior::variables(fit) if names differ; AR(1) models
  #      add ar[...] columns but the fixed/random effect names are unchanged)
  draws_all  <- posterior::as_draws_df(fit)
  b1_pop     <- draws_all[["b_dist_signed"]]
  b2_pop     <- draws_all[["b_Idist_signedE2"]]
  
  sp_coef_df <- purrr::map_dfr(sp_list, function(sp) {
    re_b1_col <- paste0("r_species[", sp, ",dist_signed]")
    re_b2_col <- paste0("r_species[", sp, ",Idist_signedE2]")
    
    if (!re_b1_col %in% names(draws_all)) {
      warning("[", model_label, "] Column not found: ", re_b1_col,
              "\n  Run posterior::variables(fit) to check exact names.")
      return(NULL)
    }
    
    b1_sp <- b1_pop + draws_all[[re_b1_col]]
    b2_sp <- b2_pop + draws_all[[re_b2_col]]
    
    tibble::tibble(
      Species        = sp,
      n_obs          = sum(model_s$species == sp),
      beta1_median   = round(median(b1_sp), 2),
      beta1_CI_low   = round(quantile(b1_sp, 0.025), 2),
      beta1_CI_high  = round(quantile(b1_sp, 0.975), 2),
      beta2_median   = round(median(b2_sp), 2),
      beta2_CI_low   = round(quantile(b2_sp, 0.025), 2),
      beta2_CI_high  = round(quantile(b2_sp, 0.975), 2),
      P_beta1_neg    = round(mean(b1_sp < 0), 3)   # P(β₁ < 0): ACH support probability
    )
  })
  
  sp_coef_df <- sp_coef_df %>%
    dplyr::mutate(P_beta1_neg = round(P_beta1_neg * 100, 1)) %>%
    dplyr::left_join(ach_s %>%
                       dplyr::select(species, rho_adj, p_adj, ACH) %>%
                       dplyr::mutate(species = as.character(species)),
                     by = c("Species" = "species")
    ) %>%
    dplyr::rename(
      Spearman_rho  = rho_adj,
      p_FDR         = p_adj,
      ACH_support   = ACH
    ) %>%
    dplyr::arrange(Spearman_rho)
  
  # ── Panel B: community-level Bayesian curve + rug ──────────────────────
  p_B <- ggplot() +
    geom_ribbon(data  = comm_curve,
                aes(x = dist_std, ymin = y_low, ymax = y_high),
                fill = "#C05B3C", alpha = 0.20) +
    geom_line(data  = comm_curve,
              aes(x = dist_std, y = y_med, colour = "#C05B3C"),
              linewidth = 1.2) +
    geom_rug(data = rug_df,
             aes(x = dist_std_plot),
             side = "b", colour = "grey40", alpha = 0.25, linewidth = 0.3,
             length = unit(0.04, "npc")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(1.5, 5.5), breaks = seq(2, 5, 1)) +
    labs(
      x     = "Standardised distance (0–1)",
      y     = expression(log[10](abundance + 1)),
      title = "B"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.position = "none",
      axis.line       = element_blank(),
      plot.title      = element_text(face = "bold", size = 14)
    )
  
  # ── Panel C: species-level spaghetti ────────────────────────────────────
  p_C <- ggplot(sp_pred_df %>%
                  dplyr::left_join(sp_coef_df %>% dplyr::select(Species, P_beta1_neg),
                                   by = c("species" = "Species")),
                aes(x = dist_std, y = y_med,
                    group = species, colour = P_beta1_neg)
  ) +
    geom_line(linewidth = 0.55, alpha = 0.85) +
    scale_colour_gradient2(
      limits   = c(0,100),
      low      = "#FF6F61",
      mid      = "grey80",
      high     = "#38D5BE",
      midpoint = 50,
      name     = expression("Prob"~"(%)"),
      guide    = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 4)
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(1.5, 5.5), breaks = seq(1, 5, 1)) +
    labs(
      x   = "Standardised distance (0–1)",
      y   = expression(log[10](abundance + 1)),
      title = "C"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = element_blank(),
      legend.position = "inside",
      legend.justification = c(0.05, 0.98),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # ── Combine & save (A + B + C for THIS model, kept as a reference figure;
  #    final manuscript panel selection is done by the user) ───────────────
  fig2_model <- (p_A + p_B + p_C) +
    patchwork::plot_layout(widths = c(1, 1, 1)) &
    theme(plot.margin = margin(4, 4, 4, 4))
  
  ggplot2::ggsave(
    file.path(OUT_DIR, sprintf("Fig2_final_%s_%s.pdf", FIG2_SETTING, model_label)),
    fig2_model, width = 12, height = 4.5, dpi = 600
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, sprintf("Fig2_final_%s_%s.tiff", FIG2_SETTING, model_label)),
    fig2_model, width = 12, height = 4.5, dpi = 300
  )
  message(sprintf("✓ Figure 2 [%s] saved (A + B + C layout)", model_label))
  
  # ── Table S6 for THIS model ──────────────────────────────────────────────
  b1_pop_sum <- tibble::tibble(
    Species       = "Community (fixed effect)",
    n_obs         = nrow(model_s),
    beta1_median  = round(median(b1_pop), 2),
    beta1_CI_low  = round(quantile(b1_pop, 0.025), 2),
    beta1_CI_high = round(quantile(b1_pop, 0.975), 2),
    beta2_median  = round(median(b2_pop), 2),
    beta2_CI_low  = round(quantile(b2_pop, 0.025), 2),
    beta2_CI_high = round(quantile(b2_pop, 0.975), 2),
    P_beta1_neg   = round(mean(b1_pop < 0) * 100, 1),
    Spearman_rho  = NA_real_,
    p_FDR         = NA_real_,
    ACH_support   = NA_character_
  )
  
  table_s6 <- dplyr::bind_rows(b1_pop_sum, sp_coef_df) %>%
    dplyr::mutate(
      beta1_95CI  = sprintf("[%.2f, %.2f]", beta1_CI_low, beta1_CI_high),
      beta2_95CI  = sprintf("[%.2f, %.2f]", beta2_CI_low, beta2_CI_high),
    ) %>%
    dplyr::select(
      Species, n_obs,
      beta1_median, beta1_95CI,
      beta2_median, beta2_95CI,
      P_beta1_neg,
      Spearman_rho, p_FDR, ACH_support
    )
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Table_S6")
  
  hs <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill         = "#D9E1F2",
    border         = "Bottom",
    wrapText       = TRUE
  )
  
  openxlsx::writeDataTable(
    wb, "Table_S6",
    x           = table_s6,
    tableStyle  = "TableStyleLight2",
    headerStyle = hs
  )
  
  out_xlsx <- file.path(OUT_DIR, sprintf("TableS6_Bayesian_sp_%s_%s.xlsx", FIG2_SETTING, model_label))
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message(sprintf("✓ Table S6 [%s] saved: %s", model_label, out_xlsx))
  
  list(
    model_label = model_label,
    comm_curve  = comm_curve,
    sp_pred_df  = sp_pred_df,
    sp_coef_df  = sp_coef_df,
    table_s6    = table_s6,
    panel_B     = p_B,
    panel_C     = p_C,
    figure      = fig2_model
  )
}


# ── 5. Run pipeline for both Bayesian models ─────────────────────────────────
res_base <- build_bayes_outputs(bayes_fits,     "Base", has_ar1 = FALSE)
res_ar1  <- build_bayes_outputs(bayes_fits_ar1, "AR1",  has_ar1 = TRUE)


# ── 6. Quick summary print for both models ────────────────────────────────────
cat("\n── Table S6 preview [Base model] ──\n")
print(res_base$table_s6, n = 30)

cat("\n── Table S6 preview [AR1 model] ──\n")
print(res_ar1$table_s6, n = 30)

cat(sprintf(
  "\n✓ Panel A is model-independent; Panel B/C generated separately for Base and AR1.\n  Combine manually (e.g. p_A + res_base$panel_B + res_ar1$panel_C) for the final 3-panel figure.\n"
))
