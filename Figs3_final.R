## ── FigS3_from_files.R ─────────────────────────────────────────────────────
##  Sensitivity of the abundance–distance relationship to temporal
##  autocorrelation correction — built directly from saved output files.
##
##  Input files (unchanged column structure):
##    ACH_sp_results_corrected.csv   — species x setting: rho_adj, p_adj,
##                                      p_eff_adj, boot_ci_low/high, ACH, ...
##    TableS6_Bayesian_sp_E5_Base.xlsx — species-level beta1 (base model)
##    TableS6_Bayesian_sp_E5_AR1.xlsx  — species-level beta1 (AR(1) model)
##    LOO_AR1_comparison.csv         — setting-level elpd_diff, se_diff
## ────────────────────────────────────────────────────────────────────────────

pkgs <- c("readr", "readxl", "dplyr", "tidyr", "stringr", "ggplot2", "patchwork")
invisible(lapply(pkgs, library, character.only = TRUE))

INPUT_DIR      <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology/ACH/Output_260805" 
OUTPUT_DIR     <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology/ACH/Output_260805" 
FIG3_SETTING <- "E5"
ALPHA       <- 0.05

ACH_cols <- c(
  "Supported"       = "#38D5BE",
  "Opposite"        = "#FF6F61",
  "Not significant" = "#BDBDBD"
)

## ── 1. Load ──────────────────────────────────────────────────────────────
ach_res <- read_csv(file.path(INPUT_DIR, "ACH_sp_results_corrected.csv"), show_col_types = FALSE)
loo_tab <- read_csv(file.path(INPUT_DIR, "LOO_AR1_comparison.csv"), show_col_types = FALSE)
base_df <- read_excel(file.path(INPUT_DIR, "TableS6_Bayesian_sp_E5_Base.xlsx"), sheet = "Table_S6")
ar1_df  <- read_excel(file.path(INPUT_DIR, "TableS6_Bayesian_sp_E5_AR1.xlsx"), sheet = "Table_S6")

## ── 2. Panel A data: species x method significance (setting E5) ──────────
ach_res <- ach_res %>%
  dplyr::filter(setting == FIG3_SETTING) %>%
  dplyr::select(species, rho_adj, pval, p_eff, p_boot) %>%
  dplyr::mutate(p_sig    = case_when(pval <= ALPHA & rho_adj < 0 ~ "Supported", 
                                     pval <= ALPHA & rho_adj > 0 ~ "Opposite",
                                     TRUE ~ "Not significant"),
                eff_sig  = case_when(pval <= ALPHA & rho_adj < 0 ~ "Supported", 
                                     pval <= ALPHA & rho_adj > 0 ~ "Opposite",
                                     TRUE ~ "Not significant"),
                boot_sig = case_when(pval <= ALPHA & rho_adj < 0 ~ "Supported", 
                                     pval <= ALPHA & rho_adj > 0 ~ "Opposite",
                                     TRUE ~ "Not significant")
                ) %>%
  dplyr::arrange(rho_adj) %>%
  dplyr::mutate(species = factor(species, levels = species))

ach_res_long <- ach_res %>%
  dplyr::rename(Raw = pval, ESS = p_eff, MBB = p_boot) %>%
  pivot_longer(Raw:MBB, names_to = "method", values_to = "pvalue") %>%
  dplyr::mutate(method = factor(method, levels = c("MBB","ESS","Raw")))

p_S3A <- ggplot(ach_res_long, aes(x = species, y = method, fill = pvalue)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", pvalue)), colour = "white", angle = 90, size = 4) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans = "sqrt") +
  labs(x = "Species ordered by rho",y = NULL, title = NULL) +
  coord_fixed(ratio = 3) +
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_blank(),
    axis.text.x  = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5), 
    axis.ticks   = element_blank(),
    legend.position = "right",
    plot.title   = element_text(face = "bold", size = 16)
  )

## ── 3. Panel B data: base vs AR1 species-level beta1 ──────────────────────
parse_ci <- function(x) {
  m <- stringr::str_match(x, "\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]")
  tibble::tibble(lo = as.numeric(m[, 2]), hi = as.numeric(m[, 3]))
}

# Load Table S6 from base model
base_sp <- base_df %>%
  dplyr::filter(Species != "Community (fixed effect)") %>%
  dplyr::bind_cols(parse_ci(.$beta1_95CI)) %>%
  dplyr::select(Species, beta1_median_base = beta1_median, beta1_lo_base = lo, beta1_hi_base = hi)

# Load Table S6 from AR(1) model
ar1_sp <- ar1_df %>%
  dplyr::filter(Species != "Community (fixed effect)") %>%
  dplyr::bind_cols(parse_ci(.$beta1_95CI)) %>%
  dplyr::select(Species, beta1_median_ar1 = beta1_median, beta1_lo_ar1 = lo, beta1_hi_ar1 = hi)

beta_comp <- base_sp %>%
  dplyr::inner_join(ar1_sp, by = "Species") %>%
  dplyr::left_join(ach_res %>% dplyr::select(species, eff_sig) %>% dplyr::rename(Species = species),
                   by = "Species") %>%
  dplyr::mutate(ACH = eff_sig) %>%
  dplyr::select(-eff_sig)

r_val <- round(cor(beta_comp$beta1_median_base, beta_comp$beta1_median_ar1), 2)
lims  <- range(c(beta_comp$beta1_lo_base, beta_comp$beta1_hi_base,
                 beta_comp$beta1_lo_ar1,  beta_comp$beta1_hi_ar1))

p_S3B <- ggplot(beta_comp, aes(x = beta1_median_base, y = beta1_median_ar1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_errorbar(aes(ymin = beta1_lo_ar1, ymax = beta1_hi_ar1), width = 0, colour = "grey70", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = beta1_lo_base, xmax = beta1_hi_base), height = 0, colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = ACH), shape = 21, size = 5, colour = "black", stroke = 0.3) +
  scale_fill_manual(values = ACH_cols, name = "ACH support\n(ESS)", drop = FALSE) +
  coord_equal(xlim = lims, ylim = lims) +
  annotate("text", x = lims[1], y = lims[2], label = paste0("r = ", r_val),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = expression(beta[1]~"(base)"),
       y = expression(beta[1]~"(AR(1))"),
       title = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black"),
        axis.line = element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        plot.margin = margin(t = 5.5, r = 40, b = 5.5, l = 5.5)
        )

## ── 4. Panel C data: LOO elpd_diff, restricted to settings actually used ──
SETTINGS <- c("E1","E2","E4","E5","E7","E8")

p_S3C <- loo_tab %>%
  dplyr::filter(setting %in% SETTINGS) %>%
  dplyr::mutate(setting = factor(setting, levels = SETTINGS)) %>%
  ggplot(aes(x = elpd_diff, y = rev(setting))) +
  #geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = elpd_diff - 1.96 * se_diff, xmax = elpd_diff + 1.96 * se_diff),
                 height = 0, colour = "#2166AC", linewidth = 0.9) +
  geom_point(colour = "#2166AC", size = 5) +
  scale_x_continuous(limits = c(-800, -300), breaks = seq(-800, -300, by = 100)) +
  labs(x = expression(Delta*"ELPD (AR(1)"~-~"base)"), y = NULL, title = NULL) +
  coord_fixed(ratio = 130) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black"),
        axis.line = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 20)
        )

## ── 5. Combine & export ────────────────────────────────────────────────────
p_S3 <- p_S3A / (p_S3B | p_S3C) + 
  plot_layout(widths = c(1, 1), heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A", 
                  theme = theme(plot.tag = element_text(face = "bold", size = 16, colour = "black")))

ggsave(file.path(OUT_DIR, "FigS3_autocorr_sensitivity.pdf"), p_S3, width = 9, height = 9, dpi = 600)
ggsave(file.path(OUT_DIR, "FigS3_autocorr_sensitivity.tiff"), p_S3, width = 9, height = 9, dpi = 600)
