# ============================================================
#  Functional Traits vs. ACH Support
#  Correlation among functional traits for checking robustness 
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ── 0. Packages ───────────────────────────────────────────────
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr", "stringr",
  "Hmisc",  # rcorr()
  "ggplot2", "pheatmap"
)

if (length(pkgs)) install.packages(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. Paths ──────────────────────────────────────────────────
INPUT_DIR  <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology"
RHO_FILE   <- "ACH/Output_260502/ACH_sp_results.csv"
FT_FILE    <- "dino_functraits_origin.xlsx"
NICHE_FILE <- "OMI/Output_260506/OMI_results_260506.xlsx"

WORK_DIR <- file.path(INPUT_DIR, "ACH")
OUT_DIR  <- file.path(WORK_DIR, paste0("Output_", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 2. Global settings ────────────────────────────────────────
METHOD      <- "E2"
CONT_TRAITS <- c("log_ESD", "log_Biovolume", "log_Speed_max", "Marginality", "Niche_breadth")
BIN_TRAITS  <- c("Spincule", "Colony", "Cyst", "Toxin")
CAT_TRAITS  <- "Trophic_type"

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
    log_Biovolume = log10(Biovolume),
    log_Speed_max = log10(Speed_max)
  ) %>%
  dplyr::select(-ESD, -Biovolume, -Speed_max) %>%
  dplyr::filter(
    Species %in% rho_df$Species,
    !(Species == "Nsci" & Trophic_type == "eSNCM"),
    !(Species == "Poly" & Trophic_type == "CM")
  )

# ── 4. Calculate correlation between traits ────────────────────
# 1) continuous x continuous traits: Spearman
cont_df <- ft_df %>%
  mutate(Trophic_num = recode(Trophic_type, CM=1, pSNCM=2, HET=3, OPA=4)) %>%
  dplyr::select(CONT_TRAITS, BIN_TRAITS, Trophic_num) %>%
  rename(Trophic_type = Trophic_num)

trait_cols <- c("log_Biovolume", "log_Speed_max", "Niche_breadth", "Marginality", 
                "Spincule", "Colony", "Cyst", "Toxin", "Trophic_type")
trait_labs <- c("log\u2081\u2080Biovolume (\u03BCm)", "log\u2081\u2080Speed_max (\u03BCm s\u207B\u00B9)", "Niche breadth", "Marginality",
                "Spincule (0,1)", "Colony (0,1)", "Cyst (0,1)", "Toxin(0,1)", "Trophic type(0-4)")

mat <- cont_df %>%
  dplyr::select(all_of(trait_cols)) %>%
  rename_with(~ trait_labs) %>%
  as.matrix()

# Spearman's correlation
cor_res  <- rcorr(mat, type = "spearman")
rho      <- cor_res$r
pval     <- cor_res$P

diag(rho) <- NA

# Significant label
sig_label <- function(p) {case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )}
sig_mat <- matrix(sig_label(as.vector(pval)),
                  nrow = nrow(pval),
                  dimnames = dimnames(pval))
diag(sig_mat) <- ""

htmap <- pheatmap(rho, display_numbers = matrix(paste0(round(rho, 2), sig_mat), nrow = nrow(rho)),
                  number_color = "black", fontsize_number = 12,
                  color = colorRampPalette(c("#2166ac","white","#d6604d"))(100),
                  breaks = seq(-1, 1, length.out = 101),
                  clustering_distance_rows = as.dist(1 - abs(rho)),
                  clustering_distance_cols = as.dist(1 - abs(rho)),
                  clustering_method = "average",
                  border_color = "black",
                  na_col = "#d0d0d0",
                  treeheight_row = 40,
                  treeheight_col = 40,
                  fontsize_row = 12,
                  fontsize_col = 12,
                  angle_col = 90
                  )

ggsave(file.path(OUT_DIR, "ft_cor_heatmap.png"),
       htmap, width = 8, height = 8, dpi = 300)

# 2) continuous x binary traits: Mann-Whitney U
mw_res <- lapply(CONT_TRAITS, function(cont) {
  lapply(BIN_TRAITS, function(bin) {
    sub <- ft_df %>% dplyr::select(all_of(cont), all_of(bin)) %>% drop_na()
    g   <- split(sub[[cont]], sub[[bin]])
    if (length(g) < 2) return(NULL)
    mw  <- wilcox.test(g[["0"]], g[["1"]], exact = FALSE)
    data.frame(Cont = cont, Bin = bin, MW_p = mw$p.value)
  }) %>% bind_rows()
}) %>% bind_rows()

# 3) continuous x category traits: Kruskal-Wallis
kw_res <- lapply(CONT_TRAITS, function(cont) {
  sub <- ft_df %>% dplyr::select(all_of(cont), all_of(CAT_TRAITS)) %>% drop_na()
  kw  <- kruskal.test(as.formula(paste(cont, "~", CAT_TRAITS)), data = sub)
  data.frame(Cont = cont, Cat = "Trophic_type", KW_p = kw$p.value)
}) %>% bind_rows()

# 4) binary x category traits: Fisher's exact test
fisher_res <- lapply(BIN_TRAITS, function(bin) {
  sub  <- ft_df %>% dplyr::select(all_of(CAT_TRAITS), all_of(bin)) %>% drop_na()
  tbl  <- table(sub[[CAT_TRAITS]], sub[[bin]])
  ft   <- fisher.test(tbl, simulate.p.value = TRUE)
  data.frame(Bin = bin, Cat = CAT_TRAITS, Fisher_p = ft$p.value)
}) %>% bind_rows()

# ── 5. Save results ──────────────────────────────────────────
cor_df <- cor_res$r %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to  = "Var2", values_to = "cor") %>%
  mutate(n = as.vector(cor_res$n), p = as.vector(cor_res$P)) %>%
  filter(Var1 != Var2) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Var1, Var2)), collapse = "__")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

wb <- createWorkbook()
addWorksheet(wb, "Cont_cor")
addWorksheet(wb, "MW_cont&bin")
addWorksheet(wb, "KW_cont&cat")
addWorksheet(wb, "Fisher_bin&cat")
writeData(wb, "Cont_cor", cor_df)
writeData(wb, "MW_cont&bin", mw_res)
writeData(wb, "KW_cont&cat", kw_res)
writeData(wb, "Fisher_bin&cat", fisher_res)
saveWorkbook(wb, file.path(OUT_DIR, "ft_cors.xlsx"), overwrite = TRUE)

cat("\nDone. All outputs saved to:\n", OUT_DIR, "\n")

