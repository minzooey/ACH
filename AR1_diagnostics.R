## ── AR(1) diagnostics: identifiability check for autocorrelation-corrected
##    hierarchical models ──────────────────────────────────────────────────
##
##  Purpose
##    Fig2_final.R 결과에서 AR(1) 모형의 species-level β1 가 Base 모형 대비
##    전반적으로 0으로 수축되는 현상이 관찰됨. 이것이 (a) 정당한 자기상관
##    보정에 따른 attenuation인지, 아니면 (b) dist_signed의 자체 시간적
##    자기상관과 AR(1) 잔차구조가 같은 시간축을 공유하면서 발생하는
##    비식별(non-identifiability)/구조 붕괴인지 판별하기 위한 진단.
##
##  Checks performed (per setting, Base vs AR1)
##    1. AR(1) coefficient posterior (median, 95% CI, Rhat, ESS)
##       → 0.8~0.9 이상이면 거의 unit-root 급으로 추세를 흡수하고 있다는 신호
##    2. Random-effect SD posteriors: sd(species: Intercept, dist_signed,
##       I(dist_signed^2)) — AR1에서 이 값들이 Base 대비 0 근처로 붕괴했는지
##    3. Rhat / bulk-ESS / tail-ESS for b_dist_signed, b_Idist_signedE2
##       → 수렴은 됐지만 정보량이 사실상 없는 상태(낮은 ESS, prior 근접)인지
##    4. Divergent transitions (per chain, via brms::nuts_params)
##
##  Prerequisites (must exist in environment before sourcing this file):
##    bayes_fits      — named list of brms fits (base quadratic model)
##    bayes_fits_ar1  — named list of brms fits (AR(1) residual model)
##    OUT_DIR         — output directory path
##    FIG2_SETTING    — setting to inspect in detail (e.g. "E5"); diagnostics
##                       for ALL settings shared by both lists are also computed
##
##  Required packages: dplyr, purrr, tibble, posterior, brms, ggplot2
## ─────────────────────────────────────────────────────────────────────────────

pkgs <- c("dplyr", "purrr", "tibble", "posterior", "brms", "ggplot2")
invisible(lapply(pkgs, library, character.only = TRUE))

stopifnot(
  exists("bayes_fits"),
  exists("bayes_fits_ar1"),
  exists("OUT_DIR"),
  exists("FIG2_SETTING")
)

# ── 1. Per-fit diagnostic extractor ──────────────────────────────────────────
## Pulls: AR(1) coefficient (if present), random-effect SDs for species,
## and fixed-effect Rhat/ESS — all in one summarise_draws() call so the
## Rhat/ESS numbers are computed identically for both model types.
extract_ar1_diag <- function(fit, model_label, setting) {

  vars <- posterior::variables(fit)

  # AR term naming has varied slightly across brms versions ("ar[1]", "ar1");
  # match anything starting with "ar" that is NOT "ar_diag"/etc. to be safe.
  ar_vars <- grep("^ar(\\[|$|1)", vars, value = TRUE)

  # Random-effect SDs for the species grouping factor
  sd_vars <- grep("^sd_species__", vars, value = TRUE)

  # Fixed effects of interest
  fx_vars <- intersect(c("b_dist_signed", "b_Idist_signedE2"), vars)

  target_vars <- c(ar_vars, sd_vars, fx_vars)
  if (length(target_vars) == 0) {
    warning("[", model_label, " / ", setting, "] No matching variables found.")
    return(NULL)
  }

  draws <- posterior::as_draws_df(fit)

  diag_tab <- posterior::summarise_draws(
    draws[, target_vars, drop = FALSE],
    median,
    ~ quantile2(.x, probs = c(0.025, 0.975)),
    rhat, ess_bulk, ess_tail
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(parameter = variable) %>%
    dplyr::mutate(
      model   = model_label,
      setting = setting,
      # Flag: parameter type for downstream filtering/plotting
      param_type = dplyr::case_when(
        parameter %in% ar_vars ~ "AR1_coef",
        parameter %in% sd_vars ~ "RandomEffect_SD",
        parameter %in% fx_vars ~ "FixedEffect"
      )
    )

  # Divergent transitions (AR1 models often need higher adapt_delta;
  # divergences indicate the sampler struggled with the posterior geometry —
  # relevant context if fixed effects also look collapsed/unstable)
  n_div <- tryCatch({
    np <- brms::nuts_params(fit)
    sum(np$Parameter == "divergent__" & np$Value == 1)
  }, error = function(e) NA_integer_)

  diag_tab$n_divergent <- n_div
  diag_tab
}


# ── 2. Run across all settings shared by both model lists ────────────────────
common_settings <- intersect(names(bayes_fits), names(bayes_fits_ar1))
if (length(common_settings) == 0)
  stop("No settings shared between bayes_fits and bayes_fits_ar1.")

diag_all <- purrr::map_dfr(common_settings, function(st) {
  dplyr::bind_rows(
    extract_ar1_diag(bayes_fits[[st]],     "Base", st),
    extract_ar1_diag(bayes_fits_ar1[[st]], "AR1",  st)
  )
})

# Save full diagnostic table (all settings, both models)
diag_out <- file.path(OUT_DIR, "AR1_diagnostics_allSettings.csv")
write.csv(diag_all, diag_out, row.names = FALSE)
message("✓ AR1 diagnostics (all settings) saved: ", diag_out)


# ── 3. Flag settings/parameters that look problematic ────────────────────────
## Heuristics (adjust thresholds as needed):
##   - AR1 coefficient median > 0.7   → strong absorption risk
##   - RandomEffect_SD median < 0.1   → species-level slope variance collapsed
##   - Rhat > 1.01 or ess_bulk < 400  → convergence/information concerns
##   - n_divergent > 0                → sampler geometry issues
flags <- diag_all %>%
  dplyr::mutate(
    flag_ar1_absorb   = param_type == "AR1_coef"        & median > 0.7,
    flag_sd_collapse  = param_type == "RandomEffect_SD" & median < 0.1,
    flag_convergence  = rhat > 1.01 | ess_bulk < 400,
    flag_divergent    = n_divergent > 0
  ) %>%
  dplyr::filter(flag_ar1_absorb | flag_sd_collapse | flag_convergence | flag_divergent)

cat("\n── Flagged parameters (potential AR(1) identifiability issues) ──\n")
if (nrow(flags) > 0) {
  print(flags %>%
          dplyr::select(setting, model, parameter, param_type, median,
                        rhat, ess_bulk, n_divergent,
                        flag_ar1_absorb, flag_sd_collapse,
                        flag_convergence, flag_divergent),
        n = 100)
} else {
  cat("None flagged under current thresholds.\n")
}


# ── 4. Detailed view for FIG2_SETTING ─────────────────────────────────────────
diag_focus <- diag_all %>% dplyr::filter(setting == FIG2_SETTING)

cat("\n── Diagnostics for setting:", FIG2_SETTING, "──\n")
print(diag_focus %>%
        dplyr::select(model, parameter, param_type, median, q2.5, q97.5,
                      rhat, ess_bulk, ess_tail, n_divergent),
      n = 50)

# Side-by-side comparison of random-effect SD shrinkage (Base vs AR1)
sd_compare <- diag_focus %>%
  dplyr::filter(param_type == "RandomEffect_SD") %>%
  dplyr::select(model, parameter, median, q2.5, q97.5) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(median, q2.5, q97.5))

cat("\n── Random-effect SD: Base vs AR1 (", FIG2_SETTING, ") ──\n")
print(sd_compare)


# ── 5. Posterior density plot: AR(1) coefficient + key random-effect SDs ─────
fit_ar1_focus <- bayes_fits_ar1[[FIG2_SETTING]]
draws_focus   <- posterior::as_draws_df(fit_ar1_focus)
vars_focus    <- posterior::variables(fit_ar1_focus)

plot_vars <- c(
  grep("^ar(\\[|$|1)", vars_focus, value = TRUE),
  grep("^sd_species__", vars_focus, value = TRUE)
)

if (length(plot_vars) > 0) {
  dens_df <- draws_focus %>%
    dplyr::select(dplyr::all_of(plot_vars)) %>%
    tidyr::pivot_longer(dplyr::everything(),
                         names_to = "parameter", values_to = "value")

  p_ar1_diag <- ggplot(dens_df, aes(x = value)) +
    geom_density(fill = "#1A7595", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    facet_wrap(~ parameter, scales = "free", ncol = 2) +
    labs(
      x = "Posterior draw",
      y = "Density",
      title = sprintf("AR(1) model diagnostics — %s", FIG2_SETTING),
      subtitle = "AR(1) coefficient near 1 and/or SD collapsing to 0 indicate absorption of the dist_signed signal"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggplot2::ggsave(
    file.path(OUT_DIR, sprintf("AR1_diagnostic_posteriors_%s.pdf", FIG2_SETTING)),
    p_ar1_diag, width = 8, height = 5.5, dpi = 300
  )
  message("✓ AR1 diagnostic posterior plot saved for ", FIG2_SETTING)
} else {
  message("No AR(1)/random-effect-SD variables found to plot for ", FIG2_SETTING,
          " — check posterior::variables(bayes_fits_ar1[[FIG2_SETTING]]) for exact names.")
}

message("\n✓ AR(1) diagnostics complete. Review flagged parameters above before",
        " deciding which model to report in the main text.")
