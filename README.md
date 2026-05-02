# ACH    
# Abundant-Centre Hypothesis — Jangcheon Harbor Dinoflagellates

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> Tests the **Abundant-Centre Hypothesis (ACH)** for 32 dinoflagellate taxa using a 411-day fixed-station monitoring dataset from Jangcheon Harbor, South Korea (2020–2021).  
> Compares **three niche centre definitions × three distance metrics** (9 settings) and evaluates how functional traits mediate ACH support patterns.


## Overview

The ACH predicts that species abundance peaks near the centre of the environmental niche and declines towards the periphery. This pipeline:

1. Reduces 11 environmental variables to a 2D PCA space (`dudi.pca`)
2. Builds per-species niche envelopes using three centre definitions
3. Computes distance from each observation to each centre
4. Tests the distance–abundance relationship via Spearman's ρ + BH-FDR
5. Compares which centre definition best supports ACH


## Centre definitions compared

| Centre | Definition | ID |
|--------|-----------|-----|
| **CH** | Centroid of Convex Hull vertices in PCA space | E1–E3 |
| **MVE** | Centre of Minimum Volume Ellipsoid (MASS::cov.mve) | E4–E6 |
| **DMF** | Centroid of PCA coordinates on Day(s) of Maximum Fitness | E7–E9 |

**DMF centre** is derived from phenology analysis output (`5_phenology_variables.xlsx`).  
When a species has a single DMF event, its environmental position on that day is used directly.  
When multiple events exist, the centroid of all DMF environmental positions is computed.


## Distance metrics

| Metric | Description | Settings |
|--------|-------------|----------|
| **Euclidean** | Straight-line distance to centre | E1, E4, E7 |
| **Mahalanobis** | Covariance-weighted distance to centre | E2, E5, E8 |
| **Margin** | Interior distance to envelope boundary | E3, E6, E9 |

For margin distances, Spearman ρ sign is reversed so that negative ρ always indicates ACH support (abundance decreases toward the periphery).

> **E9 (DMF-Margin)** requires ≥ 3 DMF events to form a convex hull; otherwise set to `NA`.


## Pipeline

```
Environmental PCA (dudi.pca, 11 variables)
         │
         ▼
make_envelope()  ────────────────────────────────────────────────
  Per-species:                                                   │
    CH centre  (chull, centroid of vertices)                     │
    MVE centre (MASS::cov.mve)                                   │
    DMF centre (from phenology output, date-matched to PCA)      │
         │
         ▼
calc_all_dist()  →  9 distance vectors per species (E1–E9)
         │
         ▼
test_sp_ach()    →  Spearman ρ + p-value per species × setting
         │
         ▼
BH-FDR           →  classify: Supported / Opposite / Not significant
         │
         ▼
Outputs: ACH heatmap, violin, proportion plots, GAM curves,
         centre comparison summary (Fig 5)
```


## Repository structure

```
ach-dinoflagellate/
│
├── run_ACH.R                       # Main analysis script
│
├── R/                              # Helper functions (sourced by run_ACH.R)
│   ├── utils.R                     # std01, bh_adjust, classify_ach, get_sig_label
│   ├── geometry.R                  # eucl_dist, maha_dist, p2seg_dist,
│   │                               #   margin_nr_dist, margin_mve_dist, ellipse_df
│   ├── make_envelope.R             # Per-species CH / MVE / DMF envelope builder
│   ├── calc_all_dist.R             # 9-setting distance calculator
│   ├── test_sp_ach.R               # Spearman ACH test
│   ├── build_model_df.R            # Long-format model data
│   ├── plot_sp_niche.R             # Species niche panel plot
│   └── fit_quad_lm_sp.R            # Species-specific quadratic LM curves
│
├── data/                           # Input data (not tracked by Git — see .gitignore)
│   ├── JC_envs_daily.xlsx
│   ├── JC_abundance.xlsx
│   ├── dino_5day.xlsx
│   └── 5_phenology_variables.xlsx  # Output from phenology-dinoflagellate pipeline
│
├── Output_ACH_YYMMDD/              # Auto-created output directory
│   ├── dudi_loadings.csv
│   ├── OMI_niche_parameters.csv
│   ├── DMF_centre_coverage.csv
│   ├── ACH_sp_results.csv
│   ├── ACH_centre_comparison.csv   ← 3-centre comparison summary
│   ├── ACH_model_data_long.csv
│   ├── Fig1_species_niches_all_panels.pdf
│   ├── Fig2_1_Spearman_violin.pdf
│   ├── Fig2_2_ACH_heatmap.pdf
│   ├── Fig3_ACH_support_proportion.pdf
│   ├── Fig4_predictive_models/
│   └── Fig5_centre_comparison/     ← new
│
├── .gitignore
└── README.md
```


## Requirements

```r
pkgs <- c(
  "openxlsx", "dplyr", "tidyr", "tibble", "purrr", "stringr", "lubridate",
  "ggplot2", "patchwork", "scales",
  "ade4", "MASS", "sp", "ellipse", "mgcv"
)
```

Optional (Bayesian models):
```r
bayes_pkgs <- c("brms", "posterior")
```

All required packages are installed automatically on first run.  
Tested on **R ≥ 4.2** (macOS / Linux).


## Usage

1. Clone the repository.
2. Place input files in `data/` (see structure above).
3. Open `run_ACH.R` and confirm `INPUT_DIR <- "."` points to the project root.
4. Source the script:

```r
source("run_ACH.R")
```


## Bug fixes from original version

| ID | Severity | Description |
|----|----------|-------------|
| B1 | Critical | `dino_sm7` undefined variable → `dino_sm5` |
| B2 | Critical | `make_envelope` signature / return value completely redesigned to match caller |
| B3 | Critical | `calc_all_dist` function name and argument mismatch resolved |
| B4 | Critical | `fit_quad_lm_species` → correct name `fit_quad_lm_sp` |
| B5 | Moderate | `arrange(., by = "Date")` → `arrange(Date)` |
| L1 | Logic    | CH margin now uses `p2seg_dist()` (edge distance) not vertex distance |
| L2 | Logic    | 2D convex hull uses `chull()` instead of `convhulln()` (eliminates centroid bias) |
| L3 | Logic    | Removed `library()` calls inside function bodies |
| L4 | Logic    | `if (as.numeric(level))` dead branch in `ellipse_df` removed |
| L5 | Logic    | `build_model_df` now takes `min_occ` as explicit argument |
| L7 | Logic    | Hardcoded column indices (`4:37`) replaced with name-based selection |
| L8 | Logic    | `DOY` coerced to `integer` (was `character`) |
| W1 | Warning  | Hardcoded absolute paths replaced with relative paths |


## New feature: DMF centre comparison (E7–E9)

The DMF (Day of Maximum Fitness) centre is derived from the phenology pipeline:

```
5_phenology_variables.xlsx
    └─ DMF_Date (per species, per bloom event)
           ↓ matched to dudi$li site scores
    DMF centre = colMeans(PCA scores on all DMF dates)
```

The comparison is summarised in `ACH_centre_comparison.csv` and visualised in `Fig5_centre_comparison/`:

- **Fig 5-1**: Median Spearman ρ by centre × distance metric (bar chart + Wilcoxon label)
- **Fig 5-2**: Pairwise CH vs. DMF ρ scatter (colour = MVE ρ, Euclidean only)


## Links

- Phenology pipeline: [phenology-dinoflagellate](../phenology-dinoflagellate)
- Functional trait database: see manuscript supplementary materials


## Citation

> [Author(s)]. (*in prep.*). Functional traits mediate the Abundant-Centre Hypothesis in coastal dinoflagellate communities. *Ecology Letters* / *Limnology and Oceanography*.


## License

MIT © [Author Name] — see [LICENSE](LICENSE) for details.
