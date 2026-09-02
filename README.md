# Challenges of Incorporating Rotation Information into Crop Insurance Rates

**Lawson Connor, Victor Funes-Leal, Eunchun Park**  
Department of Agricultural Economics and Agribusiness  
University of Arkansas

---

## Overview

This repository contains the replication code and manuscript files for *Challenges of Incorporating Rotation Information into Crop Insurance Rates*. The paper estimates how six-year crop rotation sequences affect the mean, variance, and skewness of **corn** yields in Illinois from 2003 to 2022, using the Just-Pope (1979) production risk framework with a three-stage estimator (FGLS variance stage) and a field-level pairs cluster bootstrap. Rotation histories are built from corn, soybean, and winter-wheat years; corn is the outcome crop. Results are used to develop a rotation score that summarizes actuarially relevant rotation history in a form suitable for crop insurance rating.

**Manuscript of record: `just_pope_rotations.tex`.** Earlier drafts (`rotations_losses.tex`, `rotations_insurance.tex`) are archived under `old/`.

The companion paper, *Crop Rotations, Transition Costs, and the Path Dependence of Farm Delinquency*, uses the yield effect estimates from this paper as inputs to a simulation of field-level operating credit dynamics.

---

## Repository structure

```
crop_rotation_and_losses/
├── code/
│   ├── corn_analysis_full.R          # Full corn pipeline: OLS, LASSO, RCI, regimes, VPD, JP moments, maps
│   ├── rotation_setup_wa.R           # Shared setup: packages, helpers, sequence set, PCA, bootstrap fn
│   ├── just_pope_bootstrap_moments.R # Three-stage B=999 bootstrap (mean/variance/skewness), sourced by the pipeline
│   ├── jp_boot_vcov.R                # Turns bootstrap draws into vcov matrices for etable()
│   ├── zvector_bootstrap_var.R       # B=999 bootstrap for the parsimonious (Z-vector) variance stage
│   ├── tables_combined.r             # Z-vector (structural-feature) table
│   ├── lasso_rotation_selection.R    # LASSO sequence selection (tab:corn_lasso)
│   ├── slx_corn_regime.r             # Rotation-regime / SLX models (tab:slx_corn_regime*)
│   ├── rci_vectorized.R, cdl_recode.R, rotation_features_multicrop.R  # helpers
│   ├── main_analysis.R               # Older corn+soy, corn-soybean-only pipeline (superseded)
│   ├── soybeans/                     # Soybean-outcome pipeline (not used by the current manuscript)
│   └── old/                          # Archived prior versions of analysis scripts
├── figures/                          # All figures saved by the analysis scripts (PNG, 300 dpi)
├── tables/                           # All tables saved by the analysis scripts (LaTeX)
├── old/                              # Archived prior manuscripts and presentations
├── Makefile                          # Builds the paper PDF from just_pope_rotations.tex
├── just_pope_rotations.tex           # Main manuscript (manuscript of record)
├── abstract.tex                      # Abstract (included in main)
├── dynamic_optimization.tex          # Appendix: dynamic optimization framework
├── bibliography.bib                  # BibTeX references
└── mplainnat.bst                     # Bibliography style file (shipped; manuscript currently sets apalike)
```

---

## Data

Data are not included in this repository. The analysis uses the following sources, all of which require separate access or download:

| Dataset | Source | Variables used |
|---|---|---|
| QDANN yield estimates | Ma et al. (2024) | `corn_yield` |
| Cropland Data Layer (CDL) | USDA-NASS | `rot_crop`, `RCI` |
| GRIDMET | Abatzoglou (2013) | `pr_*`, `tmmx_*`, `tmmn_*`, `vpd_*`, `vpdmax_*`, `soil_*`, `cGDD_*` |
| TerraClimate | Abatzoglou et al. (2018) | climate water deficit |
| gSSURGO | USDA-NRCS | `nccpi3*_mean`, `rootznaws_mean`, `soc0_100_mean` |

The current pipeline reads a single processed field-year panel (Parquet):

```
d_igis13_12_1_2025.with_rci.parquet   # Illinois corn field-year panel with rotation history, RCI, weather, soil
```

Each row is a field-tile × year observation. The path is set at the top of `code/corn_analysis_full.R`.

---

## Running the analysis

Run in a fresh R session:

```r
setwd("code")
source("corn_analysis_full.R")   # sources rotation_setup_wa.R, then just_pope_bootstrap_moments.R
```

Edit the working directory and the Parquet path near the top of `corn_analysis_full.R` before running. The script produces every corn table and figure used in `just_pope_rotations.tex`, then the spatial maps.

`code/zvector_bootstrap_var.R` can be sourced after the Z-vector models are built to attach bootstrap SEs to the parsimonious variance stage.

**Runtime**: several hours on a modern workstation, dominated by the three-stage bootstrap (B = 999 on ~1.79M observations). `just_pope_bootstrap_moments.R` writes `boot_progress_N.rds` every 100 iterations in case of interruption and saves the final result to `boot_moments.rds`.

---

## R packages required

```r
install.packages(c(
  "tidyverse", "data.table", "fixest", "broom", "arrow",
  "marginaleffects", "hdm", "dotwhisker",
  "ggfortify", "ggrepel", "patchwork", "knitr", "kableExtra",
  "sf", "usmap", "statar"
))
```

R version 4.2.0 or later is required. The analysis was developed and tested on R 4.6.0 (Windows).

---

## Tables and figures produced

### Tables (saved to `tables/` as `.tex`, `\input` by the manuscript)

| File | Label | Content |
|---|---|---|
| `corn_summary_stats.tex` | `tab:corn_summary` | Summary statistics by rotation type |
| `corn_rot.tex` | `tab:corn_rot` | Rotation sequences and corn yields (28 non-monoculture sequences) |
| `corn_lasso.tex` | `tab:corn_lasso` | LASSO-selected rotation sequences (post-selection OLS) |
| `zvector.tex` | `tab:zvector` | Structural-feature (Z-vector) mean and variance effects |
| `corn_rci.tex` | `tab:corn_rci` | Factor RCI and corn yields |
| `slx_corn_regime.tex` | `tab:slx_corn_regime` | Rotation-regime model |
| `slx_corn_regime_gaps.tex` | `tab:slx_corn_regime_gaps` | Regime model, gap terms |
| `slx_corn_regime_ctrl.tex` | `tab:slx_corn_regime_ctrl` | Regime model with soil controls |
| `slx_corn_regime_hybrid.tex` | `tab:slx_corn_regime_hybrid` | Regime × RCI hybrid |
| `corn_rot_vpd.tex` | `tab:corn_rot_vpd` | Rotation sequences × growing-season VPD/drought |
| `corn_rci_vpd.tex` | `tab:corn_rci_vpd` | RCI × July drought interaction |
| `corn_jp_moments.tex` | `tab:corn_jp_moments` | Just-Pope stage 2/3: variance (FGLS) and standardized skewness, bootstrap SEs |
| `corn_rci_jp.tex` | `tab:corn_rci_jp` | Just-Pope factor RCI: corn yield moments |

### Figures (saved to `figures/` as `.png`, 300 dpi)

| File | Content |
|---|---|
| `corn_rot_plot.png` | Response of corn yields to rotation sequences |
| `corn_var_plot.png` | Response of corn yield std dev to rotation sequences |
| `corn_jp_plot.png` | Just-Pope mean-variance decomposition (corn) |
| `corn_rci_plot.png` | Nonlinear effect of RCI on corn yield mean and variance |
| `score_yield.png` | Predicted corn yield as a function of the rotation score |
| `corn_yield_map.png` | Spatial map of corn yields, Illinois 2016 |
| `rci_map.png` | Spatial map of RCI values, Illinois 2016 |
| `nccpi_corn_map.png` | Spatial map of NCCPI corn productivity index |

---

## Building the paper

```bash
make
```

The `Makefile` compiles `just_pope_rotations.tex` with `pdflatex` and `bibtex`. Tables are pulled in via `\input{}` from `tables/`; figures via `\includegraphics{}` from `figures/`. The manuscript sets `\bibliographystyle{apalike}` and `\bibliography{bibliography.bib}`.

---

## Econometric framework

The paper implements the Just-Pope (1979) stochastic production function in a two-way fixed effects panel setting.

**Stage 1 (conditional mean)**:
$$y_{it} = \text{seq}_{it}'\tau + W_{it}'\beta + \lambda_i + \gamma_t + \varepsilon_{it}$$

**Stage 2 (conditional variance)**:
$$\hat{u}_{it}^2 = \text{seq}_{it}'\pi + W_{it}'\delta + \lambda_i + \gamma_t + \nu_{it}$$

**Stage 3 (conditional skewness)** regresses $\hat{u}_{it}^3$ on the same specification (Antle 1983).

Here $\lambda_i$ are field (`tile_field_ID`) fixed effects and $\gamma_t$ are year fixed effects, so identification of $\hat{\tau}$ comes from within-field variation in rotation sequences over time. The variance stage uses FGLS (Harvey 1976). Analytic standard errors are two-way clustered at the county (`COUNTY_FIPS`) and year levels. Because the stage-2 and stage-3 dependent variables are generated regressors, inference of record for those stages is a pairs cluster bootstrap at the field level (B = 999) that re-estimates all three stages inside each replicate.

---

## Citation

If you use this code or data, please cite:

```
Connor, L., V. Funes-Leal, and E. Park. 2025. "Challenges of Incorporating
Rotation Information into Crop Insurance Rates." Working paper, University
of Arkansas Department of Agricultural Economics and Agribusiness.
```

---

## License

Code is released under the MIT License. See `LICENSE` for details.

---

## Contact

Victor Funes-Leal — vf006@uark.edu  
Department of Agricultural Economics and Agribusiness  
University of Arkansas, Fayetteville AR 72701
