# Challenges of Incorporating Rotation Information into Crop Insurance Rates

**Lawson Connor, Victor Funes-Leal, Eunchun Park**  
Department of Agricultural Economics and Agribusiness  
University of Arkansas

---

## Overview

This repository contains the replication code and manuscript files for *Challenges of Incorporating Rotation Information into Crop Insurance Rates*. The paper estimates how six-year corn-soybean rotation sequences affect the mean and variance of corn and soybean yields in Illinois from 2009 to 2022, using the Just-Pope (1979) production risk framework with a three-stage FGLS estimator. Results are used to develop a rotation score that summarizes actuarially relevant rotation history in a form suitable for crop insurance rating.

The companion paper, *Crop Rotations, Transition Costs, and the Path Dependence of Farm Delinquency*, uses the yield effect estimates from this paper as inputs to a simulation of field-level operating credit dynamics.

---

## Repository structure

```
crop_rotation_and_losses/
├── code/
│   ├── rotation_setup.R        # Shared setup: packages, helpers, PCA, bootstrap function
│   ├── corn_analysis.R         # Full corn pipeline: OLS, RCI, VPD, JP, FGLS, maps
│   └── soy_analysis.R          # Full soy pipeline: OLS, RCI, VPD, JP, FGLS, maps
├── figures/                    # All figures saved by the analysis scripts (PNG, 300 dpi)
├── tables/                     # All tables saved by the analysis scripts (LaTeX)
├── old/                        # Archived prior versions of analysis scripts
├── Makefile                    # Builds the paper PDF from .tex sources
├── rotations_losses.tex        # Main manuscript
├── abstract.tex                # Abstract (included in main)
├── dynamic_optimization.tex    # Appendix: dynamic optimization framework
├── bibliography.bib            # BibTeX references
└── mplainnat.bst               # Bibliography style file
```

---

## Data

Data are not included in this repository. The analysis uses the following sources, all of which require separate access or download:

| Dataset | Source | Variables used |
|---|---|---|
| QDANN yield estimates | Ma et al. (2024) | `corn_yield`, `soy_yield` |
| SCYM yield estimates | Lobell et al. (2015) | `corn_yield`, `soy_yield` |
| Cropland Data Layer (CDL) | USDA-NASS | `rot_crop`, `RCI` |
| GRIDMET | Abatzoglou (2013) | `pr_*`, `tmmx_*`, `tmmn_*`, `vpd_*`, `soil_*`, `cGDD_*` |
| TerraClimate | Abatzoglou et al. (2018) | `pdsi_mean` |
| gSSURGO | USDA-NRCS | `nccpi3all_mean`, `rootznaws_mean`, `soc0_100_mean` |

The processed analysis files are stored at:
```
corn_rci_il_long.csv   # Corn field-year panel with RCI and all covariates
soy_rci_il_long.csv    # Soybean field-year panel with RCI and all covariates
```

Each row is a field-tile × year observation. Both files share the same non-yield columns. The corn file contains `corn_yield`; the soy file contains `soy_yield`.

---

## Running the analysis

The analysis is split into two independent scripts that share a common setup file. Run them in a fresh R session in order:

```r
# Step 1 — Corn analysis (produces all corn tables, figures, and bootstrap)
source("code/corn_analysis.R")

# Step 2 — Soy analysis (produces all soy tables, figures, and bootstrap)
# Run in a fresh R session to avoid memory conflicts with corn objects
source("code/soy_analysis.R")
```

Both scripts source `rotation_setup.R` automatically. Set the working directory and file paths in `rotation_setup.R` before running.

**Runtime**: each script takes approximately 2–4 hours on a modern workstation, dominated by the FGLS bootstrap (B = 499 iterations on ~1.5M observations). Bootstrap results are saved incrementally every 100 iterations to `boot_progress_N.rds` in case of interruption, and the final result is saved to `boot_corn.rds` / `boot_soy.rds`.

---

## R packages required

```r
install.packages(c(
  "tidyverse", "data.table", "fixest", "broom",
  "marginaleffects", "hdm", "dotwhisker",
  "ggfortify", "ggrepel", "sf", "usmap", "statar"
))
```

R version 4.2.0 or later is required. The analysis was developed and tested on R 4.6.0 (Windows).

---

## Tables and figures produced

### Tables (saved to `tables/` as `.tex`)

| File | Label | Content |
|---|---|---|
| `corn_rot.tex` | `tab:corn_rot` | Rotation patterns and corn yields (OLS) |
| `soy_rot.tex` | `tab:soy_rot` | Rotation patterns and soybean yields (OLS) |
| `corn_rci.tex` | `tab:corn_rci` | Rotational Complexity Index and corn yields |
| `soy_rci.tex` | `tab:soy_rci` | Rotational Complexity Index and soybean yields |
| `corn_rot_vpd.tex` | `tab:corn_rot_vpd` | Rotation × VPD interaction, corn |
| `soy_rot_vpd.tex` | `tab:soy_rot_vpd` | Rotation × VPD interaction, soy |
| `corn_rci_vpd.tex` | `tab:corn_rci_vpd` | RCI × VPD interaction, corn |
| `soy_rci_vpd.tex` | `tab:soy_rci_vpd` | RCI × VPD interaction, soy |
| `corn_jp_mean.tex` | `tab:corn_jp_mean` | Just-Pope stage 1: corn yield mean |
| `corn_jp_var.tex` | `tab:corn_jp_var` | Just-Pope stage 2: corn yield variance (OLS) |
| `soy_jp_mean.tex` | `tab:soy_jp_mean` | Just-Pope stage 1: soybean yield mean |
| `soy_jp_var.tex` | `tab:soy_jp_var` | Just-Pope stage 2: soybean yield variance (OLS) |
| `corn_rci_jp.tex` | `tab:corn_rci_jp` | Just-Pope factor RCI: corn yield moments |
| `soy_rci_jp.tex` | `tab:soy_rci_jp` | Just-Pope factor RCI: soybean yield moments |

### Figures (saved to `figures/` as `.png`, 300 dpi)

| File | Content |
|---|---|
| `rot_pca_plot.png` | PCA biplot of corn-soy rotation features (Section 4) |
| `corn_rot_plot.png` | Response of corn yields to rotation sequences |
| `soy_rot_plot.png` | Response of soybean yields to rotation sequences |
| `corn_rci_plot.png` | Response of corn yields to RCI levels |
| `soy_rci_plot.png` | Response of soybean yields to RCI levels |
| `corn_var_plot.png` | Response of corn yield std dev to rotation sequences |
| `soy_var_plot.png` | Response of soybean yield std dev to rotation sequences |
| `corn_coeff_plot.png` | Corn mean vs variance coefficient scatter |
| `soy_coeff_plot.png` | Soybean mean vs variance coefficient scatter |
| `corn_jp_plot.png` | Just-Pope mean-variance decomposition (corn) |
| `soy_jp_plot.png` | Just-Pope mean-variance decomposition (soy) |
| `rci_plot.png` | Nonlinear effect of RCI on corn yield mean and variance |
| `corn_yield_map.png` | Spatial map of corn yields, Illinois 2016 |
| `soy_yield_map.png` | Spatial map of soybean yields, Illinois 2016 |
| `rci_map.png` | Spatial map of RCI values, Illinois 2016 |
| `nccpi_corn_map.png` | Spatial map of NCCPI corn productivity index |
| `nccpi_soy_map.png` | Spatial map of NCCPI soybean productivity index |

---

## Building the paper

```bash
make
```

The `Makefile` compiles `rotations_losses.tex` using `pdflatex` and `bibtex`. Tables are included via `\input{}` from the `tables/` directory; figures are included via `\includegraphics{}` from the `figures/` directory.

---

## Econometric framework

The paper implements the Just-Pope (1979) stochastic production function in a two-way fixed effects panel setting:

**Stage 1 (conditional mean)**:
$$y_{it} = \text{seq}_{it}'\tau + W_{it}'\beta + \lambda_k + \gamma_t + \varepsilon_{it}$$

**Stage 2 (conditional variance)**:
$$\hat{u}_{it}^2 = \text{seq}_{it}'\pi + W_{it}'\delta + \lambda_k + \gamma_t + \nu_{it}$$

where $\lambda_k$ are county (FIPS) fixed effects and $\gamma_t$ are year fixed effects. The variance stage uses three-stage FGLS following Harvey (1976) and Saha, Havenner, and Talpaz (1997). Inference uses a pairs cluster bootstrap at the field level (B = 499) wrapping all three stages to account for the generated-regressor problem.

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
