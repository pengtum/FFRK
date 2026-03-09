# FFRK: Focal-Feature Regression Kriging

This repository contains the implementation of **FFRK (Focal-Feature Regression Kriging)**, baseline/comparison methods, experiment data, paper figures, and reproducible results.

## 1. Repository Structure

```text
.
├── Code/
│   ├── FFRK/                         # Core FFRK R package
│   ├── Fig6_Evaluation_14models/     # 14-model benchmark (Fig. 6)
│   ├── FFRK_Ablation/                # FFRK ablation experiments
│   ├── FFRKvsFFR/                    # FFRK vs FFR comparison
│   ├── Fig7&8_Sensitivity_KQ/        # Sensitivity analysis of k and q (Fig. 7-8)
│   ├── Tab2_Evaluiation_HIER/        # Table 2 (HIER) evaluation
│   ├── Tab2_Evaluiation_STK/         # Table 2 (STK) evaluation
│   ├── FigS1-S3_Scatter/             # Supplementary scatter plots (S1-S3)
│   ├── Figx_Variograms/              # Variogram analysis
│   ├── Figx_Variance/                # Kriging variance analysis and visualization
│   ├── Fig5_Global/                  # Global prediction visualization (Fig. 5)
│   └── VIF_P-value_Feature Importance/
├── Data/                             # Input datasets
├── Results/                          # Generated outputs (csv/jpg/txt, etc.)
├── Figures/                          # Final paper figures
└── readme.md
```

## 2. Requirements

- R >= 4.1
- Recommended OS: macOS / Linux / Windows

Main R packages used in this project:

`devtools`, `this.path`, `caret`, `sp`, `gstat`, `automap`, `geosimilarity`, `doParallel`, `progress`, `rpart`, `randomForest`, `dplyr`, `ggplot2`, `tidyr`, `readxl`, `patchwork`, `ggpubr`, `gridExtra`, `car`

You can install most dependencies with:

```r
install.packages(c(
  "devtools","this.path","caret","sp","gstat","automap","geosimilarity",
  "doParallel","progress","rpart","randomForest","dplyr","ggplot2","tidyr",
  "readxl","patchwork","ggpubr","gridExtra","car","readr","MASS","scales"
))
```

## 3. Quick Start

From the repository root:

```bash
Rscript Code/Fig6_Evaluation_14models/Fig6_Evaluation_14models.R
```

This script calls core functions in `Code/FFRK/` and writes benchmark outputs to:

`Results/Fig6_Evaluation_14models/`

## 4. Main Reproducible Experiment Entrypoints

- `Rscript Code/Fig6_Evaluation_14models/Fig6_Evaluation_14models.R`  
  Output: `Results/Fig6_Evaluation_14models/`
- `Rscript Code/FFRK_Ablation/FFRK_Ablation.R`  
  Output: `Results/FFRK_Ablation/`
- `Rscript Code/FFRKvsFFR/FFRKvsFFR.R`  
  Output: `Results/FFRKvsFFR/FFRKvsFFR.csv`
- `Rscript "Code/Fig7&8_Sensitivity_KQ/Fig7&8_Sensitivity_KQ.R"`  
  Output: `Results/Fig7&8_Sensitivity_KQ/`
- `Rscript Code/Tab2_Evaluiation_HIER/Tab2_Evaluiation_HIER.R`  
  Output: `Results/Tab2_Evaluiation_HIER/`
- `Rscript Code/Tab2_Evaluiation_STK/Tab2_Evaluiation_STK.R`  
  Output: `Results/Tab2_Evaluiation_STK/`
- `Rscript Code/FigS1-S3_Scatter/FigS1-S3_Scatter.R`  
  Output: `Results/Scatter/`

## 5. Data Description

- `Data/data_Cu.csv`, `Data/data_Pb.csv`, `Data/data_Zn.csv`: modeling inputs for Cu, Pb, and Zn.
- Main fields include: target variable (e.g., `Cu_ppm`), coordinates (`DLONG`, `DLAT`), and environmental covariates (`Dlith`, `Dfault`, `Slope`, `Water`, `NDVI`, `MainRd`, `Road`, `SOC`, `pH`).
- `Data/Global/`: inputs related to global mapping and pre-trained model objects.

## 6. Where FFRK Is Implemented

Core method code is in `Code/FFRK/R/`:

- `feature_extract.R`: geo-self feature extraction (IDW, quantile features, similarity feature)
- `modeling.R`: training/evaluation workflows for OK/UK/ML/RK/FFRK/HIER
- `prepare_data.R`: data preprocessing (including optional log-transform)
- `metrics.R`: R2/MAE/RMSE metrics

## 7. Notes

- `Results/` and `Figures/` already include most outputs corresponding to the manuscript.
- For reproducibility, keep random seeds fixed (the scripts use `20250718`) and use consistent package versions.

## 8. Citation

If this repository is useful for your research, please cite the corresponding FFRK paper: Luo, P., Wu, Y., Song, Y., Focal-Feature Regression Kriging. Geographical Analysis
