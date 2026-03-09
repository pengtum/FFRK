# FFRK Project Overview

This repository organizes all codes, data, figures, and results related to the FFRK project.

## Directory Structure

- `Code/`: Model implementations, experiment scripts, and analysis workflows.  
  - `FFRK/`: R package implementing the core FFRK method.  
  - `FFRK_Ablation/`: Codes for FFRK ablation experiments.  
  - `FFRKvsFFR/`: Codes for comparison experiments between FFRK and FFR.  
  - `Fig4+5/`: ArcGIS Pro project packages for Figures 4 and 5.  
  - `Fig5_Global/`: Codes for global prediction visualization in Figure 5.  
  - `Fig6_Evaluation_14models/`: Evaluation of 14 models (OK + UK + RK + ML + FFRK) and codes related to Figure 6.  
  - `Fig7&8_Sensitivity_KQ/`: Codes for K & Q sensitivity analysis.  
  - `FigS1-S3_Scatter/`: Codes for scatter plots of Supplementary Figures S1–S3.  
  - `Figx_Variance/`: Kriging variance analysis and visualization.  
  - `Figx_Variograms/`: Variogram modeling and plotting.  
  - `Tab2_Evaluiation_HIER/`: Codes for Table 2 (HIER) evaluation.  
  - `Tab2_Evaluiation_STK/`: Codes for Table 2 (STK) evaluation.  
  - `VIF_P-value_Feature Importance/`: Feature importance, VIF, and p-value analysis for the 9 covariates.  

- `Data/`: Input datasets for experiments.  
  - `Global/`: Input files for global prediction visualization.  
    - `rds/`: RDS files used for generating global predictions.  
    - `Global.xls`: Input table for global prediction visualization.  

  - `data_Cu.csv`: Input data for copper (Cu).  
  - `data_Pb.csv`: Input data for lead (Pb).  
  - `data_Zn.csv`: Input data for zinc (Zn).  
  - `Cu_Result.csv`: Baseline results for copper (Cu).  
  - `Cu_Hier_Result.csv`: HIER model results for copper (Cu).  
  - `Cu_STK_Result.csv`: STK model results for copper (Cu).  

- `Figures/`: Finalized figures.  
  - `Figures3_revised/`: Revised Figure 3.  
  - `Figures5_revised/`: Revised Figure 5.  
  - `FiguresS1-S3_revised/`: Revised Supplementary Figures S1–S3.  
  - `OKUK.jpg`: Global prediction maps for OK and UK.  

- `Results/`: Experimental outputs.  
  - `Endogeneity/`: Endogeneity analysis results.  
  - `FFRK_Ablation/`: Ablation experiment results.  
  - `FFRKvsFFR/`: Comparison results between FFRK and FFR.  
  - `Fig6_Evaluation_14models/`: Evaluation results of the 14 models (OK + UK + RK + ML + FFRK).  
  - `Fig7&8_Sensitivity_KQ/`: K & Q sensitivity analysis results.  
  - `Global/`: Global prediction results.  
  - `Hier_Evaluation/`: HIER evaluation results.  
  - `Scatter/`: Scatter plot results.  
  - `Tab2_Evaluiation_HIER/`: Table 2 (HIER) results.  
  - `Tab2_Evaluiation_STK/`: Table 2 (STK) results.  
  - `Universal Kriging(UK)/`: UK results.  
  - `VIF_P-value_Feature Importance/`: Feature importance, VIF, and p-value results for the 9 covariates.  
  - `Variance/`: Kriging variance results.  
  - `Variogram/`: Variogram analysis results.
