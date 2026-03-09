rm(list = ls()); gc()

library(car)          
library(randomForest)  
library(ggplot2)    
library(dplyr)
library(tidyr)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(workspace_dir, "Results", "VIF_P-value_Feature Importance")

target_elements <- c("Cu", "Pb", "Zn")
predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")

x_field <- "DLONG"
y_field <- "DLAT"

load_and_clean <- function(filepath, val_col) {
  df <- read.csv(filepath)
  df[[val_col]] <- log(df[[val_col]])
  df <- na.omit(df)
  return(df)
}

final_stats_summary <- data.frame()

for (elem in target_elements) {
  data_dir <- file.path(workspace_dir,"Data", "0_Origin")
  
  data_file <- file.path(data_dir, paste0("data_", elem, ".csv"))
  
  target_var <- paste0(elem, "_ppm")
  
  data <- load_and_clean(data_file, target_var)
  
  Y <- data[[target_var]]
  X <- data[, predictor_fields]
  
  # VIF
  f_str <- paste(target_var, "~", paste(predictor_fields, collapse = "+"))
  lm_model <- lm(as.formula(f_str), data = data)
  
  vif_values <- car::vif(lm_model)
  
  # P-value)
  lm_sum <- summary(lm_model)
  p_values <- lm_sum$coefficients[-1, "Pr(>|t|)"]
  coef_est <- lm_sum$coefficients[-1, "Estimate"]
  
  # Feature Importance
  set.seed(123)
  rf_model <- randomForest(x = X, y = Y, importance = TRUE, ntree = 500)
  imp_data <- importance(rf_model)
  rf_imp <- imp_data[, "%IncMSE"]
  

  vars <- names(vif_values)
  
  elem_stats <- data.frame(
    Element = elem,
    Variable = vars,
    VIF = as.numeric(vif_values[vars]),
    LM_P_Value = as.numeric(p_values[vars]),
    LM_Coefficient = as.numeric(coef_est[vars]),
    RF_Importance_IncMSE = as.numeric(rf_imp[vars])
  )
  
  elem_stats$Significant_LM <- ifelse(elem_stats$LM_P_Value < 0.05, "Yes", "No")
  
  final_stats_summary <- rbind(final_stats_summary, elem_stats)
  
  # plot
  p <- ggplot(elem_stats, aes(x = reorder(Variable, RF_Importance_IncMSE), y = RF_Importance_IncMSE)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(
      title = paste0("Feature Importance for ", elem, " (Random Forest)"),
      subtitle = "Metric: % Increase in MSE",
      x = "Covariates",
      y = "Importance (%IncMSE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(filename = file.path(out_dir, paste0("VIF_P-value_Feature Importance", elem, ".png")), p, width = 6, height = 4, bg = "white")
}


num_cols <- sapply(final_stats_summary, is.numeric)
final_stats_summary[num_cols] <- round(final_stats_summary[num_cols], 4)

write.csv(final_stats_summary,file.path(out_dir, "VIF_P-value_Feature Importance.csv"), row.names = FALSE)
