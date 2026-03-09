rm(list = ls()); gc()

library(parallel)
library(doParallel)
library(caret)
library(devtools)
library(progress)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
load_all(file.path(workspace_dir, "Code", "FFRK"))
out_dir <- file.path(workspace_dir, "Results", "Fig6_Evaluation_14models")

target_elements <- c("Cu", "Pb", "Zn")

x_field <- "DLONG"
y_field <- "DLAT"


predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")

train_test_split_ratio <- 0.8
reproducible_seed <- 20250718


cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


all_results_list <- list() 
all_elements_metrics <- data.frame()

for (elem in target_elements) {
  
  cat(paste0(">>> Processing: ", elem, "\n"))
  
  data_file <- paste0("data_", elem, ".csv") 
  target_var <- paste0(elem, "_ppm")
  
  data <- prepare_data(
    data_file,
    value_field = target_var,
    x_field = x_field,
    y_field = y_field,
    log_transform = TRUE
  )
  
  range(data[[target_var]], na.rm = TRUE)
  
  features <- getfeature(
    xo = data[[x_field]],
    yo = data[[y_field]],
    zo = data[[target_var]],
    xp = data[[x_field]],
    yp = data[[y_field]],
    k_neighbors = 15, 
    num_quantiles = 11
  )
  
  data <- cbind(data, features)
  
  formulas <- generate_formula(
    yvar = target_var,
    features = predictor_fields,
    geo_features = names(features)
  )
  
  
  base_methods_caret <- c(
    lm = "lm",
    rf = "rf",
    rpart = "rpart",
    svmRadial = "svmRadial"
  )
  
  models_list <- c("ok", # Ordinary Kriging 
                   "uk", # Universal Kriging 
                   names(base_methods_caret),  # ML：lm, rf, rpart, svmRadial
                   paste0("rk_", names(base_methods_caret)),
                   paste0("FFRK_", names(base_methods_caret)))  
  
  
  train_control <- trainControl(
    method = "CV", 
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  
  cat(paste0(">>> Starting training & evaluation: [", elem, "]...\n"))
  
  res <- run_train_test_evaluation(
    data = data,
    models_list = models_list,
    target_var = target_var,
    x_field = x_field,
    y_field = y_field,
    formula_main = formulas$formula_main,
    formula_geo = formulas$formula_geo,
    train_control = train_control,
    split_ratio = train_test_split_ratio, 
    seed = reproducible_seed             
  )
  
  write.csv(res$detailed_results, file.path(out_dir, paste0("Evaluation_14models_", elem, ".csv")), row.names = FALSE)
  
  temp_res <- res$detailed_results
  temp_res$Element <- elem
  all_results_list[[elem]] <- temp_res
  
  current_metrics <- res$metrics_summary
  current_metrics$Element <- elem 
  
  all_elements_metrics <- rbind(all_elements_metrics, current_metrics)
  print(all_elements_metrics)
}


stopCluster(cl)
registerDoSEQ()


