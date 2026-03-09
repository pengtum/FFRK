rm(list = ls()); gc()

library(parallel)
library(doParallel)
library(caret)
library(devtools)
library(progress)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))

data_dir <- file.path(workspace_dir, "Data")
out_dir  <- file.path(workspace_dir, "Results", "FFRK_Ablation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load_all(file.path(workspace_dir, "Code", "FFRK"))

x_field <- "DLONG"
y_field <- "DLAT"
target_var <- "Cu_ppm"

predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")

train_test_split_ratio <- 0.8
reproducible_seed <- 20250718

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


target_elements <- c("Cu", "Pb", "Zn")

for (elem in target_elements) {
  
  data_file <- file.path(data_dir, paste0("data_", elem, ".csv"))
  target_var <- paste0(elem, "_ppm")
  
  cat("Running element:", elem, "\n")
  cat("  data_file :", data_file, "\n")
  cat("  target_var:", target_var, "\n")
  
  data <- prepare_data(
    data_file,
    value_field = target_var,
    x_field = x_field,
    y_field = y_field,
    log_transform = TRUE
  )

  
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
  
  models_list <- c(paste0("FFRK_", names(base_methods_caret)))
  
  train_control <- trainControl(
    method = "CV",
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  

  features_full <- features
  ablation_list <- list(
    "IDW"                   = "f.idw",
    "IDW_Quantiles"          = c("f.idw", grep("^f\\.sda", names(features_full), value = TRUE)),
    "IDW_Similarity"         = c("f.idw", "f.gos"),
    "Quantiles_Similarity"   = c(grep("^f\\.sda", names(features_full), value = TRUE), "f.gos")
  )
  
  ablation_results <- list()
  
  for (scenario in names(ablation_list)) {
    geo_feats <- ablation_list[[scenario]]
    
    data_base   <- data[, c(x_field, y_field, target_var, predictor_fields)]
    data_ablate <- cbind(data_base, features_full[, geo_feats, drop = FALSE])
    
    formulas_ablate <- generate_formula(
      yvar         = target_var,
      features     = predictor_fields,
      geo_features = geo_feats
    )
    
    set.seed(reproducible_seed)
    res_ablate <- run_train_test_evaluation(
      data          = data_ablate,
      models_list   = models_list,
      target_var    = target_var,
      x_field       = x_field,
      y_field       = y_field,
      formula_main  = formulas_ablate$formula_main,
      formula_geo   = formulas_ablate$formula_geo,
      train_control = train_control,
      split_ratio   = train_test_split_ratio,
      seed          = reproducible_seed
    )
    
    metrics          <- res_ablate$metrics_summary
    metrics$scenario <- scenario
    ablation_results[[scenario]] <- metrics
  }
  
  ablation_df <- do.call(rbind, ablation_results)
  rownames(ablation_df) <- paste0(ablation_df$scenario, "_", ablation_df$model)
  ablation_df$scenario <- NULL
  
  print(ablation_df)
  
  output_file <- file.path(out_dir, paste0(elem, "_ablation_metrics.csv"))
  write.csv(ablation_df, file = output_file, row.names = TRUE)
}

stopCluster(cl)

rm(list = ls()); gc()
