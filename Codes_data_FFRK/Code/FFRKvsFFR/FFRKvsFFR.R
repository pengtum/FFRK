rm(list = ls()); gc()

library(caret)
library(parallel)
library(doParallel)
library(devtools)
library(automap)
library(gstat)
library(sp)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
load_all(file.path(workspace_dir, "Code", "FFRK"))
out_dir <- file.path(workspace_dir, "Results", "FFRKvsFFR")
data_dir <- file.path(workspace_dir,"Data")

get_model_info <- function(model_type, formula_main = NULL, formula_geo = NULL) {
  base_methods <- c(lm = "lm", rf = "rf", rpart = "rpart", svmRadial = "svmRadial")
  
  if (model_type == "ok") return(list(type = "ok", caret_method = NA, use_formula = NULL))
  if (model_type == "uk") return(list(type = "uk", caret_method = NA, use_formula = NULL))
  
  # FFRK
  if (grepl("^FFRK_", model_type)) {
    base <- sub("^FFRK_", "", model_type)
    return(list(type = "FFRK", caret_method = base_methods[base], use_formula = formula_geo))
  } 
  # FFR
  if (grepl("^FFR_", model_type)) {
    base <- sub("^FFR_", "", model_type)
    return(list(type = "ml_geo", caret_method = base_methods[base], use_formula = formula_geo))
  }
  
  stop("Unknown model type: ", model_type)
}

run_fold <- function(train_data, test_data, model_type, fold_i,
                     target_var, x_field, y_field,
                     train_control, formula_main, formula_geo, seed = 20250718) {
  
  train_data <- na.omit(train_data); test_data  <- na.omit(test_data)
  if(nrow(train_data)==0 || nrow(test_data)==0) return(list(o=NA, p=NA))
  
  test_data_df <- as.data.frame(test_data)
  o <- test_data[[target_var]]
  
  train_sp <- train_data; test_sp <- test_data
  coordinates(train_sp) <- as.formula(paste("~", x_field, "+", y_field))
  coordinates(test_sp)  <- as.formula(paste("~", x_field, "+", y_field))
  
  info <- get_model_info(model_type, formula_main, formula_geo)
  
  current_seed <- seed + fold_i + sum(as.integer(charToRaw(model_type)))
  set.seed(current_seed)
  
  p <- rep(NA, length(o))
  
  if (info$type == "ml_geo") {
    # FFR
    fit <- train(info$use_formula, data=train_data, method=info$caret_method, trControl=train_control)
    p <- predict(fit, newdata=test_data_df)
    
  } else if (info$type == "FFRK") {
    # FFRK
    fit <- train(info$use_formula, data=train_data, method=info$caret_method, trControl=train_control)
    ml_pred <- predict(fit, newdata=test_data_df)
    
    pred_train <- predict(fit, newdata=train_data)
    train_data$res <- train_data[[target_var]] - pred_train
    
    train_res_sp <- train_data
    coordinates(train_res_sp) <- as.formula(paste("~", x_field, "+", y_field))
    
    vario_fit <- fit_variogram_robust(train_res_sp, res ~ 1)
    
    if(!is.null(vario_fit)) {
      pred_krig_obj <- krige(res ~ 1, train_res_sp, test_sp, model=vario_fit)
      p <- ml_pred + pred_krig_obj$var1.pred
    } else {
      p <- ml_pred
    }
  }
  
  return(list(o = o, p = p, DLONG = test_data[[x_field]], DLAT = test_data[[y_field]]))
}

run_train_test_evaluation <- function(data, models_list, target_var, x_field, y_field,
                                      formula_main, formula_geo, train_control,
                                      split_ratio = 0.8, seed = 20250718) {
  
  set.seed(seed)
  train_idx <- createDataPartition(data[[target_var]], p = split_ratio, list = FALSE)
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]
  
  metrics_summary <- data.frame()
  
  for (m in models_list) {
    cat(paste("  Model:", m, "...\n"))
    res <- run_fold(train_data, test_data, m, 1, target_var, x_field, y_field,
                    train_control, formula_main, formula_geo)
    
    if (!all(is.na(res$p))) {
      metrics_summary <- rbind(metrics_summary, data.frame(
        model = m,
        R2    = metric_r2(res$o, res$p),
        MAE   = metric_mae(res$o, res$p),
        RMSE  = metric_rmse(res$o, res$p)
      ))
    }
  }
  return(list(metrics_summary = metrics_summary))
}


target_elements <- c("Cu", "Pb", "Zn")
predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")
x_field <- "DLONG"; y_field <- "DLAT"

# FFRK vs FFR
models_to_test <- c(
  "FFRK_lm", "FFR_lm",            # Linear Regression
  "FFRK_rpart", "FFR_rpart",      # Decision Tree
  "FFRK_rf", "FFR_rf",            # Random Forest
  "FFRK_svmRadial", "FFR_svmRadial" # SVM
)

cl <- makeCluster(detectCores() - 1); registerDoParallel(cl)
comparison_results <- data.frame()

for (elem in target_elements) {
  cat(paste0("\n>>> Analyzing Element: ", elem, " <<<\n"))
  
  data_file <- file.path(data_dir, paste0("data_", elem, ".csv"))
  target_var <- paste0(elem, "_ppm")
  
  data <- prepare_data(data_file, value_field = target_var, x_field = x_field, y_field = y_field, log_transform = TRUE)
  features <- getfeature(xo=data[[x_field]], yo=data[[y_field]], zo=data[[target_var]], 
                         xp=data[[x_field]], yp=data[[y_field]], k_neighbors = 15, num_quantiles = 11)
  data_model <- na.omit(cbind(data, features))
  
  formulas <- generate_formula(target_var, predictor_fields, names(features))
  
  res <- run_train_test_evaluation(
    data = data_model,
    models_list = models_to_test,
    target_var = target_var,
    x_field = x_field,
    y_field = y_field,
    formula_main = formulas$formula_main,
    formula_geo = formulas$formula_geo,
    train_control = trainControl(method="CV", number=5, verboseIter=FALSE, allowParallel=TRUE),
    split_ratio = 0.8,
    seed = 20250718
  )
  
  curr <- res$metrics_summary
  curr$Element <- elem
  comparison_results <- rbind(comparison_results, curr)
}

stopCluster(cl); registerDoSEQ()

comparison_results <- comparison_results[order(comparison_results$Element), ]
print(comparison_results[, c("Element", "model", "R2", "MAE", "RMSE")])
write.csv(comparison_results, file.path(out_dir, "FFRKvsFFR.csv"), row.names = FALSE)