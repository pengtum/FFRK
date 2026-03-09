rm(list = ls()); gc()

library(sp)
library(gstat)
library(automap)
library(rpart)
library(rpart.plot)
library(dplyr)
library(caret)
library(doParallel)
library(ggplot2)
library(progress)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(workspace_dir, "Data")

output_dir <- file.path(workspace_dir, "Results", "Tab2_Evaluiation_STK")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

target_elements <- c("Cu","Pb","Zn")
element_to_file <- list(
  Cu = "data_Cu.csv",
  Pb = "data_Pb.csv",
  Zn = "data_Zn.csv"
)

k_values <- 2:10

R2 <- function(o, p) {
  idx <- !is.na(o) & !is.na(p)
  if (sum(idx) < 2) return(NA_real_)
  1 - sum((o[idx] - p[idx])^2) / sum((o[idx] - mean(o[idx], na.rm = TRUE))^2)
}
MAE <- function(actual, predicted) mean(abs(actual - predicted), na.rm = TRUE)
RMSE <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))

find_tree_with_k_leaves <- function(formula, data, k) {
  initial_tree <- rpart(formula, data = data, method = "anova",
                        cp = 0.0001, minsplit = 2, minbucket = 1)
  if (length(unique(initial_tree$where)) <= k) return(initial_tree)
  cp_table <- initial_tree$cptable
  best_cp <- cp_table[which.min(abs(cp_table[, "nsplit"] - (k - 1))), "CP"]
  prune(initial_tree, cp = best_cp)
}

run_stk_evaluation <- function(train_data, test_data, k_val,
                               target_column, formula_OK, formula_partition) {
  
  hier_tree <- find_tree_with_k_leaves(formula_partition, data = train_data, k = k_val)
  train_data$region <- as.factor(predict(hier_tree, newdata = train_data, type = "vector"))
  test_data$region  <- as.factor(predict(hier_tree, newdata = test_data, type = "vector"))
  
  o <- test_data[[target_column]]
  p <- rep(NA, nrow(test_data))
  
  train_sp_global <- train_data
  coordinates(train_sp_global) <- ~ DLONG + DLAT
  train_sp_global$region <- train_data$region
  
  global_vario_fit <- tryCatch(
    autofitVariogram(formula_OK, train_sp_global, model = c("Sph")),
    error = function(e) NULL
  )
  
  unique_train_regions <- unique(train_data$region)
  
  for (r in unique_train_regions) {
    train_r <- train_data[train_data$region == r, ]
    test_r_idx <- which(test_data$region == r)
    
    if (length(test_r_idx) == 0) next
    
    model_to_use <- NULL
    
    if (nrow(train_r) >= 10) {
      train_r_sp <- train_r
      coordinates(train_r_sp) <- ~ DLONG + DLAT
      
      local_vario_fit <- tryCatch(
        autofitVariogram(formula_OK, train_r_sp, model = c("Sph")),
        error = function(e) NULL
      )
      
      if (!is.null(local_vario_fit$var_model)) {
        model_to_use <- local_vario_fit$var_model
      }
    }
    
    if (is.null(model_to_use) && !is.null(global_vario_fit$var_model)) {
      model_to_use <- global_vario_fit$var_model
      cat(paste0("Info: k=", k_val, ", Region '", r,
                 "' used global model due to insufficient data or fit failure.\n"))
    }
    
    # kriging
    if (!is.null(model_to_use)) {
      test_r_sp <- test_data[test_r_idx, ]
      coordinates(test_r_sp) <- ~ DLONG + DLAT
      
      train_region_sp <- train_sp_global[train_sp_global$region == r, ]
      if (nrow(train_region_sp) < 3) {
        train_region_sp <- train_sp_global
      }
      
      pred <- krige(formula_OK, train_region_sp, test_r_sp, model = model_to_use)
      p[test_r_idx] <- pred$var1.pred
    }
  }
  
  new_region_idx <- which(!(test_data$region %in% unique_train_regions))
  if (length(new_region_idx) > 0 && !is.null(global_vario_fit$var_model)) {
    test_new_sp <- test_data[new_region_idx, ]
    coordinates(test_new_sp) <- ~ DLONG + DLAT
    
    pred_new <- krige(formula_OK, train_sp_global, test_new_sp, model = global_vario_fit$var_model)
    p[new_region_idx] <- pred_new$var1.pred
    
    cat(paste0("Info: k=", k_val, ", ", length(new_region_idx),
               " test points in new regions predicted by global model.\n"))
  }
  
  data.frame(o = o, p = p, region = test_data$region)
}

set.seed(20250718)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

for (target_element in target_elements) {
  cat(paste0("\n>>> Analyzing Element: ", target_element, " <<<\n"))
  
  target_column <- paste0(target_element, "_ppm")
  formula_OK <- as.formula(paste(target_column, "~ 1"))
  formula_partition <- as.formula(paste(target_column, "~ DLONG + DLAT"))
  
  data_file <- file.path(data_dir, element_to_file[[target_element]])
  data <- read.csv(data_file)
  
  data[[target_column]] <- log(data[[target_column]])
  data <- na.omit(data)
  data$DLONG <- as.numeric(data$DLONG)
  data$DLAT  <- as.numeric(data$DLAT)
  
  train_idx <- createDataPartition(data[[target_column]], p = 0.8, list = FALSE)
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]

  all_results <- data.frame()
  
  pb <- progress_bar$new(
    total = length(k_values),
    format = paste0("[", target_element, "] [:bar] :percent [k=:k_val] eta: :eta")
  )
  
  for (k in k_values) {
    pb$tick(tokens = list(k_val = k))
    
    res_df <- run_stk_evaluation(
      train_data, test_data, k_val = k,
      target_column = target_column,
      formula_OK = formula_OK,
      formula_partition = formula_partition
    )
    
    res_df$k <- k
    res_df$model <- "STK"
    all_results <- rbind(all_results, res_df)
  }
  
  overall_metrics <- all_results %>%
    group_by(k, model) %>%
    summarise(
      overall_R2   = R2(o, p),
      overall_MAE  = MAE(o, p),
      overall_RMSE = RMSE(o, p),
      .groups = "drop"
    )
  
  region_metrics <- all_results %>%
    group_by(k, model, region) %>%
    summarise(
      region_R2   = R2(o, p),
      region_MAE  = MAE(o, p),
      region_RMSE = RMSE(o, p),
      n_points    = n(),
      .groups     = "drop"
    )
  
  region_instability <- region_metrics %>%
    group_by(k, model) %>%
    summarise(
      R2_sd    = sd(region_R2, na.rm = TRUE),
      RMSE_sd  = sd(region_RMSE, na.rm = TRUE),
      min_R2   = min(region_R2, na.rm = TRUE),
      max_R2   = max(region_R2, na.rm = TRUE),
      .groups  = "drop"
    )
  
  combined_metrics <- merge(overall_metrics, region_instability, by = c("k", "model")) %>%
    arrange(k)
  
  print(overall_metrics, row.names = FALSE)
  
  write.csv(all_results,
            file.path(output_dir, paste0(target_element, "_STK_predictions.csv")),
            row.names = FALSE)
  write.csv(region_metrics,
            file.path(output_dir, paste0(target_element, "_STK_region_metrics.csv")),
            row.names = FALSE)
  write.csv(combined_metrics,
            file.path(output_dir, paste0(target_element, "_STK_metrics.csv")),
            row.names = FALSE)
  
}

stopCluster(cl)

rm(list = ls())
gc()
