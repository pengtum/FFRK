rm(list = ls()); gc()

library(sp)
library(gstat)
library(automap)
library(randomForest)
library(rpart)
library(rpart.plot)
library(dplyr)
library(caret)
library(doParallel)
library(ggplot2)
library(progress)

target_elements <- c("Cu","Pb","Zn")

element_to_file <- list(
  Cu = "data_Cu.csv",
  Pb = "data_Pb.csv",
  Zn = "data_Zn.csv"
)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(workspace_dir,"Data")
out_dir  <- file.path(workspace_dir,"Results","Tab2_Evaluiation_HIER") 

predictors <- c("Dlith", "Dfault", "Slope", "Water", "NDVI", "MainRd", "Road", "SOC", "pH")

base_methods_caret <- c(lm = "lm",
                        rf = "rf",
                        rpart = "rpart",
                        svmRadial = "svmRadial")
models_list <- paste0("hier_", names(base_methods_caret))

k_values <- 2:10 
train_control <- trainControl(method = "CV", allowParallel = TRUE)

R2 <- function(o, p) {
  idx <- !is.na(o) & !is.na(p)
  if (sum(idx) < 2) return(NA)
  1 - sum((o[idx] - p[idx])^2) / sum((o[idx] - mean(o[idx], na.rm = TRUE))^2)
}
MAE <- function(actual, predicted) mean(abs(actual - predicted), na.rm = TRUE)
RMSE <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))

find_tree_with_k_leaves <- function(formula, data, k) {
  initial_tree <- rpart(formula, data = data, method = "anova",
                        cp = 0.0001, minsplit = 2, minbucket = 1)
  if(length(unique(initial_tree$where)) <= k) return(initial_tree)
  cp_table <- initial_tree$cptable
  best_cp <- cp_table[which.min(abs(cp_table[,"nsplit"] - (k-1))), "CP"]
  prune(initial_tree, cp = best_cp)
}

run_hierarchical_model <- function(train_data, test_data, model_info, k_val, namey, formula_partition) {
  
  hier_tree <- find_tree_with_k_leaves(formula_partition, data = train_data, k = k_val)
  
  train_data$region <- as.factor(predict(hier_tree, newdata = train_data, type = "vector"))
  test_data$region  <- as.factor(predict(hier_tree, newdata = test_data, type = "vector"))
  
  p_plot <- ggplot() +
    geom_point(data = train_data, aes(x=DLONG, y=DLAT, color=region), alpha=0.5, shape=1) +
    geom_point(data = test_data,  aes(x=DLONG, y=DLAT, color=region), alpha=0.8, shape=16) +
    ggtitle(paste("k =", k_val)) +
    theme_bw() + theme(legend.position = "none")
  print(p_plot)
  
  o <- test_data[[namey]]
  p <- rep(NA, nrow(test_data))
  
  global_model <- train(model_info$use_formula, data = train_data,
                        method = model_info$caret_method, trControl = train_control)
  
  unique_train_regions <- unique(train_data$region)
  
  for(r in unique_train_regions) {
    train_r <- train_data[train_data$region == r, ]
    
    if(nrow(train_r) < 10) {
      model_fit <- global_model
    } else {
      model_fit <- train(model_info$use_formula, data = train_r,
                         method = model_info$caret_method, trControl = train_control)
    }
    
    test_r_idx <- which(test_data$region == r)
    if(length(test_r_idx) > 0) {
      p[test_r_idx] <- predict(model_fit, newdata = test_data[test_r_idx, ])
    }
  }
  
  new_region_idx <- which(!(test_data$region %in% unique_train_regions))
  if (length(new_region_idx) > 0) {
    p[new_region_idx] <- predict(global_model, newdata = test_data[new_region_idx, ])
    cat("Info: Found", length(new_region_idx),
        "test points in regions not present in the training set. Used global model for prediction.\n")
  }
  
  list(o = o, p = p)
}

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

for (elem in target_elements) {
  
  namey <- paste0(elem, "_ppm")
  
  data_file <- file.path(data_dir, element_to_file[[elem]])
  data <- read.csv(data_file)
  
  formula_ml <- as.formula(paste(namey, "~", paste(predictors, collapse = "+")))
  formula_partition <- as.formula(paste(namey, "~ DLONG + DLAT"))
  
  data[[namey]] <- log(data[[namey]])
  data <- na.omit(data)
  data$DLONG <- as.numeric(data$DLONG)
  data$DLAT  <- as.numeric(data$DLAT)
  
  set.seed(20250718)
  train_idx <- createDataPartition(data[[namey]], p = 0.8, list = FALSE)
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]
  
  metrics <- data.frame()
  all_results <- data.frame()
  
  total_runs <- length(k_values) * length(models_list)
  pb <- progress_bar$new(
    total = total_runs,
    format = paste0("[", elem, "] [:bar] :percent [k=:k_val, model=:model] eta: :eta")
  )
  
  for (k in k_values) {
    for (m in models_list) {
      pb$tick(tokens = list(k_val = k, model = m))
      
      base_m <- sub("^hier_", "", m)
      model_info <- list(
        type = "hier",
        caret_method = base_methods_caret[base_m],
        use_formula = formula_ml
      )
      
      res <- run_hierarchical_model(
        train_data, test_data, model_info, k_val = k,
        namey = namey, formula_partition = formula_partition
      )
      
      metrics <- rbind(metrics, data.frame(
        element = elem,
        k = k,
        model = m,
        R2 = R2(res$o, res$p),
        MAE = MAE(res$o, res$p),
        RMSE = RMSE(res$o, res$p)
      ))
      
      all_results <- rbind(all_results, data.frame(
        element = elem, k = k, model = m, o = res$o, p = res$p
      ))
    }
  }
  
  print(metrics, row.names = FALSE)
  
  write.csv(metrics,
            file = file.path(out_dir, paste0(namey, "_Hierarchical_Metrics.csv")),
            row.names = FALSE)
  
  write.csv(all_results,
            file = file.path(out_dir, paste0(namey, "_Hierarchical_Predictions.csv")),
            row.names = FALSE)
}

stopCluster(cl)

rm(list = ls()); gc()
