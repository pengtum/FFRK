rm(list = ls()); gc()

library(parallel)
library(doParallel)
library(caret)
library(devtools)
library(progress)

target_elements <- c("Cu", "Pb", "Zn")

x_field <- "DLONG"
y_field <- "DLAT"

predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")

train_test_split_ratio <- 0.8
reproducible_seed <- 20250718

k_values      <- c(15, 20, 30, 50)
num_quantiles <- seq(2, 21, by = 1)


script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(workspace_dir, "Data")
ffrk_dir <- file.path(workspace_dir, "Code", "FFRK")
load_all(ffrk_dir)

out_dir <- file.path(workspace_dir, "Results", "Fig7&8_Sensitivity_KQ")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

base_methods_caret <- c(
  lm        = "lm",
  rf        = "rf",
  rpart     = "rpart",
  svmRadial = "svmRadial"
)

models_list  <- paste0("FFRK_", names(base_methods_caret))

train_control <- trainControl(
  method        = "CV",
  verboseIter   = FALSE,
  allowParallel = TRUE
)

metrics_all <- data.frame()


for (elem in target_elements) {

  data_path<- file.path(data_dir, paste0("data_", elem, ".csv"))
  
  target_var <- paste0(elem, "_ppm")
  
  raw_data <- prepare_data(
    data_path,
    value_field   = target_var,
    x_field       = x_field,
    y_field       = y_field,
    log_transform = TRUE
  )

  
  for (k in k_values) {
    for (q in num_quantiles) {
      
      if (k < q) {
        cat(sprintf("⚠️ 跳过 k=%d, q=%d (因为 k < q)\n", k, q))
        next
      }
      
      cat("\n==============================\n")
      cat(sprintf(">>> 正在处理 k = %d, q = %d\n", k, q))
      cat("==============================\n")
      
      features <- getfeature(
        xo            = raw_data[[x_field]],
        yo            = raw_data[[y_field]],
        zo            = raw_data[[target_var]],
        xp            = raw_data[[x_field]],
        yp            = raw_data[[y_field]],
        k_neighbors   = k,
        num_quantiles = q
      )
      
      data <- cbind(raw_data, features)
      
      formulas <- generate_formula(
        yvar         = target_var,
        features     = predictor_fields,
        geo_features = names(features)
      )
      
      res <- run_train_test_evaluation(
        data          = data,
        models_list   = models_list,
        target_var    = target_var,
        x_field       = x_field,
        y_field       = y_field,
        formula_main  = formulas$formula_main,
        formula_geo   = formulas$formula_geo,
        train_control = train_control,
        split_ratio   = train_test_split_ratio,
        seed          = reproducible_seed
      )
      
      m <- res$metrics_summary
      m$k <- k
      m$q <- q
      metrics_all <- rbind(metrics_all, m)
    }
  }
  
  rownames(metrics_all) <- NULL
  print(metrics_all)
  
  write.csv(
    metrics_all,
    file = file.path(out_dir, paste0(elem,"_k_q_sensitivity_metrics.csv")),
    row.names = FALSE
  )
  
}



stopCluster(cl)
registerDoSEQ()
rm(list = ls()); gc()
