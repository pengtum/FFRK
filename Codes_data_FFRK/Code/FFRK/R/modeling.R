#' Fit Variogram with Robust Handling
#'
#' This function attempts to fit a variogram model using `autofitVariogram`.
#' It catches errors and warnings and returns `NULL` if fitting fails.
#'
#' @param vario_data A SpatialPointsDataFrame or gstat-compatible object containing data.
#' @param formula A formula specifying the response and predictors, e.g., `z ~ 1`.
#'
#' @return A fitted variogram model (object of class `variogramModel`) if successful, otherwise NULL.
#'
#' @examples
#' \dontrun{
#' vm <- fit_variogram_robust(my_spatial_data, z ~ 1)
#' if(!is.null(vm)) print(vm)
#' }
#'
#' @export
fit_variogram_robust <- function(vario_data, formula) {
  fit_auto <- tryCatch({
    autofitVariogram(formula, vario_data, model = c("Sph", "Exp"),    # 不要用 Gau
                     kappa = c(0.5,1.0,1.5))

  }, error = function(e) {
    NULL
  }, warning = function(w) {
    NULL
  })

  if (!is.null(fit_auto) && !is.null(fit_auto$var_model)) {
    print(fit_auto)
    plot(fit_auto, cutoff = 100)
    return(fit_auto$var_model)
  }
  return(NULL)
}

#' Run a Single Fold of Spatial Modeling
#'
#' This function handles one fold of spatial modeling using ML, RK, FFRK, or hierarchical methods.
#'
#' @param train_data A data.frame containing training data.
#' @param test_data A data.frame containing testing data.
#' @param model_type A string indicating the model type (e.g., "rf", "rk_rf", "FFRK_rf", "hier_rf").
#' @param fold_i Integer. The current fold index.
#' @param target_var Character. The name of the target variable.
#' @param x_field Character. The name of the X/longitude field.
#' @param y_field Character. The name of the Y/latitude field.
#' @param train_control A caret trainControl object.
#' @param formula_main Formula for standard models (ML, RK, hier).
#' @param formula_geo Formula for FFRK models.
#'
#' @return A list with observed values, predicted values, coordinates and variogram parameters.
#'
#' @importFrom rpart rpart
#' @importFrom gstat krige
#' @importFrom automap autofitVariogram
#'
#' @export
run_fold <- function(train_data, test_data, model_type, fold_i,
                                   target_var, x_field, y_field,
                                   train_control,
                                   formula_main, formula_geo) {

  train_data <- na.omit(train_data)
  test_data  <- na.omit(test_data)

  test_coords <- test_data[, c(x_field, y_field)]

  if(nrow(train_data)==0 || nrow(test_data)==0) {
    return(list(o=rep(NA, nrow(test_data)),
                p=rep(NA, nrow(test_data)),
                p_var=rep(NA, nrow(test_data)),
                vario_model_str=rep(NA, nrow(test_data)),
                vario_psill_str=rep(NA, nrow(test_data)),
                vario_range_str=rep(NA, nrow(test_data)),
                vario_kappa_str=rep(NA, nrow(test_data))))
  }

  o <- test_data[[target_var]]

  if (length(o) == 0) {
    warning("Test data has no target values, skipping fold.")
    return(list(o=NA, p=NA,
                vario_model_str=NA,
                vario_psill_str=NA,
                vario_range_str=NA,
                vario_kappa_str=NA))
  }

  train_sp <- train_data
  test_sp  <- test_data
  coordinates(train_sp) <- stats::as.formula(paste("~", x_field, "+", y_field))
  coordinates(test_sp)  <- stats::as.formula(paste("~", x_field, "+", y_field))

  info <- get_model_info(model_type)

  vario_model_str <- NA
  vario_psill_str <- NA
  vario_range_str <- NA
  vario_kappa_str <- NA

  p <- rep(NA, length(o))
  p_var <- rep(NA, length(o))
  if(info$type == "ok") {
    # OK
    var_formula <- stats::as.formula(paste(target_var, "~ 1"))
    res_auto <- fit_variogram_robust(train_sp, var_formula)
    if(!is.null(res_auto)) {
      vario_model_str <- paste(res_auto$model, collapse=";")
      vario_psill_str <- paste(round(res_auto$psill,4), collapse=";")
      vario_range_str <- paste(round(res_auto$range,4), collapse=";")
      vario_kappa_str <- paste(round(res_auto$kappa,4), collapse=";")

      pred.krig <- krige(var_formula, train_sp, test_sp, model=res_auto)
      p <- pred.krig$var1.pred
      p_var <- pred.krig$var1.var

      saveRDS(res_auto, paste0("ok_variogram_", Sys.Date(), ".rds"))
    } else {
      warning(paste(model_type, "variogram fitting failed, predictions set to NA."))
    }

  } else if(info$type == "uk") {
    # UK
    # Formula: Target ~ x + y + x^2 + y^2 + x*y

    uk_formula_str <- paste(target_var, "~",
                            x_field, "+", y_field, "+",
                            "I(", x_field, "^2) +",
                            "I(", y_field, "^2) +",
                            "I(", x_field, "*", y_field, ")")
    var_formula <- stats::as.formula(uk_formula_str)

    res_auto <- fit_variogram_robust(train_sp, var_formula)

    if(!is.null(res_auto)) {
      vario_model_str <- paste(res_auto$model, collapse=";")
      vario_psill_str <- paste(round(res_auto$psill,4), collapse=";")
      vario_range_str <- paste(round(res_auto$range,4), collapse=";")
      vario_kappa_str <- paste(round(res_auto$kappa,4), collapse=";")

      pred.krig <- krige(var_formula, train_sp, test_sp, model=res_auto)
      p <- pred.krig$var1.pred
      p_var <- pred.krig$var1.var

      # saveRDS(res_auto, paste0("uk_variogram_", Sys.Date(), ".rds"))
    } else {
      warning(paste(model_type, "variogram fitting failed, predictions set to NA."))
    }

  } else if(info$type == "ml") {
    # ML
    model_fit <- train(formula_main, data=train_data,
                       method=info$caret_method, trControl=train_control, na.action=na.omit)
    p <- predict(model_fit, newdata=test_data)
    p_var <- rep(0, length(p))

    final_model <- train(formula_main, data=data, method=info$caret_method, trControl=train_control, na.action=na.omit)
    saveRDS(final_model, paste0("final_model_", info$caret_method, "_", Sys.Date(), ".rds"))


  } else if(info$type == "rk") {
    # rk
    model_fit <- train(formula_main, data=train_data,
                       method=info$caret_method, trControl=train_control, na.action=na.omit)

    pred_train <- predict(model_fit, newdata=train_data)
    train_data$res <- train_data[[target_var]] - pred_train
    Aimp <- train_data$res
    train_res_sp <- train_data
    coordinates(train_res_sp) <- stats::as.formula(paste("~", x_field, "+", y_field))
    proj4string(train_res_sp) <- CRS("+proj=longlat +datum=WGS84")
    train_res_sp <- spTransform(train_res_sp, CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs"))
    proj4string(test_sp) <- CRS("+proj=longlat +datum=WGS84")
    test_sp <- spTransform(test_sp, CRS(proj4string(train_res_sp)))
    summary(train_res_sp$res)
    vario_fit <- fit_variogram_robust(train_res_sp, res ~ 1)

    if(!is.null(vario_fit)) {
      pred.krig <- krige(res~1, train_res_sp, test_sp, model=vario_fit)
      p <- predict(model_fit, newdata=test_data) + pred.krig$var1.pred
      p_var <- pred.krig$var1.var

      saveRDS(model_fit, paste0("rk_ml_model_", Sys.Date(), ".rds"))
      saveRDS(vario_fit, paste0("rk_variogram_", Sys.Date(), ".rds"))
    } else {
      warning(paste(model_type, "variogram fitting failed, using only ML prediction."))
      p <- predict(model_fit, newdata=test_data)
    }

  } else if(info$type == "FFRK") {
    # FFRK
    model_fit <- train(formula_geo, data=train_data,
                       method=info$caret_method, trControl=train_control, na.action=na.omit)
    pred_train <- predict(model_fit, newdata=train_data)
    train_data$res <- train_data[[target_var]] - pred_train
    train_res_sp <- train_data
    saveRDS(train_res_sp,
            file = paste0("train_res_", info$type, ".rds"))
    coordinates(train_res_sp) <- stats::as.formula(paste("~", x_field, "+", y_field))

    vario_fit <- fit_variogram_robust(train_res_sp, res ~ 1)
    if(!is.null(vario_fit)) {
      pred.krig <- krige(res ~ 1, train_res_sp, test_sp, model=vario_fit)
      p <- predict(model_fit, newdata=test_data) + pred.krig$var1.pred
      p_var <- pred.krig$var1.var

      # saveRDS(model_fit, paste0("ffrk_ml_model_", info$caret_method, "_fold", fold_i, "_", Sys.Date(), ".rds"))
      # saveRDS(vario_fit, paste0("ffrk_variogram_", info$caret_method, "_fold", fold_i, "_", Sys.Date(), ".rds"))
    } else {
      warning(paste(model_type, "variogram fitting failed, using only ffRK ML prediction."))
      p <- predict(model_fit, newdata=test_data)
    }

    if (fold_i == 1) {
      final_ffrk_model <- train(formula_geo, data=data,  # data 为全局完整数据
                               method=info$caret_method, trControl=train_control, na.action=na.omit)
      pred_train_full <- predict(final_ffrk_model, newdata=data)
      data$res <- data[[target_var]] - pred_train_full
      data_sp <- data
      coordinates(data_sp) <- stats::as.formula(paste("~", x_field, "+", y_field))
      final_vario_fit <- fit_variogram_robust(data_sp, res ~ 1)

      data_for_final_model <- data

      final_ffrk_model <- train(formula_geo, data=data_for_final_model,  # data 为全局完整数据
                               method=info$caret_method, trControl=train_control, na.action=na.omit)
      pred_train_full <- predict(final_ffrk_model, newdata=data_for_final_model)

      data_for_final_model$res <- data_for_final_model[[target_var]] - pred_train_full

      data_sp <- data_for_final_model
      coordinates(data_sp) <- stats::as.formula(paste("~", x_field, "+", y_field))
      final_vario_fit <- fit_variogram_robust(data_sp, res ~ 1)

      if(!is.null(final_vario_fit)) {
        saveRDS(final_ffrk_model, paste0("final_ffrk_ml_model_", info$caret_method, "_", Sys.Date(), ".rds"))
        saveRDS(final_vario_fit, paste0("final_ffrk_variogram_", info$caret_method, "_", Sys.Date(), ".rds"))
      } else {
        warning("Final FFRK variogram fitting failed, saving only ML model.")
        saveRDS(final_ffrk_model, paste0("final_ffrk_ml_model_", info$caret_method, "_", Sys.Date(), ".rds"))
      }
    }
  }


  else if(info$type == "hier") {
    # hier
    unique_regions <- unique(train_data$region)
    p <- numeric(nrow(test_data))

    for(r in unique_regions) {
      train_r <- train_data[ train_data$region == r, ]
      if(nrow(train_r) == 0) next

      model_fit <- train(formula_main, data=train_r,
                         method=info$caret_method, trControl=train_control, na.action=na.omit)

      saveRDS(model_fit, paste0("hier_model_region_", r, "_", Sys.Date(), ".rds"))

      test_r_idx <- which(test_data$region == r)
      if(length(test_r_idx) > 0) {
        p[test_r_idx] <- predict(model_fit, newdata=test_data[test_r_idx, ])
      }
    }
  }

  return(list(
    o = o,
    p = p,
    p_var = p_var,
    DLONG = test_coords[[x_field]],
    DLAT  = test_coords[[y_field]],
    vario_model_str  = rep(vario_model_str, length(o)),
    vario_psill_str  = rep(vario_psill_str, length(o)),
    vario_range_str  = rep(vario_range_str, length(o)),
    vario_kappa_str  = rep(vario_kappa_str, length(o))
  ))
}

#' Run Cross-Validation Workflow
#'
#' Runs repeated k-fold cross-validation for given models and data, including
#' region-based partitioning and result aggregation.
#'
#' @param data A data.frame with all input data.
#' @param models_list A character vector of model identifiers.
#' @param num_repeats Integer, number of repetitions.
#' @param k_folds Integer, number of folds.
#' @param target_var Character, target variable name.
#' @param x_field Character, X/longitude field.
#' @param y_field Character, Y/latitude field.
#' @param formula_main Formula for standard models.
#' @param formula_geo Formula for FFRK models.
#' @param train_control caret trainControl object.
#'
#' @return A list with final_results (data.frame) and cv_results (list).
#'
#' @export
run_cross_validation <- function(data,
                                 models_list,
                                 num_repeats,
                                 k_folds,
                                 target_var,
                                 x_field,
                                 y_field,
                                 formula_main,
                                 formula_geo,
                                 train_control) {
  # 初始化存储
  init <- initialize_cv_storage(models_list, num_repeats, k_folds)
  final_results <- init$final_results
  cv_results <- init$cv_results
  pb <- init$pb

  for (rep in seq_len(num_repeats)) {
    cat("Starting repetition:", rep, "\n")
    hier_tree <- rpart(stats::as.formula(paste(target_var, "~", x_field, "+", y_field)),
                       data = data, method = "anova")
    data$region <- as.factor(hier_tree$where)

    set.seed(123 + rep)

    folds <- createFolds(data[[target_var]], k = k_folds, list = TRUE, returnTrain = FALSE)

    for (m in models_list) {
      all_o <- c(); all_p <- c()
      for (fold_i in seq_len(k_folds)) {
        test_idx <- folds[[fold_i]]
        test_data <- data[test_idx, ]
        train_data <- data[-test_idx, ]

        res <- run_fold(
          train_data, test_data,
          model_type = m,
          fold_i = fold_i,
          target_var = target_var,
          x_field = x_field,
          y_field = y_field,
          train_control = train_control,
          formula_main = formula_main,
          formula_geo = formula_geo
        )

        all_o <- c(all_o, res$o)
        all_p <- c(all_p, res$p)

        fold_results <- data.frame(
          model = rep(m, length(res$o)),
          fold = rep(fold_i, length(res$o)),
          o = res$o,
          p = res$p,
          DLONG = res[[x_field]],
          DLAT = res[[y_field]],
          autofit_model = res$vario_model_str,
          autofit_psill = res$vario_psill_str,
          autofit_range = res$vario_range_str,
          autofit_kappa = res$vario_kappa_str,
          stringsAsFactors = FALSE
        )
        final_results <- rbind(final_results, fold_results)

        pb$tick(tokens = list(
          model_name = m,
          iteration = fold_i,
          total_fold = k_folds,
          rep = rep,
          total_rep = num_repeats
        ))
      }
      cv_results[[m]] <- rbind(cv_results[[m]],
                               data.frame(
                                 R2 = metric_r2(all_o, all_p),
                                 MAE = metric_mae(all_o, all_p),
                                 RMSE = metric_rmse(all_o, all_p)
                               ))
    }
  }

  return(list(final_results = final_results, cv_results = cv_results))
}

#' Run a Single Train-Test Split Evaluation
#'
#' Replaces k-fold cross-validation with a single, fixed train-test split.
#' It partitions the data, iterates through a list of models, runs each model
#' using run_fold, and aggregates the results. The random seed for partitioning
#' is fixed to ensure reproducibility, matching the logic from the user's Code 2.
#'
#' @param data A data.frame with all input data.
#' @param models_list A character vector of model identifiers.
#' @param target_var Character, target variable name.
#' @param x_field Character, X/longitude field.
#' @param y_field Character, Y/latitude field.
#' @param formula_main Formula for standard models.
#' @param formula_geo Formula for FFRK models.
#' @param train_control caret trainControl object.
#' @param split_ratio Numeric, the proportion of data to be used for training (e.g., 0.8 for 80%).
#' @param seed Integer, the random seed for data partitioning to ensure reproducibility.
#'
#' @return A list with final_results (a data.frame of detailed predictions) and
#'         metrics_summary (a data.frame of aggregated performance metrics).
#'
#' @importFrom caret createDataPartition
#' @importFrom rpart rpart
#'
#' @export
run_train_test_evaluation <- function(data,
                                      models_list,
                                      target_var,
                                      x_field,
                                      y_field,
                                      formula_main,
                                      formula_geo,
                                      train_control,
                                      split_ratio = 0.8,
                                      seed = 20250718) { # 使用与代码2一致的种子

  cat("Initializing single train-test evaluation...\n")
  hier_tree <- rpart::rpart(stats::as.formula(paste(target_var, "~", x_field, "+", y_field)),
                            data = data, method = "anova")
  data$region <- as.factor(hier_tree$where)

  cat(sprintf("Partitioning data with a %.0f/%.0f split. Seed: %d\n", split_ratio*100, (1-split_ratio)*100, seed))
  set.seed(seed)
  train_idx <- caret::createDataPartition(data[[target_var]], p = split_ratio, list = FALSE)
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]

  final_results <- data.frame()
  metrics_summary <- data.frame(model=character(), R2=numeric(), MAE=numeric(), RMSE=numeric(), stringsAsFactors=FALSE)

  pb <- progress::progress_bar$new(
    total = length(models_list),
    format = "  Evaluating [:bar] :percent | Model: :model | Time: :elapsed"
  )

  cat("Starting model training and evaluation...\n")
  for (m in models_list) {
    pb$tick(tokens = list(model = m))

    res <- run_fold(
      train_data = train_data,
      test_data = test_data,
      model_type = m,
      fold_i = 1,
      target_var = target_var,
      x_field = x_field,
      y_field = y_field,
      train_control = train_control,
      formula_main = formula_main,
      formula_geo = formula_geo
    )

    o_vals <- unlist(res$o)
    p_vals <- unlist(res$p)
    if (length(o_vals) > 0 && length(p_vals) > 0) {
      metrics <- data.frame(
        model = m,
        R2    = metric_r2(o_vals, p_vals),
        MAE   = metric_mae(o_vals, p_vals),
        RMSE  = metric_rmse(o_vals, p_vals)
      )
      metrics_summary <- rbind(metrics_summary, metrics)

      var_vals <- if(is.null(res$p_var)) rep(NA, length(o_vals)) else unlist(res$p_var)

      fold_results <- data.frame(
        model = rep(m, length(o_vals)),
        fold = rep(1, length(o_vals)),
        o = o_vals,
        p = p_vals,
        variance = var_vals,
        DLONG = res$DLONG,
        DLAT = res$DLAT,
        autofit_model = res$vario_model_str,
        autofit_psill = res$vario_psill_str,
        autofit_range = res$vario_range_str,
        autofit_kappa = res$vario_kappa_str,
        stringsAsFactors = FALSE
      )
      final_results <- rbind(final_results, fold_results)
    } else {
      cat(paste("\nWarning: Model", m, "produced no valid results. Skipping.\n"))
    }
  }

  cat("\nEvaluation complete.\n")
  return(list(detailed_results = final_results, metrics_summary = metrics_summary))
}
