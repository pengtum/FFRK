#' Get Model Meta Information
#'
#' Returns modeling type and method name based on model identifier.
#'
#' @param model_type A string (e.g., "rf", "FFRK_rf")
#' @param formula_main Formula for base models (e.g., Pb_ppm ~ x1 + x2 + ...)
#' @param formula_geo  Formula for FFRK models (e.g., Pb_ppm ~ f.idw + f.sda.0 + ...)
#'
#' @return A list with type, caret method, and applicable formula
#'
#' @export
get_model_info <- function(model_type,
                           formula_main = NULL,
                           formula_geo = NULL) {
  base_methods <- c(lm = "lm", rf = "rf", rpart = "rpart", svmRadial = "svmRadial")

  if (model_type == "ok") {
    return(list(type = "ok", caret_method = NA, use_formula = NULL))
  } else if (model_type == "uk") {  # <--- 新增这部分
    return(list(type = "uk", caret_method = NA, use_formula = NULL))
  } else if (model_type %in% names(base_methods)) {
    return(list(type = "ml", caret_method = base_methods[model_type], use_formula = formula_main))
  } else if (grepl("^rk_", model_type)) {
    base <- sub("^rk_", "", model_type)
    return(list(type = "rk", caret_method = base_methods[base], use_formula = formula_main))
  } else if (grepl("^FFRK_", model_type)) {
    base <- sub("^FFRK_", "", model_type)
    return(list(type = "FFRK", caret_method = base_methods[base], use_formula = formula_geo))
  } else if (grepl("^hier_", model_type)) {
    base <- sub("^hier_", "", model_type)
    return(list(type = "hier", caret_method = base_methods[base], use_formula = formula_main))
  }
  stop("Unknown model type: ", model_type)
}

#' Generate Standard Formulas
#'
#' @param yvar Target variable (e.g., "Pb_ppm")
#' @param features Vector of covariates (for ML/RK/Hier)
#' @param geo_features Vector of geo-self features (for FFRK)
#'
#' @return A list with `formula_main` and `formula_geo`
#' @export
generate_formula <- function(yvar, features, geo_features) {
  f_main <- as.formula(paste(yvar, "~", paste(features, collapse = "+")))
  f_geo  <- as.formula(paste(yvar, "~", paste(geo_features, collapse = "+")))
  list(formula_main = f_main, formula_geo = f_geo)
}
