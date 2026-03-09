#' Get Base ML Methods Mapping for caret
#'
#' @return A named character vector mapping method aliases to caret method names.
#' @export
get_base_methods_caret <- function() {
  c(
    lm    = "lm",
    rf    = "rf",
    rpart = "rpart",
    svmRadial = "svmRadial"
  )
}

#' Generate Model Type List Based on Base Methods
#'
#' @param methods A named vector from get_base_methods_caret().
#' @param include_grk Logical. Whether to include GRK-prefixed models.
#'
#' @return A character vector of model types (e.g., "rf", "GRK_rf")
#' @export
generate_model_list <- function(methods = get_base_methods_caret(), include_grk = TRUE) {
  base <- names(methods)
  if (include_grk) {
    return(c(base, paste0("GRK_", base)))
  }
  return(base)
}

#' Get Default caret::trainControl Object
#'
#' @return A caret trainControl object with default 5-fold CV and parallel support.
#' @export
get_default_train_control <- function() {
  trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
}
