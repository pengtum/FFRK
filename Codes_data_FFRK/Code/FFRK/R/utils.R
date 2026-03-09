#' Initialize Result Storage and Progress Bar
#'
#' Creates empty result containers and a progress bar for cross-validation loops.
#'
#' @param models_list Vector of model names.
#' @param num_repeats Number of repetitions.
#' @param k_folds Number of folds.
#'
#' @return A list with elements: `final_results`, `cv_results`, and `pb` (progress bar object).
#'
#' @importFrom progress progress_bar
#'
#' @export
initialize_cv_storage <- function(models_list, num_repeats, k_folds) {
  final_results <- data.frame(
    model = character(),
    fold  = integer(),
    o     = numeric(),
    p     = numeric(),
    autofit_model = character(),
    autofit_psill = character(),
    autofit_range = character(),
    autofit_kappa = character(),
    stringsAsFactors = FALSE
  )

  cv_results <- lapply(models_list, function(x) {
    data.frame(R2 = numeric(), MAE = numeric(), RMSE = numeric())
  })
  names(cv_results) <- models_list

  total_steps <- num_repeats * length(models_list) * k_folds
  pb <- progress_bar$new(
    format = paste0(
      "  Progress [:bar] :percent ",
      "[Model: :model_name] ",
      "[Fold: :iteration/:total_fold] ",
      "[Repeat: :rep/:total_rep] ",
      "Elapsed: :elapsed Remaining: :eta"
    ),
    total = total_steps,
    clear = TRUE,
    width = 100
  )

  return(list(final_results = final_results,
              cv_results = cv_results,
              pb = pb))
}
