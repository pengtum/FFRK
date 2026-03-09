#' Calculate Coefficient of Determination (R-squared)
#'
#' Computes the proportion of variance explained by predictions.
#'
#' @param observed A numeric vector of observed values.
#' @param predicted A numeric vector of predicted values.
#'
#' @return A numeric value representing R².
#'
#' @examples
#' metric_r2(c(1,2,3), c(1.1,1.9,3.2))
#'
#' @export
metric_r2 <- function(observed, predicted) {
  valid <- !is.na(observed) & !is.na(predicted)
  o <- observed[valid]
  p <- predicted[valid]
  if(length(o) < 2) return(NA)
  1 - sum((o - p)^2) / sum((o - mean(o))^2)
}

#' Calculate Mean Absolute Error (MAE)
#'
#' Computes the mean of absolute differences between observed and predicted values.
#'
#' @param observed A numeric vector of observed values.
#' @param predicted A numeric vector of predicted values.
#'
#' @return A numeric value representing MAE.
#'
#' @examples
#' metric_mae(c(1,2,3), c(1.1,1.9,3.2))
#'
#' @export
metric_mae <- function(observed, predicted) {
  valid <- !is.na(observed) & !is.na(predicted)
  mean(abs(observed[valid] - predicted[valid]), na.rm = TRUE)
}

#' Calculate Root Mean Square Error (RMSE)
#'
#' Computes the root mean of squared differences between observed and predicted values.
#'
#' @param observed A numeric vector of observed values.
#' @param predicted A numeric vector of predicted values.
#'
#' @return A numeric value representing RMSE.
#'
#' @examples
#' metric_rmse(c(1,2,3), c(1.1,1.9,3.2))
#'
#' @export
metric_rmse <- function(observed, predicted) {
  valid <- !is.na(observed) & !is.na(predicted)
  sqrt(mean((observed[valid] - predicted[valid])^2, na.rm = TRUE))
}

#' Calculate Mean Evaluation Metrics
#'
#' This function calculates the mean R², MAE, and RMSE from a data.frame
#' containing these metrics (e.g., `cv_results[[m]]`).
#'
#' @param df A data.frame with columns `R2`, `MAE`, and `RMSE`.
#'
#' @return A named numeric vector with mean values for R², MAE, and RMSE.
#'
#' @examples
#' \dontrun{
#' calc_metrics(data.frame(R2 = c(0.8, 0.9),
#'                         MAE = c(1.2, 1.0),
#'                         RMSE = c(2.0, 1.8)))
#' }
#'
#' @export
calc_metrics <- function(df) {
  c(
    R2_mean   = mean(df$R2,   na.rm = TRUE),
    MAE_mean  = mean(df$MAE,  na.rm = TRUE),
    RMSE_mean = mean(df$RMSE, na.rm = TRUE)
  )
}

