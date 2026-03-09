#' Extract Geo-self Features for Spatial Prediction
#'
#' This function generates spatial features (IDW, quantiles, and similarity-based indicators)
#' for each prediction location based on observed x, y, z values. It includes:
#' - Inverse Distance Weighted (IDW) prediction
#' - Local spatial distributional attributes (21 quantiles)
#' - Similarity-based feature using geosimilarity
#'
#' @param xo A numeric vector of longitude (or x) coordinates of observation points.
#' @param yo A numeric vector of latitude (or y) coordinates of observation points.
#' @param zo A numeric vector of observed values (e.g., log(Pb)).
#' @param xp A numeric vector of longitude coordinates of prediction points.
#' @param yp A numeric vector of latitude coordinates of prediction points.
#' @param k_neighbors Integer. The number of nearest neighbors to use for quantile and similarity features. Default is 15.
#' @param num_quantiles Integer. The number of quantiles to generate. Default is 11.
#'
#' @return A data.frame with extracted features:
#' - `f.idw`: IDW predicted value
#' - `f.sda.0` to `f.sda.1`: 21 quantile-based features
#' - `f.gos`: geosimilarity-based similarity prediction
#'
#' @examples
#' # Assuming xo, yo, zo, xp, yp are defined
#' # features <- getfeature(xo, yo, zo, xp, yp, k_neighbors = 15, num_quantiles = 11)
#'
#' @import sp
#' @importFrom gstat idw
#' @import geosimilarity
#'
#' @export
getfeature <- function(xo, yo, zo, xp, yp, k_neighbors = 15, num_quantiles = 11){
  m <- length(zo)
  n <- length(xp)

  quantile_probs <- seq(0, 1, length.out = num_quantiles)

  coordinates <- data.frame("xi" = xo, "yi" = yo)
  coords_obs <- coordinates(coordinates)
  values <- zo
  obs <- SpatialPointsDataFrame(coords = coordinates, data.frame(values))

  # 1. IDW prediction
  f.idw <- numeric(n)
  for (i in 1:n){
    pred.i <- data.frame("xi" = xp[i], "yi" = yp[i])
    coordinates(pred.i) <- ~xi + yi
    non_matching_indices <- which(!(coords_obs[, "xi"] == pred.i$xi &
                                      coords_obs[, "yi"] == pred.i$yi))
    filtered_obs <- obs[non_matching_indices, ]

    idw.result <- idw(formula = values ~ 1, locations = filtered_obs, newdata = pred.i)
    f.idw[i] <- idw.result$var1.pred
  }

  # 2. Quantile-based features
  f.sda <- data.frame(matrix(NA, nrow = n, ncol = num_quantiles))
  names(f.sda) <- paste("f.sda", round(quantile_probs, 4), sep = ".")
  for (i in 1:n){
    distances <- sqrt((xo - xp[i])^2 + (yo - yp[i])^2)
    non_zero_indices <- which(distances > 0)
    sorted_indices <- non_zero_indices[order(distances[non_zero_indices])]
    k_top <- sorted_indices[1:k_neighbors]
    f.sda[i,] <- quantile(zo[k_top], probs = quantile_probs, na.rm = TRUE)
  }

  # 3. Quantiles for observations
  f.sda.obs <- data.frame(matrix(NA, nrow = m, ncol = num_quantiles))
  names(f.sda.obs) <- paste("f.sda", round(quantile_probs, 4), sep = ".")
  for (i in 1:m){
    distances <- sqrt((xo - xo[i])^2 + (yo - yo[i])^2)
    non_zero_indices <- which(distances > 0)
    sorted_indices <- non_zero_indices[order(distances[non_zero_indices])]
    k_top <- sorted_indices[1:k_neighbors]
    f.sda.obs[i,] <- quantile(zo[k_top], probs = quantile_probs, na.rm = TRUE)
  }

  # 4. Similarity-based feature (Geo-Similarity)
  f.gos <- numeric(n)
  d3 <- data.frame(zo, f.sda.obs)
  for (i in 1:n){
    pred.i <- data.frame("xi" = xp[i], "yi" = yp[i])
    non_matching_indices <- which(!(coords_obs[, "xi"] == pred.i$xi &
                                      coords_obs[, "yi"] == pred.i$yi))
    filtered_d3 <- d3[non_matching_indices, ]
    predictor_columns <- names(filtered_d3)[-1]
    formula_str <- paste("zo ~", paste(predictor_columns, collapse = " + "))
    formula <- as.formula(formula_str)

    d4 <- f.sda[i, , drop = FALSE]
    g2 <- gos(formula, data = filtered_d3, newdata = d4, kappa = 0.05)
    f.gos[i] <- g2$pred
  }

  out <- data.frame(f.idw, f.sda, f.gos)
  return(out)
}
