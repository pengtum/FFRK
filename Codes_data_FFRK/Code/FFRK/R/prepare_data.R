#' Prepare Input Data for Spatial Modeling
#'
#' This function reads a CSV file, logs the target variable (if needed),
#' handles missing values, and ensures spatial coordinate fields are numeric.
#'
#' @param filepath A character string. Path to the CSV file to be loaded.
#' @param value_field A character string. Column name of the target variable (e.g., "Pb_ppm").
#' @param x_field A character string. Column name of the longitude or X coordinate (e.g., "DLONG").
#' @param y_field A character string. Column name of the latitude or Y coordinate (e.g., "DLAT").
#' @param log_transform Logical. Whether to log-transform the target variable. Default is TRUE.
#'
#' @return A cleaned data.frame with log-transformed target (if specified) and numeric coordinates.
#'
#' @examples
#' \dontrun{
#' data <- prepare_data("data_Pb.csv", value_field = "Pb_ppm",
#'                      x_field = "DLONG", y_field = "DLAT")
#' }
#'
#' @export
prepare_data <- function(filepath,
                         value_field,
                         x_field,
                         y_field,
                         log_transform = TRUE) {

  df <- read.csv(filepath)

  # Check columns
  if (!(value_field %in% names(df))) stop("Target value field not found.")
  if (!(x_field %in% names(df))) stop("X coordinate field not found.")
  if (!(y_field %in% names(df))) stop("Y coordinate field not found.")



  # Log-transform the target variable (if needed)
  if (log_transform) {
    df[[value_field]] <- log(df[[value_field]])
  }

  # Drop rows with missing values
  df <- na.omit(df)

  # Convert coordinates to numeric
  df[[x_field]] <- as.numeric(df[[x_field]])
  df[[y_field]] <- as.numeric(df[[y_field]])

  return(df)
}
