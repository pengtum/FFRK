rm(list = ls()); gc()

library(automap)
library(gstat)
library(sp)
library(caret)
library(parallel)
library(doParallel)
library(devtools)
library(ggplot2)
library(gridExtra)
library(grid)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(workspace_dir, "Data")
ffrk_dir <- file.path(workspace_dir, "Code", "FFRK")
out_dir <- file.path(workspace_dir, "Results", "Variograms")

target_elements <- c("Cu", "Pb", "Zn") 

x_field <- "DLONG"
y_field <- "DLAT"

predictor_fields <- c("Dlith","Dfault","Slope","Water","NDVI","MainRd","Road","SOC","pH")

reproducible_seed <- 20250718

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

plots_ok <- list()
plots_rk <- list()
plots_ffrk <- list()
summary_data <- data.frame() 

create_simple_variogram <- function(auto_fit_obj) {
  exp_data <- auto_fit_obj$exp_var
  line_data <- variogramLine(auto_fit_obj$var_model, maxdist = max(exp_data$dist))
  
  y_max <- max(exp_data$gamma)
  y_min <- min(exp_data$gamma)
  
  limit_setting <- c(0, y_max * 1.05)
  
  p <- ggplot() +
    geom_line(data = line_data, aes(x = dist, y = gamma), color = "black", size = 1) +
    geom_point(data = exp_data, aes(x = dist, y = gamma), color = "black", shape = 16, size = 2.5) +
    
    scale_y_continuous(limits = limit_setting, expand = c(0,0)) +
    
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(), 
      axis.text = element_text(size = 9, color = "black"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )
  return(p)
}

get_variogram_params <- function(obj, elem_name, method_name) {
  model <- obj$var_model
  nugget <- model$psill[1]
  psill <- model$psill[2]
  total_sill <- nugget + psill
  range_val <- model$range[2]
  model_type <- as.character(model$model[2])
  ratio <- nugget / total_sill
  
  return(data.frame(
    Element = elem_name,
    Method = method_name,
    Model_Type = model_type,
    Nugget = round(nugget, 5),
    Partial_Sill = round(psill, 5),
    Total_Sill = round(total_sill, 5),
    Range = round(range_val, 4),
    Nugget_Sill_Ratio = round(ratio, 4)
  ))
}

for (elem in target_elements) {
  cat(paste0(">>> Processing: ", elem, "\n"))
  
  data_file <- paste0("data_", elem, ".csv") 
  target_var <- paste0(elem, "_ppm")
  
  data <- prepare_data(data_file, value_field = target_var, x_field = x_field, y_field = y_field, log_transform = TRUE)
  
  features <- getfeature(xo = data[[x_field]], yo = data[[y_field]], zo = data[[target_var]], 
                         xp = data[[x_field]], yp = data[[y_field]], k_neighbors = 15, num_quantiles = 11)
  
  data_model <- cbind(data, features)
  data_model <- na.omit(data_model)
  
  formulas <- generate_formula(target_var, predictor_fields, names(features))
  
  data_sp <- data_model
  coordinates(data_sp) <- as.formula(paste("~", x_field, "+", y_field))
  
  f_ok <- as.formula(paste(target_var, "~ 1"))
  vario_ok <- autofitVariogram(f_ok, data_sp, model = c("Sph", "Exp"))
  
  plots_ok[[elem]] <- create_simple_variogram(vario_ok)
  summary_data <- rbind(summary_data, get_variogram_params(vario_ok, elem, "OK (Raw Data)"))
  
  set.seed(reproducible_seed)
  fit_rk <- train(formulas$formula_main, data = data_model, method = "rf", 
                  trControl = trainControl(method = "none"), 
                  tuneGrid = data.frame(mtry = floor(length(predictor_fields)/3)), na.action = na.omit)
  
  res_rk <- data_model[[target_var]] - predict(fit_rk, data_model)
  data_sp_rk <- data_sp
  data_sp_rk$residuals <- res_rk
  
  vario_rk <- autofitVariogram(residuals ~ 1, data_sp_rk, model = c("Sph", "Exp"))
  
  plots_rk[[elem]] <- create_simple_variogram(vario_rk)
  summary_data <- rbind(summary_data, get_variogram_params(vario_rk, elem, "RK-RF (Residuals)"))
  
  set.seed(reproducible_seed)
  fit_ffrk <- train(formulas$formula_geo, data = data_model, method = "rf", 
                    trControl = trainControl(method = "none"), 
                    tuneGrid = data.frame(mtry = floor(ncol(features)/3)), na.action = na.omit)
  
  res_ffrk <- data_model[[target_var]] - predict(fit_ffrk, data_model)
  data_sp_ffrk <- data_sp
  data_sp_ffrk$residuals <- res_ffrk
  
  vario_ffrk <- autofitVariogram(residuals ~ 1, data_sp_ffrk, model = c("Sph", "Exp"))
  
  plots_ffrk[[elem]] <- create_simple_variogram(vario_ffrk)
  summary_data <- rbind(summary_data, get_variogram_params(vario_ffrk, elem, "FFRK-RF (Residuals)"))
}

stopCluster(cl)
registerDoSEQ()

header_ok <- textGrob("OK", gp = gpar(fontsize = 12, fontface = "bold"))
header_rk <- textGrob("RK (RF)", gp = gpar(fontsize = 12, fontface = "bold"))
header_ffrk <- textGrob("FFRK (RF)", gp = gpar(fontsize = 12, fontface = "bold"))

label_cu <- textGrob("Cu", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
label_pb <- textGrob("Pb", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
label_zn <- textGrob("Zn", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))

common_y_label <- textGrob("Semivariance", rot = 90, gp = gpar(fontsize = 14))
common_x_label <- textGrob("Distance (degrees)", gp = gpar(fontsize = 14))

widths_config <- c(1, 10, 10, 10)

row_cu <- arrangeGrob(label_cu, plots_ok[["Cu"]], plots_rk[["Cu"]], plots_ffrk[["Cu"]], ncol = 4, widths = widths_config)
row_pb <- arrangeGrob(label_pb, plots_ok[["Pb"]], plots_rk[["Pb"]], plots_ffrk[["Pb"]], ncol = 4, widths = widths_config)
row_zn <- arrangeGrob(label_zn, plots_ok[["Zn"]], plots_rk[["Zn"]], plots_ffrk[["Zn"]], ncol = 4, widths = widths_config)

header_row <- arrangeGrob(nullGrob(), header_ok, header_rk, header_ffrk, ncol = 4, widths = widths_config)

main_grid <- arrangeGrob(
  header_row,
  row_cu,
  row_pb,
  row_zn,
  nrow = 4,
  heights = c(1, 10, 10, 10)
)

final_plot <- arrangeGrob(main_grid, left = common_y_label, bottom = common_x_label)

jpeg(file.path(out_dir, "Variograms.jpg"), 
     width = 16,
     height = 9.33,
     units = "in",
     res = 600,
     quality = 100
)

grid.draw(final_plot)

dev.off()


summary_data <- summary_data[order(summary_data$Element, summary_data$Method), ]


print(head(summary_data, 10))

write.csv(summary_data, file.path(out_dir, "Variogram_Parameters.csv"), row.names = FALSE)