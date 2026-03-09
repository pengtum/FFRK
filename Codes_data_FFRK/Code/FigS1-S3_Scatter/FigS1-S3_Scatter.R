rm(list = ls()); gc()

library(dplyr)
library(readr)
library(ggplot2)
library(MASS)       
library(ggpubr)     
library(readxl)    

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(workspace_dir,"Data")
out_dir <- file.path(workspace_dir, "Results", "Scatter")

target_element <- "Zn"  

old_data_path <- file.path(data_dir, paste0(target_element, "_Result.csv"))
data_old <- read_csv(old_data_path)

new_data_path <- file.path(data_dir, paste0(target_element, "_Hier_Result.csv"))
data_new <- read_csv(new_data_path) %>%
  filter(k == 5)

stk_data_path <- file.path(data_dir, paste0(target_element, "_STK_Result.csv"))
data_stk <- read_csv(stk_data_path) %>%
  filter(k == 5)

rename_model <- function(df, model_col = "model") {
  df %>%
    mutate({{model_col}} := dplyr::recode(
      .data[[model_col]],
      # ------- FFRK -------
      "FFRK_lm"         = "FFRK_LM",
      "FFRK_rpart"      = "FFRK_DT",
      "FFRK_rf"         = "FFRK_RF",
      "FFRK_svmRadial"  = "FFRK_SVM",
      # ------- RK -------
      "rk_lm"          = "GRK_LM",
      "rk_rpart"       = "GRK_DT",
      "rk_rf"          = "GRK_RF",
      "rk_svmRadial"   = "GRK_SVM",
      # ------- hier -------
      "hier_lm"        = "Stratified_LM",
      "hier_rpart"     = "Stratified_DT",
      "hier_rf"        = "Stratified_RF",
      "hier_svmRadial" = "Stratified_SVM",
      # ------- ML -------
      "lm"             = "LM",
      "rpart"          = "DT",
      "rf"             = "RF",
      "svmRadial"      = "SVM"

    ))
}

data_old <- rename_model(data_old, "model")
data_new <- rename_model(data_new, "model")
data_stk <- rename_model(data_stk, "model") 
data_stk <- data_stk %>%
  mutate(
    o = as.numeric(o),
    p = as.numeric(p)
  )


data_old <- data_old %>% mutate(source = "old")
data_new <- data_new %>% mutate(source = "new")
data_stk <- data_stk %>% mutate(source = "stk")

all_data <- bind_rows(data_old, data_new, data_stk)

new_levels <- c(
  # 1) FFRK
  "FFRK_LM", "FFRK_DT", "FFRK_RF", "FFRK_SVM",
  # 2) GRK
  "GRK_LM", "GRK_DT", "GRK_RF", "GRK_SVM",
  # 3) ML
  "LM", "DT", "RF", "SVM",
  # 4) Stratified
  "Stratified_LM", "Stratified_DT", "Stratified_RF", "Stratified_SVM",
  # 5) others
  "STK"
)

all_data <- all_data %>%
  mutate(model = factor(model, levels = new_levels))

calculate_density <- function(x, y, n = 100) {
  valid <- is.finite(x) & is.finite(y)
  
  if(sum(valid) == 0) {
    return(rep(NA, length(x)))
  }
  
  dens <- MASS::kde2d(x[valid], y[valid], n = n)
  density <- rep(NA, length(x))
  ix <- findInterval(x[valid], dens$x)
  iy <- findInterval(y[valid], dens$y)
  density[valid] <- dens$z[cbind(ix, iy)]
  
  return(density)
}


all_data <- all_data %>%
  group_by(model) %>%
  mutate(density = calculate_density(o, p, n = 100)) %>%
  ungroup()

if(all(c("o","p") %in% names(all_data))) {
  x_range <- range(all_data$o, na.rm = TRUE)
  y_range <- range(all_data$p, na.rm = TRUE)
  x_padding <- diff(x_range) * 0.01
  y_padding <- diff(y_range) * 0.01
  
  x_lim <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_lim <- c(min(0.65, y_range[1] - y_padding),
             max(-0.1, y_range[2] + y_padding))
} else {
  x_lim <- NULL
  y_lim <- NULL
}
plot <- ggplot(all_data, aes(x = o, y = p, color = density)) +
  geom_point(size = 1, alpha = 1) +
  coord_cartesian(
    xlim = if(!is.null(x_lim)) x_lim else NULL,
    ylim = if(!is.null(y_lim)) y_lim else NULL
  ) +
  scale_color_viridis_c(option = "viridis", name = "# of data") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", size = 0.55) +
  facet_wrap(~ model, ncol = 4) +
  labs(x = "Observed", y = "Predicted") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = rgb(157/255, 171/255, 200/255), color = "black"),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    text = element_text(family = "Times New Roman", size = 14)
  )


print(plot)

ggsave(
  filename = file.path(out_dir, paste0("Scatter_", target_element, ".jpg")),
  plot = plot,
  width = 9,
  height = 10,
  dpi = 600,
  create.dir = TRUE
)

