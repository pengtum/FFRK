library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

script_dir <- dirname(this.path::this.path())
workspace_dir <- normalizePath(file.path(script_dir, "..", ".."))
file_path <- file.path(workspace_dir, "Data", "Global", "Global.xls")

data_zn <- read_excel(file_path, sheet = 1)
data_cu <- read_excel(file_path, sheet = 2)
data_pb <- read_excel(file_path, sheet = 3)

base_theme <- theme_minimal() +
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

make_tile_plot <- function(df) {
  cols_pred <- names(df)[4:22]
  
  df_long <- df %>%
    mutate(across(all_of(cols_pred), ~ suppressWarnings(as.numeric(.)))) %>%
    pivot_longer(
      cols = all_of(cols_pred),
      names_to = "method",
      values_to = "pred"
    ) %>%
    mutate(method = factor(method, levels = cols_pred))
  
  ggplot(df_long, aes(x = DLONG, y = DLAT, fill = pred)) +
    geom_tile() +
    facet_wrap(~ method, ncol = 4) +
    scale_fill_viridis_c(option = "viridis", na.value = "grey70") +
    base_theme
}

p_zn <- make_tile_plot(data_zn)
p_cu <- make_tile_plot(data_cu)
p_pb <- make_tile_plot(data_pb)

p_zn
p_cu
p_pb
