library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(scales)

full_data <- do.call(rbind, all_results_list)
kriging_data <- full_data %>%
  filter(grepl("^(ok|uk|rk_|FFRK_)", model)) %>%
  filter(!is.na(variance))

variance_stats <- kriging_data %>%
  group_by(Element, model) %>%
  summarise(
    Min = min(variance, na.rm = TRUE),
    Max = max(variance, na.rm = TRUE),
    Mean = mean(variance, na.rm = TRUE),
    Diff = max(variance, na.rm = TRUE) - min(variance, na.rm = TRUE), # 极差
    
    Status = ifelse(Diff < 0.0001, "×", "✔"),
    .groups = 'drop'
  )

print(as.data.frame(variance_stats))

write.csv(variance_stats, "All_Kriging_Variance_Stats.csv", row.names = FALSE)
target_models <- c("ok", "rk_rf", "FFRK_rf")
plot_data <- full_data %>% 
  filter(model %in% target_models) %>%
  mutate(
    model_label = factor(model, levels = target_models, labels = c("OK", "RK (RF)", "FFRK (RF)")),
    Element = factor(Element, levels = c("Cu", "Pb", "Zn"))
  )

common_palette <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")


create_single_panel <- function(df_input, elem_name, mod_label) {
  
  sub_df <- df_input %>% filter(Element == elem_name, model_label == mod_label)
  sub_df$PanelTitle <- mod_label
  
  v_min <- min(sub_df$variance, na.rm = TRUE)
  v_max <- max(sub_df$variance, na.rm = TRUE)
  v_diff <- v_max - v_min
  
  p <- ggplot(sub_df, aes(x = DLONG, y = DLAT)) +
    coord_fixed() + 
    theme_bw() +
    facet_wrap(~ PanelTitle) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold", margin = margin(b = 5)),
      axis.title = element_blank(), 
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    ) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = NULL)
  
  if (v_diff < 0.00001) {
    label_text <- sprintf("%.4f", v_min)
    
    p <- p + 
      geom_point(aes(color = "Constant"), size = 1.5, alpha = 1) +
      scale_color_manual(
        name = "Variance",
        values = c("Constant" = "#2c7bb6"), 
        labels = label_text
      ) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
  } else {
    plot_limits <- c(v_min, v_max)
    plot_breaks <- seq(v_min, v_max, length.out = 5)
    
    p <- p + 
      geom_point(aes(color = variance), size = 1.5, alpha = 1) +
      scale_color_gradientn(
        colors = common_palette,
        limits = plot_limits,
        breaks = plot_breaks,
        labels = scales::number_format(accuracy = 0.0001), 
        name = "Variance"
      ) +
      theme(
        legend.key.height = unit(0.6, "cm"), 
        legend.key.width = unit(0.2, "cm")
      )
  }
  
  return(p)
}

rows_list <- list()

for (e in c("Cu", "Pb", "Zn")) {
  
  p1 <- create_single_panel(plot_data, e, "OK")
  p2 <- create_single_panel(plot_data, e, "RK (RF)")
  p3 <- create_single_panel(plot_data, e, "FFRK (RF)")
  
  row_label <- wrap_elements(panel = textGrob(e, rot = 90, gp = gpar(fontsize = 14, fontface = "bold")))
  
  row_combined <- row_label + p1 + p2 + p3 + plot_layout(widths = c(1, 8, 8, 8))
  
  rows_list[[e]] <- row_combined
}

final_grid <- (rows_list[["Cu"]] / rows_list[["Pb"]] / rows_list[["Zn"]])

jpeg("Variance.jpg", 
     width = 8.5,        
     height = 7,     
     units = "in",      
     res = 600,         
     quality = 100      
)

print(final_grid)

dev.off() 
