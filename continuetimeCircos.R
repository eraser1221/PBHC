library(circlize)
library(tidyverse)

data <- read.csv("contine_time.csv")

df_long <- data %>%
  pivot_longer(cols = everything(), names_to = "Site", values_to = "Seconds")

df_processed <- df_long %>%
  mutate(Minutes = Seconds / 60) %>%
  mutate(Time_Bin = cut(
    Minutes, 
    breaks = c(0, 4, 8, 12, Inf), 
    labels = c("0-3 min", "4-7 min", "8-11 min", ">12 min"), 
    right = FALSE 
  )) %>%
  filter(!is.na(Time_Bin))

plot_data <- df_processed %>%
  group_by(Time_Bin, Site) %>%
  summarise(Value = sum(Minutes), .groups = 'drop')

adj_matrix <- plot_data %>%
  pivot_wider(names_from = Site, values_from = Value, values_fill = 0) %>%
  column_to_rownames("Time_Bin") %>%
  as.matrix()

sites_order <- c("LS", "NMZ", "YZ", "MT") 
times_order <- c("0-3 min", "4-7 min", "8-11 min", ">12 min")

full_order <- c(sites_order, times_order)

time_colors <- c(
  "0-3 min"  = "#A6CEE3",  
  "4-7 min"  = "#B2DF8A",  
  "8-11 min" = "#CAB2D6",  
  ">12 min"  = "#FDBF6F"   
)

site_colors <- c(
  "LS"  = "#66C2A5", 
  "NMZ" = "#FC8D62", 
  "YZ"  = "#8DA0CB", 
  "MT"  = "#E78AC3"
)
grid_col <- c(time_colors, site_colors)

circos.clear()
par(mar = c(1, 1, 1, 1))

gaps <- c(rep(2, 3), 15, rep(2, 3), 15)

circos.par(gap.after = gaps, start.degree = 90) 

chordDiagram(
  adj_matrix,
  order = full_order,         
  grid.col = grid_col,        
  transparency = 0.3,
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(track.height = 0.1)
)

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  
  circos.text(
    mean(xlim), 
    ylim[1] + 0.3, 
    sector.name, 
    facing = "clockwise", 
    niceFacing = TRUE, 
    adj = c(0, 0.5),
    cex = 0.8
  )
}, bg.border = NA)

circos.clear()

