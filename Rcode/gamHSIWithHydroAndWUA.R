library(mgcv)
library(terra)
library(readr)


df_data <- read_csv("hydro_with_porpoise.csv")

df_data$habitat_type <- "run"  
df_data$habitat_type[df_data$depth > 15 & df_data$speed < 0.5] <- "deep_pool"
df_data$habitat_type[df_data$speed > 1.0] <- "riffle"
df_data$habitat_type <- factor(df_data$habitat_type, 
                               levels = c("run", "deep_pool", "riffle"))

table(df_data$habitat_type)

# GAM
gam_global <- gam(
  porpoise_present ~ 
    s(depth, k = 5) + 
    s(speed, k = 5) + 
    ti(depth, speed) +
    s(vorticity, k = 5) +
    habitat_type,
  data = df_data,
  family = binomial,
  method = "REML"
)

summary(gam_global)

library(mgcv)
library(ggplot2)
library(gratia)
library(patchwork)
library(dplyr)
library(scico) 


theme_elegant <- function() {
  theme_minimal(base_size = 13, base_family = "sans") %+replace%
    theme(
      
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      
      axis.line = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks = element_line(color = "grey30", linewidth = 0.5),
      axis.text = element_text(color = "grey30", size = 11),
      axis.title = element_text(color = "black", size = 13, face = "bold"),
      
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 15, hjust = 0, color = "grey10"),
      plot.subtitle = element_text(size = 12, color = "grey50", hjust = 0),
      
      plot.margin = margin(15, 15, 15, 15)
    )
}


col_depth <- "#4A90E2"  
col_speed <- "#E76F51"  
col_vort  <- "#2A9D8F"  

sm_df <- smooth_estimates(gam_global) |> add_confint()

# (1) s(depth)
p_depth <- sm_df |>
  filter(.smooth == "s(depth)") |> 
  ggplot(aes(x = depth, y = .estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = col_depth, alpha = 0.2) + 
  geom_line(color = col_depth, linewidth = 1.2) +
  labs(x = "Depth (m)", y = "Partial effect", title = "a. Effect of Depth") +
  theme_elegant()

# (2) s(speed)
p_speed <- sm_df |>
  filter(.smooth == "s(speed)") |>
  ggplot(aes(x = speed, y = .estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = col_speed, alpha = 0.2) +
  geom_line(color = col_speed, linewidth = 1.2) +
  labs(x = "Speed (m/s)", y = "Partial effect", title = "b. Effect of Speed") +
  theme_elegant()

# (3)  s(vorticity)
p_vort <- sm_df |>
  filter(.smooth == "s(vorticity)") |>
  ggplot(aes(x = vorticity, y = .estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = col_vort, alpha = 0.2) +
  geom_line(color = col_vort, linewidth = 1.2) +
  labs(x = expression(paste("Vorticity (s"^{-1}, ")")), y = "Partial effect", title = "c. Effect of Vorticity") +
  theme_elegant()

# (4)  ti(depth, speed)
p_inter <- sm_df |>
  filter(.smooth == "ti(depth,speed)") |>
  ggplot(aes(x = depth, y = speed, fill = .estimate)) + # 注意 fill = .estimate
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = .estimate), color = "white", alpha = 0.5, linewidth = 0.4) +
  scale_fill_scico(palette = "bamako", direction = 1, name = "Effect") + 
  labs(
    x = "Depth (m)", 
    y = "Speed (m/s)", 
    title = "d. Depth × Speed Interaction"
  ) +
  theme_elegant() +
  theme(panel.grid.major = element_blank()) 

# (5) core habitat defination
coefs <- summary(gam_global)$p.table
param_df <- data.frame(
  term = c("Run (Ref.)", "Deep Pool", "Riffle"),
  est = c(0, coefs["habitat_typedeep_pool", "Estimate"], coefs["habitat_typeriffle", "Estimate"]),
  se = c(0, coefs["habitat_typedeep_pool", "Std. Error"], coefs["habitat_typeriffle", "Std. Error"])
)

p_habitat <- ggplot(param_df, aes(x = est, y = reorder(term, est))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = est - 1.96*se, xmax = est + 1.96*se), 
                 height = 0, linewidth = 1, color = "#6C7A89") +
  geom_point(size = 4, shape = 21, fill = "white", color = "#6C7A89", stroke = 1.5) +
  labs(
    x = "Coefficient Estimate (Log-odds)", 
    y = NULL, 
    title = "e. Habitat Type Effects"
  ) +
  theme_elegant()


final_plot <- (p_depth | p_speed | p_vort) / 
  (p_inter | p_habitat) +
  plot_layout(heights = c(1, 1.2)) + 
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", family = "sans", margin = margin(b=15))
    )
  )

print(final_plot)


# 3. with environmental Layers
final_tif <- rast("6n_Environmental_Variables_Full.tif")

# 4. define habitat_type
pred_data <- as.data.frame(final_tif, xy = TRUE, na.rm = FALSE)

pred_data$habitat_type <- "run"
pred_data$habitat_type[pred_data$depth > 15 & pred_data$speed < 0.5] <- "deep_pool"
pred_data$habitat_type[pred_data$speed > 1.0] <- "riffle"
pred_data$habitat_type <- factor(pred_data$habitat_type, 
                                 levels = c( "run", "deep_pool", "riffle"))

# 5.  HSI
pred_data$hsi <- predict(gam_global, newdata = pred_data, type = "response")

hsi_raw <- rast(pred_data[, c("x", "y", "hsi")], 
                type = "xyz", 
                crs = crs(final_tif))

hsi_raw <- mask(hsi_raw, final_tif$depth)

hsi_smooth <- focal(hsi_raw, w = 5, fun = "mean", na.rm = TRUE)

sandbar_data <- read.csv("sandbar_coords.csv")
sandbar_data$poly_id <- cumsum(c(0, sandbar_data$flag[-nrow(sandbar_data)] == 0)) + 1
geom_matrix <- as.matrix(sandbar_data[, c("poly_id", "x", "y")])
sandbar <- vect(geom_matrix, type = "polygons", crs = crs(final_tif))
hsi_final <- mask(hsi_smooth, sandbar, inverse = TRUE)

plot(hsi_final, 
     col = colorRampPalette(c("lightgray", "cyan", "yellow", "red"))(100),
     main = "HSI (2024-6-20 20:00)")

writeRaster(hsi_final, "10d.tif", overwrite = TRUE)


cell_area_map <- cellSize(hsi_final, unit = "km")
wua_raster <- hsi_final * cell_area_map
total_wua_df <- global(wua_raster, fun = "sum", na.rm = TRUE)
WUA_value <- total_wua_df[1, 1]

print(paste("WUA:", round(WUA_value, 2), "km^2"))
print(paste("ha:", round(WUA_value * 100, 2), "ha"))

# ================================
# 4. WUA plot
# ================================

wua_high <- wua_raster
wua_mid  <- wua_raster
wua_low  <- wua_raster

values(wua_high)[values(hsi_final) < 0.6] <- NA
values(wua_mid)[values(hsi_final) >= 0.6 | values(hsi_final) < 0.3] <- NA
values(wua_low)[values(hsi_final) >= 0.3] <- NA

wua_high_val <- global(wua_high, "sum", na.rm = TRUE)[1,1]
wua_mid_val  <- global(wua_mid,  "sum", na.rm = TRUE)[1,1]
wua_low_val  <- global(wua_low,  "sum", na.rm = TRUE)[1,1]

wua_df <- data.frame(
  Category = factor(
    c("HSI > 0.6", "0.3 < HSI ≤ 0.6", "HSI ≤ 0.3"),
    levels = c("HSI > 0.6", "0.3 < HSI ≤ 0.6", "HSI ≤ 0.3")
  ),
  WUA_km2 = c(wua_high_val, wua_mid_val, wua_low_val),
  Color = c("#2E8B57", "#C8B400", "#D4826A")
)

par(mar = c(5, 5.5, 3, 2), family = "sans", bg = "white")

bp <- barplot(
  wua_df$WUA_km2,
  names.arg = wua_df$Category,
  col       = wua_df$Color,
  border    = NA,
  ylim      = c(0, max(wua_df$WUA_km2) * 1.25),
  ylab      = expression("WUA (km"^2*")"),
  xlab      = "Habitat Suitability Category",
  main      = "Weighted Usable Area by HSI Class\n2024-10-10 8:00",
  cex.axis  = 1.0,
  cex.lab   = 1.1,
  cex.main  = 1.2,
  las       = 1,
  yaxt      = "n"
)

axis(2, las = 1, cex.axis = 1.0)

text(
  x = bp,
  y = wua_df$WUA_km2,
  labels = sprintf("%.2f", wua_df$WUA_km2),
  pos = 3,
  cex = 0.9,
  font = 2
)

abline(h = pretty(c(0, max(wua_df$WUA_km2))), 
       col = "grey85", lty = "dashed", lwd = 0.8)

box(bty = "l")



