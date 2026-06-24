 
library(scatterplot3d)
library(reshape2)

df <- read.csv("d_n.csv", header = TRUE, check.names = FALSE)
df1 <- melt(df, id.vars = "Time", measure.vars = c("LS","NMZ","YZ","MT"),
            variable.name = "Group", value.name = "Value")


get_period <- function(t) {
 
  if ((t >= 0 & t <= 6) | (t >= 18 & t <= 23)) {
    return("Night")
  } else {
    return("Day")
  }
}

df1$Period <- sapply(df1$Time, get_period)

# ==========================================
# 1. 
# ==========================================
y_map <- c(LS=1, NMZ=2, YZ=3, MT=4)
group_colors <- c(LS="#F5A623", NMZ="#D0021B", YZ="#BD10E0", MT="#4A90E2")

# 
night_color_wall  <- adjustcolor("grey85", alpha.f = 0.6) 
night_color_floor <- adjustcolor("grey90", alpha.f = 0.6) 

# 
max_val <- max(df1$Value, na.rm=TRUE) * 1.25 # 调高一点，防止星号写不下
x_start <- 0
x_end   <- 23 
y_limit_max <- 4.5

# ==========================================
# 2. 
# ==========================================
par(mar = c(5, 5, 2, 4)) 

p <- scatterplot3d(
  x = df1$Time, 
  y = as.numeric(y_map[df1$Group]), 
  z = df1$Value, 
  type = 'n', xlim = c(x_start, x_end), ylim = c(0.5, 4.5), zlim = c(0, max_val), 
  angle = 50, scale.y = 1.3, box = FALSE, grid = FALSE, axis = FALSE, 
  xlab = "", ylab = "", zlab = "" 
)

# ==========================================
# 3. 
# ==========================================
# 3.1
draw_night_zone <- function(xs, xe) {
  x_f <- c(xs, xe, xe, xs); y_f <- c(0.5, 0.5, y_limit_max, y_limit_max); z_f <- c(0,0,0,0)
  c_f <- p$xyz.convert(x_f, y_f, z_f)
  polygon(c_f$x, c_f$y, col = night_color_floor, border = NA)
  
  x_w <- c(xs, xe, xe, xs); y_w <- rep(y_limit_max, 4); z_w <- c(0, 0, max_val, max_val)
  c_w <- p$xyz.convert(x_w, y_w, z_w)
  polygon(c_w$x, c_w$y, col = night_color_wall, border = NA)
}
draw_night_zone(0, 6)
draw_night_zone(18, 23)

# 3.2 
p$points3d(c(0, 0), c(0.5, 0.5), c(0, max_val), type='l', col="black") # Z
z_at <- pretty(c(0, max_val), n=5)
for(z in z_at) {
  p$points3d(c(0, -0.5), c(0.5, 0.5), c(z, z), type='l', col="black")
  text(p$xyz.convert(-1, 0.5, z), labels=z, pos=2, cex=0.8)
}
text(p$xyz.convert(-4, 2.5, max_val/2), "Counts of acoustic sightings", srt=90, cex=1)

p$points3d(c(0, 23), c(0.5, 0.5), c(0, 0), type='l', col="black") # X
x_at <- seq(0, 23, 4)
for(x in x_at) {
  p$points3d(c(x, x), c(0.5, 0.4), c(0, 0), type='l', col="black") 
  text(p$xyz.convert(x, 0.2, 0), labels=x, cex=0.8, pos=1)
}
text(p$xyz.convert(11.5, 0, -max_val*0.1), "Time (Hour)", pos=1)

# 3.3 
for(z in z_at[-1]) {
  p$points3d(c(0, 23), c(y_limit_max, y_limit_max), c(z, z), type='l', col="grey80", lty="dotted")
}

# ==========================================
# 4. 
# ==========================================
plot_order <- c("MT", "YZ", "NMZ", "LS") 

for (g in plot_order) {
  data_g <- df1[df1$Group == g,]
  if(nrow(data_g) < 2) next 
  
  # --- Step A: Wilcoxon Test ---
  
  groups_present <- unique(data_g$Period)
  sig_label <- "" 
  
  if (length(groups_present) == 2) {
    
    vals_day <- data_g$Value[data_g$Period == "Day"]
    vals_night <- data_g$Value[data_g$Period == "Night"]
    
    
    if(length(vals_day) > 0 && length(vals_night) > 0) {
     
      res <- try(wilcox.test(vals_day, vals_night), silent=TRUE)
      
      if (!inherits(res, "try-error")) {
        p_val <- res$p.value
        
        if(p_val < 0.001) { sig_label <- "***" }
        else if(p_val < 0.01) { sig_label <- "**" }
        else if(p_val < 0.05) { sig_label <- "*" }
        else { sig_label <- "ns" } 
      }
    }
  }
  
  # --- Step B: ---
  sp <- spline(data_g$Time, data_g$Value, n = 200)
  idx <- sp$x >= 0 & sp$x <= 23
  sp$x <- sp$x[idx]; sp$y <- sp$y[idx]
  sp$y[sp$y < 0] <- 0
  
  col_line <- group_colors[g]
  cur_y <- y_map[g]
  
  p$points3d(sp$x, rep(cur_y, length(sp$x)), sp$y, type = 'l', lwd = 3.5, col = col_line)
  p$points3d(c(0, 23), c(cur_y, cur_y), c(0, 0), type = 'l', lwd=1, col=adjustcolor(col_line, 0.2))
  
  max_idx <- which.max(sp$y)
  p$points3d(sp$x[max_idx], cur_y, sp$y[max_idx], pch=19, cex=0.6, col=col_line)
  t_pt <- p$xyz.convert(sp$x[max_idx], cur_y, sp$y[max_idx] + max_val*0.05)
  text(t_pt$x, t_pt$y, labels=round(sp$y[max_idx]), cex=0.6, col="grey30")
  
  # --- Step C:  ---
  p$points3d(c(23, 23.5), c(cur_y, cur_y), c(0, 0), type='l', col="black") 
  
  t_loc <- p$xyz.convert(24, cur_y, 0)
  
  final_label <- paste(g, sig_label) 
  text(t_loc$x, t_loc$y, labels=final_label, col=col_line, font=2, cex=0.9, pos=4)
}


legend("topleft", legend = c("Night Time", "ns: no significant","* : p<0.05", "** : p<0.01"), 
       fill = c(night_color_floor, NA, NA), border=NA, bty = "n", cex = 0.8)
