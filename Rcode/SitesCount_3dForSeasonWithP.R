
library(scatterplot3d)
library(reshape2)

df <- read.csv("d_s.csv", header = TRUE, check.names = FALSE)
df1 <- melt(df, id.vars = "Month", measure.vars = c("LS","NMZ","YZ","MT"),
            variable.name = "Group", value.name = "Value")

get_water_period <- function(m) {
  if (m >= 4 & m <= 8) {
    return("Rising")   
  } else if (m >= 9 & m <= 11) {
    return("Falling")  
  } else {
    return("Flat")    
  }
}

df1$Period <- sapply(df1$Month, get_water_period)

df1$Period <- factor(df1$Period, levels = c("Flat", "Rising", "Falling"))

# ==========================================
# 1. 
# ==========================================
y_map <- c(LS=1, NMZ=2, YZ=3, MT=4)
group_colors <- c(LS="#F5A623", NMZ="#D0021B", YZ="#BD10E0", MT="#4A90E2")

col_rising_floor <- adjustcolor("#87CEFA", alpha.f = 0.3) 
col_rising_wall  <- adjustcolor("#87CEFA", alpha.f = 0.3)

col_falling_floor <- adjustcolor("#FFD700", alpha.f = 0.2)
col_falling_wall  <- adjustcolor("#FFD700", alpha.f = 0.2)


max_val <- max(df1$Value, na.rm=TRUE) * 1.3 
x_start <- 1
x_end   <- 12
y_limit_max <- 4.5

# ==========================================
# 2. 
# ==========================================
par(mar = c(5, 5, 2, 5)) 

p <- scatterplot3d(
  x = df1$Month, 
  y = as.numeric(y_map[df1$Group]), 
  z = df1$Value, 
  type = 'n',              
  xlim = c(x_start, x_end), ylim = c(0.5, 4.5), zlim = c(0, max_val), 
  angle = 50, scale.y = 1.3,           
  box = FALSE, grid = FALSE, axis = FALSE,            
  xlab = "", ylab = "", zlab = "" 
)

# ==========================================
# 3. 
# ==========================================


draw_zone <- function(xs, xe, col_f, col_w) {

  x_f <- c(xs, xe, xe, xs); y_f <- c(0.5, 0.5, y_limit_max, y_limit_max); z_f <- c(0,0,0,0)
  c_f <- p$xyz.convert(x_f, y_f, z_f)
  polygon(c_f$x, c_f$y, col = col_f, border = NA)
  

  x_w <- c(xs, xe, xe, xs); y_w <- rep(y_limit_max, 4); z_w <- c(0, 0, max_val, max_val)
  c_w <- p$xyz.convert(x_w, y_w, z_w)
  polygon(c_w$x, c_w$y, col = col_w, border = NA)
}

# 1. Rising
draw_zone(4, 8, col_rising_floor, col_rising_wall)

# 2. Falling
draw_zone(9, 11, col_falling_floor, col_falling_wall)


# ==========================================
# 4. 
# ==========================================
# Z
p$points3d(c(x_start, x_start), c(0.5, 0.5), c(0, max_val), type='l', col="black")
z_at <- pretty(c(0, max_val), n=5)
for(z in z_at) {
  p$points3d(c(x_start, x_start-0.3), c(0.5, 0.5), c(z, z), type='l', col="black")
  text(p$xyz.convert(x_start-0.5, 0.5, z), labels=z, pos=2, cex=0.8)
}
text(p$xyz.convert(x_start-2, 2.5, max_val/2), "Counts of acoustic sightings", srt=90, cex=1)

# X
p$points3d(c(x_start, x_end), c(0.5, 0.5), c(0, 0), type='l', col="black", lwd=1.5)
x_at <- c(1, 4, 8, 12)
for(x in x_at) {
  p$points3d(c(x, x), c(0.5, 0.4), c(0, 0), type='l', col="black") 
  text(p$xyz.convert(x, 0.2, 0), labels=x, cex=0.8, pos=1)
}
text(p$xyz.convert((x_start+x_end)/2, 0, -max_val*0.1), "Month", pos=1)

p$points3d(c(x_end, x_end), c(0.5, 4.5), c(0, 0), type='l', col="black")

for(z in z_at[-1]) {
  p$points3d(c(x_start, x_end), c(y_limit_max, y_limit_max), c(z, z), type='l', col="grey80", lty="dotted")
}

# ==========================================
# 5. 
# ==========================================
plot_order <- c("MT", "YZ", "NMZ", "LS") 

for (g in plot_order) {
  data_g <- df1[df1$Group == g,]
  if(nrow(data_g) < 3) next 
  
  # --- Step A: Kruskal-Wallis ---
  
  if(length(unique(data_g$Period)) > 1) {
    res <- try(kruskal.test(Value ~ Period, data = data_g), silent=TRUE)
    
    sig_label <- "" 
    if (!inherits(res, "try-error")) {
      p_val <- res$p.value
      
      if(p_val < 0.001) { sig_label <- "***" }
      else if(p_val < 0.01) { sig_label <- "**" }
      else if(p_val < 0.05) { sig_label <- "*" }
      else { sig_label <- "ns" }
    }
  } else {
    sig_label <- "" 
  }
  
  # --- Step B: ---
  sp <- spline(data_g$Month, data_g$Value, n = 200)
  idx <- sp$x >= x_start & sp$x <= x_end
  sp$x <- sp$x[idx]; sp$y <- sp$y[idx]
  if(length(sp$x)==0) next
  sp$y[sp$y < 0] <- 0
  
  cur_y <- y_map[g]
  col_line <- group_colors[g]
  
  
  p$points3d(sp$x, rep(cur_y, length(sp$x)), sp$y, type = 'l', lwd = 3.5, col = col_line)
 
  p$points3d(c(x_start, x_end), c(cur_y, cur_y), c(0, 0), type = 'l', lwd=1, col=adjustcolor(col_line, 0.2))
  
  max_idx <- which.max(sp$y)
  p$points3d(sp$x[max_idx], cur_y, sp$y[max_idx], pch=19, cex=0.6, col=col_line)
  txt_pt <- p$xyz.convert(sp$x[max_idx], cur_y, sp$y[max_idx] + max_val*0.05)
  text(txt_pt$x, txt_pt$y, labels=round(sp$y[max_idx]), cex=0.6, col="grey30")
  
  # --- Step C:  ---
  
  p$points3d(c(x_end, x_end+0.2), c(cur_y, cur_y), c(0, 0), type='l', col="black")
  
  label_text <- paste(g, sig_label)
  t_loc <- p$xyz.convert(x_end + 0.5, cur_y, 0)
  text(t_loc$x, t_loc$y, labels=label_text, col=col_line, font=2, cex=0.9, pos=4)
}

# ==========================================
# 6. 
# ==========================================
legend("topleft", 
       legend = c("Rising (4-8)", "Falling (9-11)", 
                  "Flat (12-3)", 
                  
                 "ns: no significant", "* : p<0.05", "** : p<0.01"), 
       fill = c(col_rising_floor, col_falling_floor, "white", NA, NA, NA), 
       border = c(NA, NA, "grey80", NA, NA, NA), 
       bty = "n", cex = 0.8, inset=0.02)

