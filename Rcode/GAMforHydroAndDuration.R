library(mgcv) 
library(tidyverse) 
library(gratia) 

df_clean <- read.csv("Matched_Porpoise_Hydro_Data.csv") 

if(!"log_duration" %in% names(df_clean)) { 
  df_clean$log_duration <- log(df_clean$duration) 
  print("log_duration") 
} 

if("vorticity" %in% names(df_clean)) { 
  df_clean$vorticity_scaled <- as.numeric(scale(df_clean$vorticity)) 
  print("vorticity_scaled") 
} 

df_clean$site <- as.factor(df_clean$site) 
df_clean$site <- relevel(df_clean$site, ref = "MT") 


gam_reml <- gam(
  log_duration ~ site + 
    s(depth, bs = "cr", k = 10) + 
    s(speed, by = site, bs = "cr", k = 10) + 
    s(vorticity_scaled, bs = "cr", k = 10) + 
    ti(depth, speed, bs = "cr", k = 10),
  data = df_clean,      
  family = gaussian(),   
  method = "REML"
)

summary(gam_reml)

draw(gam_reml, residuals = TRUE)
