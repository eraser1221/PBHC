library(tidyverse)

data <- read.csv("contine_time.csv")

df_long <- data %>%
  pivot_longer(cols = everything(), names_to = "Site", values_to = "Seconds") %>%
  mutate(Minutes = Seconds / 60) %>%
  mutate(Time_Bin = cut(
    Minutes, 
    breaks = c(0, 4, 8, 12, Inf), 
    labels = c("0-3 min", "4-7 min", "8-11 min", ">12 min"), 
    right = FALSE 
  )) %>%
  filter(!is.na(Time_Bin))

contingency_table <- table(df_long$Time_Bin, df_long$Site)
chisq_res <- chisq.test(contingency_table)
print("Observed Counts:")
print(contingency_table)

residuals_df <- as.data.frame(as.table(chisq_res$stdres))
colnames(residuals_df) <- c("Time_Bin", "Site", "Residual")

residuals_df$Site <- factor(residuals_df$Site, levels = c("LS", "NMZ", "YZ", "MT"))


library(corrplot)

res_matrix <- chisq_res$stdres
res_matrix_ordered <- res_matrix[, c("LS", "NMZ", "YZ", "MT")]

corrplot(res_matrix_ordered, 
         is.corr = FALSE, 
         method = "circle",
         tl.col = "black", 
         tl.srt = 0,    
         cl.pos = "r",  
         title = "Standardized Residuals (Ordered)",
         mar = c(0,0,2,0))
