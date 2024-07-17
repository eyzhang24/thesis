library(tidyverse)
library(latex2exp)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

# read in data
comb <- read_rds("nhanes_data/processed_data.RDS")

cor_mat <- cor(comb[, 2:19], method = "spearman")
cor_mat[lower.tri(cor_mat)] <- NA
melt_cor <- reshape2::melt(cor_mat) |> 
  na.omit() |> 
  mutate(label = ifelse(value == 1, NA, round(value, 2))) #|> 
  # na.omit()
cor_orig <- melt_cor |> 
  ggplot(aes(x = Var1, y = Var2, fill = label)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    midpoint = 0.5, na.value = "grey") +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = TeX(r"( Spearman's $\rho$ )")) +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1), 
    panel.grid.major.x = element_line(color = "grey85",
                                      linewidth = 0.25,
                                      linetype = 2), 
    panel.border = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cor_orig

cor_mat <- cor(select(comb, LBXF08LA, LBX074LA, LBX194LA, LBXPCBLA, 
                      LBXF04LA, LBXD05LA, LBXF03LA), method = "spearman")
cor_mat[lower.tri(cor_mat)] <- NA
melt_cor <- reshape2::melt(cor_mat) |> 
  na.omit() |> 
  mutate(label = ifelse(value == 1, NA, round(value, 2))) #|> 
# na.omit()
cor_orig <- melt_cor |> 
  ggplot(aes(x = Var1, y = Var2, fill = label)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    midpoint = 0.5, na.value = "grey") +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = TeX(r"( Spearman's $\rho$ )")) +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1), 
    panel.grid.major.x = element_line(color = "grey85",
                                      linewidth = 0.25,
                                      linetype = 2), 
    panel.border = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cor_orig

