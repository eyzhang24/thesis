library(tidyverse)
library(latex2exp)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

# function to transform variable names
get_var <- function(orig) {
  result <- substr(orig, 4, 6)
  return(result)
}

########
# observed data
########

# correlation between exposures
target_first <- read_rds("nhanes_data/processed_data.RDS")

target_first2 <- target_first |> 
  rename_at(vars(2:19), ~substr(., 4, 6)) |> 
  select(2:19) 
target_first2 <- target_first2[, names(target_first2)[order(names(target_first2))]]

cor_mat <- cor(target_first2, method = "spearman")
cor_mat[lower.tri(cor_mat)] <- NA
melt_cor <- reshape2::melt(cor_mat) |> 
  mutate(label = ifelse(value == 1, NA, round(value, 2))) |> 
  na.omit()
cor_orig <- melt_cor |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = c(0, .25, .5, 0.75, 1),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = TeX(r"( Spearman's $\rho$ )")) +
  theme(
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

ggsave("nhanes_figs/ch4_corr.png", width = 6, height = 6)

# log-transformed versus original
comb_small <- read_rds("nhanes_data/processed_data.RDS")

#log-transform target data
comb_log <- comb_small |> 
  mutate(across(c(LBX074LA:LBXF08LA, LBXCOT, TELOMEAN), log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric)) |> 
  #re-order variables
  select(SEQN:LBXF08LA, #ID and pop's
         BMXBMI:LBXBAPCT, #bio
         RIDAGEYR:DMDEDUC2, #demo, gender is binary, re and educ multi
         TELOMEAN) |> 
  #gender to 0's and 1's: 0 male, 1 female
  #should really be called sex
  mutate(RIAGENDR = RIAGENDR - 1)

univ1 <- comb_small |> 
  select(2:19) |> 
  pivot_longer(cols = 1:18) |> 
  mutate(name = substr(name, 4, 6)) |> 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free", nrow = 3) +
  labs(x = "Concentration (ng/mL)")
univ2 <- comb_log |> 
  select(2:19) |> 
  pivot_longer(cols = 1:18) |> 
  mutate(name = substr(name, 4, 6)) |> 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free", nrow = 3) +
  labs(x = "Natural log concentration (log(ng/mL))")

cowplot::plot_grid(univ1, univ2, labels = "auto", nrow = 2)

ggsave("nhanes_figs/ch4_univlog.png", width = 8.5, height = 8.5)

# covariates
cov_cont <- comb_log |> 
  select(BMXBMI:RIDAGEYR) |> 
  pivot_longer(cols = 1:9) |> 
  # rename covariates
  mutate(name = case_when(
    grepl("BMXBMI", name) ~ "BMI", 
    grepl("LBXCOT", name) ~ "Cotinine levels", 
    grepl("WBCSI", name) ~ "White blood cell count", 
    grepl("LYPCT", name) ~ "Lymphocyte %", 
    grepl("MOPCT", name) ~ "Monocyte %", 
    grepl("NEPCT", name) ~ "Neutrophil %", 
    grepl("EOPCT", name) ~ "Eosinophil %", 
    grepl("BAPCT", name) ~ "Basophil %", 
    grepl("RIDAGEYR", name) ~ "Age", 
    .default = NA
  )) |> 
  mutate(name = factor(name, levels = c(
    "Age", "BMI", "Cotinine levels", "White blood cell count", "Lymphocyte %", 
    "Monocyte %", "Neutrophil %", "Eosinophil %", "Basophil %"
  ))) |> 
  ggplot(aes(x = value)) +
  geom_density() + 
  facet_wrap(~name, scales = "free", ncol = 2) 
cov_cont

df_forcovcat <- comb_log |> 
  select(RIAGENDR, RIDRETH1, DMDEDUC2) |> 
  # re-code covariates
  mutate(RIAGENDR = ifelse(RIAGENDR == 0, "Male", "Female"), 
         RIDRETH1 = case_when(
           RIDRETH1 == 1 ~ "Mexican American", 
           RIDRETH1 == 2 ~ "Other Hispanic", 
           RIDRETH1 == 3 ~ "Non-Hisp. White", 
           RIDRETH1 == 4 ~ "Non-Hisp. Black", 
           RIDRETH1 == 5 ~ "Non-Hisp. Other", 
           .default = NA
         ), 
         DMDEDUC2 = case_when(
           DMDEDUC2 == 1 | DMDEDUC2 == 2 ~ "Less than HS", 
           DMDEDUC2 == 3 ~ "HS", 
           DMDEDUC2 == 4 ~ "Some College", 
           DMDEDUC2 == 5 ~ "College", 
           .default = NA
         )) |> 
  pivot_longer(cols = 1:3) |> 
  mutate(value = factor(
    value, levels = rev(c("Male", "Female", 
                          "Mexican American", "Other Hispanic", 
                          "Non-Hisp. White", "Non-Hisp. Black", "Non-Hisp. Other", 
                          "Less than HS", "HS", "Some College", "College"))
  ), 
  name = factor(name, levels = c("RIAGENDR", "RIDRETH1", "DMDEDUC2"), 
                labels = c("Sex", "Race/Ethnicity", "Education"))) 

cov_cat <- df_forcovcat |> 
  ggplot(aes(x = value)) +
  geom_bar(stat = "count", fill = "gray") +
  geom_text(aes(label = after_stat(count)), stat = "count", 
             size = 3, hjust = 1, nudge_y = -2) +
  facet_wrap(~name, scales = "free", ncol = 1) +
  coord_flip() +
  labs(x = NULL)
cov_cat

cowplot::plot_grid(cov_cont, cov_cat, labels = "auto", nrow = 1, 
                   rel_widths = c(0.5, 0.5))
ggsave("nhanes_figs/ch4_covdist.png", width = 7.5, height = 6)

#look at correlation between race and chemicals
name_order <- substr(names(comb_log)[order(substr(names(comb_log)[2:19], 4, 6))+1], 4, 6)
comb_log |> 
  select(c(30, 2:19)) |> 
  pivot_longer(cols = 2:19, names_to = "key", values_to = "value") |> 
  mutate(key = factor(substr(key, 4, 6), levels = name_order)) |> 
  mutate(RIDRETH1 = as.factor(RIDRETH1)) |> 
  ggplot(aes(x = RIDRETH1, y = value, color = RIDRETH1)) +
  geom_boxplot() +
  scale_color_discrete(
    name = NULL, 
    labels = c("Mexican American", "Hispanic Other", "Non-Hisp. White", 
               "Non-Hisp. Black", "Non-Hisp. Other")) + 
  theme(legend.spacing.y = unit(0.25, 'cm'), 
        legend.position = "top", 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  guides(color = guide_legend(byrow = TRUE)) +
  labs(x = NULL, y = "Log concentration") +
  facet_wrap(~key, scales = "free_y", nrow = 2) 
ggsave("nhanes_figs/ch4_race_exp.png", width = 7.5, height = 4)


########
# simulated data small
########

#read small size back in
out1 <- read_rds("nhanes_sim/sim_preds_sm.RDS")

comb_sim1 <- bind_rows(out1)

#density plots for exposures
name_order <- substr(names(comb_log)[order(substr(names(comb_log)[2:19], 4, 6))+1], 4, 6)
                # "age", "bmi")
comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(1:18, sim) |> 
  pivot_longer(cols = 1:18) |>
  mutate(name = factor(substr(name, 4, 6), levels = name_order)) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(2:19) |> pivot_longer(cols = 1:18) |> 
      mutate(name = factor(substr(name, 4, 6), levels = name_order)),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free")
ggsave("nhanes_figs/ch4_univ_exp_sim.png", width = 7.5, height = 6)

# function to transform covariates 
get_cov_name <- function(name) {
  result <- case_when(
    grepl("BMXBMI", name) ~ "BMI", 
    grepl("LBXCOT", name) ~ "Cotinine levels", 
    grepl("WBCSI", name) ~ "White blood cell count", 
    grepl("LYPCT", name) ~ "Lymphocyte %", 
    grepl("MOPCT", name) ~ "Monocyte %", 
    grepl("NEPCT", name) ~ "Neutrophil %", 
    grepl("EOPCT", name) ~ "Eosinophil %", 
    grepl("BAPCT", name) ~ "Basophil %", 
    grepl("RIDAGEYR", name) ~ "Age", 
    .default = NA
  )
  return(result)
}

#density plot for covariates
cov_sim_p <- comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(BMXBMI:RIDAGEYR, sim) |> 
  pivot_longer(cols = 1:9) |>
  mutate(name = get_cov_name(name)) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(BMXBMI:RIDAGEYR) |> pivot_longer(cols = 1:9) |> 
      mutate(name = get_cov_name(name)),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free", ncol = 2)
#bar plot + violin plot for covariates
cov_sim_q <- comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(sim, RIAGENDR, RIDRETH1, DMDEDUC2) |> 
  # re-code covariates
  mutate(RIAGENDR = ifelse(RIAGENDR == 0, "Male", "Female"), 
         RIDRETH1 = case_when(
           RIDRETH1 == 1 ~ "Mexican American", 
           RIDRETH1 == 2 ~ "Other Hispanic", 
           RIDRETH1 == 3 ~ "Non-Hisp. White", 
           RIDRETH1 == 4 ~ "Non-Hisp. Black", 
           RIDRETH1 == 5 ~ "Non-Hisp. Other", 
           .default = NA
         ), 
         DMDEDUC2 = case_when(
           DMDEDUC2 == 1 | DMDEDUC2 == 2 ~ "Less than HS", 
           DMDEDUC2 == 3 ~ "HS", 
           DMDEDUC2 == 4 ~ "Some College", 
           DMDEDUC2 == 5 ~ "College", 
           .default = NA
         )) |> 
  pivot_longer(cols = 2:4) |> 
  group_by(sim, name, value) |> 
  summarize(prop = n()/250) |> 
  mutate(value = factor(
    value, levels = rev(c("Male", "Female", 
                          "Mexican American", "Other Hispanic", 
                          "Non-Hisp. White", "Non-Hisp. Black", "Non-Hisp. Other", 
                          "Less than HS", "HS", "Some College", "College"))
  ), 
  name = factor(name, levels = c("RIAGENDR", "RIDRETH1", "DMDEDUC2"), 
                labels = c("Sex", "Race/Ethnicity", "Education"))) |> 
  ggplot(aes(x = value, y = prop)) +
  geom_bar(data = df_forcovcat, aes(x = value, y = after_stat(prop), group = 1), 
           inherit.aes = FALSE,
           stat = "count", fill = "skyblue") +
  geom_violin(color = "gray30", fill = "gray", alpha = 0.25) + 
  facet_wrap(~name, scales = "free", ncol = 1) +
  coord_flip() +
  labs(x = NULL)
# cov_sim_q

cowplot::plot_grid(cov_sim_p, cov_sim_q, labels = "auto", nrow = 1, 
                   rel_widths = c(0.5, 0.5))
  
ggsave("nhanes_figs/ch4_univ_cov_sim.png", width = 7.5, height = 6)

#scatterplots???

#look at average correlation values 

cors <- out1 |> 
  purrr::map_df(\(x) {
    x <- x |> 
      rename_at(vars(1:18), ~substr(., 4, 6)) |> 
      select(1:18) 
    x <- x[, names(x)[order(names(x))]]
    cor_mat <- cor(x[, 1:18], method = "spearman")
    cor_mat[lower.tri(cor_mat)] <- NA
    melt_cor <- reshape2::melt(cor_mat) |> 
      mutate(value = ifelse(value == 1, NA, value)) |> 
      na.omit() |> 
      mutate(sim = x$sim[1])
    return(melt_cor)
  })

newbreaks <- function(lims) {
  range <- diff(lims)
  return(c(lims[1] + range/5, mean(lims), lims[2] - range/5))
}

#distribution
cors_dens <- cors |> 
  group_by(Var1, Var2) |> 
  mutate(mean = mean(value)) |> 
  ungroup() |> 
  ggplot(aes(x = value, fill = mean)) +
  geom_density() + 
  scale_x_continuous(breaks = newbreaks, 
                     labels = ~round(.x, 2)) +
  scale_y_continuous(position = "right") +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = seq(0, 1, by = 0.25),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  ggh4x::facet_grid2(fct_rev(Var2) ~ Var1, scales = "free", independent = "x", 
                     render_empty = FALSE, switch = "both") +
  labs(x = TeX(r"( Spearman's $\rho$ )")) +
  theme(strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0), 
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cors_dens

ggsave("nhanes_figs/ch4_corr_sim.png", width = 10, height = 10)

#average value
cor_sim <- cors |> 
  group_by(Var1, Var2) |> 
  summarize(value = mean(value)) |> 
  mutate(label = round(value, 2)) |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = seq(0, 1, by = 0.25),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = TeX(r"( Mean Spearman's $\rho$ )")) +
  theme(
    panel.grid.major.x = element_line(color = "grey85",
                                      linewidth = 0.25,
                                      linetype = 2), 
    panel.border = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cor_sim

ggsave("nhanes_figs/ch4_corr_avg_sim.png", width = 6, height = 6)

#put original and simulated correlation together
top_row <- cowplot::plot_grid(cor_orig, cor_sim, labels = "auto", label_size = 16, 
                   nrow = 1, scale = 0.95)
top_row
ggsave("nhanes_figs/ch4_corr_sim+orig.png", width = 10, height = 5)

#put original, simulated, and density correlation together
cowplot::plot_grid(top_row, cors_dens, labels = c("", "c"), label_size = 16, 
                   nrow = 2, rel_heights = c(5, 7), scale = c(1, 0.95))
ggsave("nhanes_figs/ch4_corr_simorigdens.png", width = 10, height = 12)

#############
# simulated data large
#############
#read small size back in
out2 <- read_rds("nhanes_sim/sim_preds_lg.RDS")

comb_sim2 <- bind_rows(out2)

one <- comb_sim2[1:2000, ]

#density plots for exposures
name_order <- substr(names(comb_log)[order(substr(names(comb_log)[2:19], 4, 6))+1], 4, 6)
# "age", "bmi")
comb_sim2 |> 
  mutate(sim = as.factor(sim)) |> 
  select(1:18, sim) |> 
  pivot_longer(cols = 1:18) |>
  mutate(name = factor(substr(name, 4, 6), levels = name_order)) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(2:19) |> pivot_longer(cols = 1:18) |> 
      mutate(name = factor(substr(name, 4, 6), levels = name_order)),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free")
ggsave("nhanes_figs/ch4_univ_exp_sim_lg.png", width = 7.5, height = 6)

#density plot for covariates
cov_sim_p2 <- comb_sim2 |> 
  mutate(sim = as.factor(sim)) |> 
  select(BMXBMI:RIDAGEYR, sim) |> 
  pivot_longer(cols = 1:9) |>
  mutate(name = get_cov_name(name)) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(BMXBMI:RIDAGEYR) |> pivot_longer(cols = 1:9) |> 
      mutate(name = get_cov_name(name)),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free", ncol = 2)
#bar plot + violin plot for covariates
cov_sim_q2 <- comb_sim2 |> 
  mutate(sim = as.factor(sim)) |> 
  select(sim, RIAGENDR, RIDRETH1, DMDEDUC2) |> 
  # re-code covariates
  mutate(RIAGENDR = ifelse(RIAGENDR == 0, "Male", "Female"), 
         RIDRETH1 = case_when(
           RIDRETH1 == 1 ~ "Mexican American", 
           RIDRETH1 == 2 ~ "Other Hispanic", 
           RIDRETH1 == 3 ~ "Non-Hisp. White", 
           RIDRETH1 == 4 ~ "Non-Hisp. Black", 
           RIDRETH1 == 5 ~ "Non-Hisp. Other", 
           .default = NA
         ), 
         DMDEDUC2 = case_when(
           DMDEDUC2 == 1 | DMDEDUC2 == 2 ~ "Less than HS", 
           DMDEDUC2 == 3 ~ "HS", 
           DMDEDUC2 == 4 ~ "Some College", 
           DMDEDUC2 == 5 ~ "College", 
           .default = NA
         )) |> 
  pivot_longer(cols = 2:4) |> 
  group_by(sim, name, value) |> 
  summarize(prop = n()/1003) |> 
  mutate(value = factor(
    value, levels = rev(c("Male", "Female", 
                          "Mexican American", "Other Hispanic", 
                          "Non-Hisp. White", "Non-Hisp. Black", "Non-Hisp. Other", 
                          "Less than HS", "HS", "Some College", "College"))
  ), 
  name = factor(name, levels = c("RIAGENDR", "RIDRETH1", "DMDEDUC2"), 
                labels = c("Sex", "Race/Ethnicity", "Education"))) |> 
  ggplot(aes(x = value, y = prop)) +
  geom_bar(data = df_forcovcat, aes(x = value, y = after_stat(prop), group = 1), 
           inherit.aes = FALSE, stat = "count", fill = "skyblue") +
  geom_violin(color = "gray30", fill = "gray", alpha = 0.25) + 
  facet_wrap(~name, scales = "free", ncol = 1) +
  coord_flip() +
  labs(x = NULL, y = "proportion")

cov_sim_q2

cowplot::plot_grid(cov_sim_p2, cov_sim_q2, labels = "auto", nrow = 1, 
                   rel_widths = c(0.5, 0.5))

ggsave("nhanes_figs/ch4_univ_cov_sim_lg.png", width = 7.5, height = 6)

#correlation plots, large
corl <- out2 |> 
  purrr::map_df(\(x) {
    x <- x |> 
      rename_at(vars(1:18), ~substr(., 4, 6)) |> 
      select(1:18) 
    x <- x[, names(x)[order(names(x))]]
    cor_mat <- cor(x[, 1:18], method = "spearman")
    cor_mat[lower.tri(cor_mat)] <- NA
    melt_cor <- reshape2::melt(cor_mat) |> 
      mutate(value = ifelse(value == 1, NA, value)) |> 
      na.omit() |> 
      mutate(sim = x$sim[1])
    return(melt_cor)
  })

newbreaks <- function(lims) {
  range <- diff(lims)
  return(c(lims[1] + range/5, mean(lims), lims[2] - range/5))
}
#distribution
cors_dens2 <- corl |> 
  group_by(Var1, Var2) |> 
  mutate(mean = mean(value)) |> 
  ungroup() |> 
  ggplot(aes(x = value, fill = mean)) +
  geom_density() + 
  scale_x_continuous(breaks = newbreaks, 
                     labels = ~round(.x, 2)) +
  scale_y_continuous(position = "right") +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = seq(0, 1, by = 0.25),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  ggh4x::facet_grid2(fct_rev(Var2) ~ Var1, scales = "free", independent = "x", 
                     render_empty = FALSE, switch = "both") +
  labs(x = TeX(r"( Spearman's $\rho$ )")) +
  theme(strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0), 
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cors_dens2

# ggsave("index/figures/ch4_corr_sim.png", width = 10, height = 7)

#average value
cor_sim2 <- corl |> 
  group_by(Var1, Var2) |> 
  summarize(value = mean(value)) |> 
  mutate(label = round(value, 2)) |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    limit = c(0, 1), breaks = seq(0, 1, by = 0.25),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = TeX(r"( Mean Spearman's $\rho$ )")) +
  theme(
    panel.grid.major.x = element_line(color = "grey85",
                                      linewidth = 0.25,
                                      linetype = 2), 
    panel.border = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
cor_sim2

# ggsave("index/figures/ch4_corr_avg_sim.png", width = 5, height = 5)

#put original and simulated correlation together
top_row2 <- cowplot::plot_grid(cor_orig, cor_sim2, labels = "auto", label_size = 16, 
                              nrow = 1, scale = 0.95)
top_row2
ggsave("nhanes_figs/ch4_corr_sim+orig.png", width = 12, height = 6)

# #put original, simulated, and density correlation together
# cowplot::plot_grid(top_row2, cors_dens2, labels = c("", "c"), label_size = 16, 
#                    nrow = 2, rel_heights = c(5, 7), scale = c(1, 0.95))
# ggsave("index/figures/ch4_corr_lg_simorigdens.png", width = 10, height = 12)


#############
#match label with equation
############
#for mlr models
equations1 <-  c(TeX("No inter", output = "character"), 
                 TeX("0.475$X_{D05}*X_{194}$"), TeX("0.24$X_{D05}*(X_{194}-1)^2$"), 
                TeX("0.5$X_{F08}*X_{F03}$"), TeX("0.17$X_{F08}*(X_{F03}-1)^2$"), 
                TeX("0.36$X_{074}*X_{194}$"), TeX("0.25$X_{074}*(X_{094}-1)^2$"), 
                TeX("0.3$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.125$X_{D05}*(X_{PCB}-1)^2*X_{194}$"), 
                TeX("0.9$X_{D05}*X_{194}$"), TeX("0.48$X_{D05}*(X_{194}-1)^2$"), 
                TeX("$X_{F08}*X_{F03}$"), TeX("0.34$X_{F08}*(X_{F03}-1)^2$"), 
                TeX("0.6$X_{074}*X_{194}$"), TeX("0.5$X_{074}*(X_{194}-1)^2$"), 
                TeX("0.72$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.25$X_{D05}*(X_{PCB}-1)^2*X_{194}$"))
names1 <- c("_base",
  "am1", "ap1", "bm1", "bp1", "cm1", "cp1", "dm1", "dp1",
  "am2", "ap2", "bm2", "bp2", "cm2", "cp2", "dm2", "dp2")
appender1 <- function(string) {
  return(equations1[match(string, names1)])
}

#for oracle models
equations <-  c(TeX("0.475$X_{D05}*X_{194}$"), TeX("0.24$X_{D05}*(X_{194}-1)^2$"), 
                TeX("0.5$X_{F08}*X_{F03}$"), TeX("0.17$X_{F08}*(X_{F03}-1)^2$"), 
                TeX("0.36$X_{074}*X_{194}$"), TeX("0.25$X_{074}*(X_{094}-1)^2$"), 
                TeX("0.3$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.125$X_{D05}*(X_{PCB}-1)^2*X_{194}$"), 
                TeX("0.9$X_{D05}*X_{194}$"), TeX("0.48$X_{D05}*(X_{194}-1)^2$"), 
                TeX("$X_{F08}*X_{F03}$"), TeX("0.34$X_{F08}*(X_{F03}-1)^2$"), 
                TeX("0.6$X_{074}*X_{194}$"), TeX("0.5$X_{074}*(X_{194}-1)^2$"), 
                TeX("0.72$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.25$X_{D05}*(X_{PCB}-1)^2*X_{194}$"))
names <- c(
  "am1", "ap1", "bm1", "bp1", "cm1", "cp1", "dm1", "dp1",
  "am2", "ap2", "bm2", "bp2", "cm2", "cp2", "dm2", "dp2")
appender <- function(string) {
  return(equations[match(string, names)])
}



###############
#mlr models
##############

#mlr model output small
mlrmods <- read_rds("nhanes_sim/_mlr/mlr_mods_sm.RDS")
rsquared1 <- mlrmods |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(mlrmods))

rsquared1 |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/mlr_sm_rsq.png", width = 7.5, height = 5)

#mlr model output large
mlrmodl <- read_rds("sim/mlr/mlr_mods_lg.RDS")
rsquared1l <- mlrmodl |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(mlrmodl))

rsquared1l |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/mlr_lg_rsq.png", width = 7.5, height = 5)

###############
#chem models
##############

#chem model output small
chem_mods <- read_rds("sim/mlr/chem_oracle_sm.RDS")
rsq_chem <- chem_mods |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_mods))

rsqsmplot <- rsq_chem |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = latex2exp::TeX("R$^2$"))
rsqsmplot
ggsave("sim/figs/chem_sm_rsq.png", width = 7, height = 5)

#chem model output large
chem_modl <- read_rds("sim/mlr/chem_oracle_lg.RDS")
rsq_cheml <- chem_modl |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_modl))

rsqlgplot <- rsq_cheml |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = latex2exp::TeX("R$^2$"))
cowplot::plot_grid(rsqsmplot, rsqlgplot, labels = "auto", nrow = 1)
ggsave("index/figures/chem_rsq.png", width = 12, height = 6)

###############
#oracle models
##############

# #small models
# oracle_mods <- read_rds("sim/oracle/oracle_mods_sm.RDS")
# 
# #extract r squared
# rsquared2 <- oracle_mods |> 
#   purrr::map_df(\(x) {
#     data.frame(
#       rsq = summary(x)$r.squared
#     )
#   }) |> 
#   mutate(name = names(oracle_mods))
# 
# rsquared2 |> 
#   ggplot(aes(x = rsq)) +
#   geom_density() + 
#   facet_wrap(~name, scales = "free_y", 
#              labeller = as_labeller(appender1, 
#                                     default = label_parsed), 
#              ncol = 4) +
#   labs(y = "Density", x = TeX("R$^2$"))
# ggsave("sim/figs/oracle_sm_rsq.png", width = 7.5, height = 5)
# 
# #extract p values
# keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
#                'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
#                'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')
# pval <- oracle_mods |> 
#   purrr::map_df(\(x) {
#     x <- summary(x)$coefficients[,4]
#     return(data.frame(pval = x[!(names(x) %in% keepnames)]))
#   }) |> 
#   mutate(name = names(oracle_mods)[101:1700])
# 
# pval |> 
#   ggplot(aes(x = pval)) +
#   geom_density() + 
#   facet_wrap(~name, scales = "free_y", 
#              labeller = as_labeller(appender, 
#                                     default = label_parsed)) +
#   labs(y = "Density", x = "p-value")
# ggsave("sim/figs/oracle_sm_pval_dist.png", width = 7.5, height = 5)
# 
# #extract power at alpha = 0.05
# oracle_sm_power <- pval |> 
#   group_by(name) |> 
#   summarize(power = sum(pval < 0.05)/n())
# 
# 
# write_csv(oracle_sm_power, "sim/tables/oracle_sm_power.csv")
# 
# #large models
# oracle_modl <- read_rds("sim/oracle/oracle_mods_lg.RDS")
# 
# #extract r squared
# rsquared2l <- oracle_modl |> 
#   purrr::map_df(\(x) {
#     data.frame(
#       rsq = summary(x)$r.squared
#     )
#   }) |> 
#   mutate(name = names(oracle_modl))
# 
# rsquared2l |> 
#   ggplot(aes(x = rsq)) +
#   geom_density() + 
#   facet_wrap(~name, scales = "free_y", 
#              labeller = as_labeller(appender1, 
#                                     default = label_parsed), 
#              ncol = 4) +
#   labs(y = "Density", x = TeX("R$^2$"))
# ggsave("sim/figs/oracle_sm_rsq.png", width = 7.5, height = 5)
# 
# #p-value
# pval_large <- oracle_modl |> 
#   purrr::map_df(\(x) {
#     x <- summary(x)$coefficients[,4]
#     return(data.frame(pval = x[!(names(x) %in% keepnames)]))
#   }) |> 
#   mutate(name = names(oracle_modl)[101:1700]) 
# 
# pval_large |> 
#   ggplot(aes(x = pval)) +
#   geom_density() + 
#   facet_wrap(~name, scales = "free_y", 
#              labeller = as_labeller(appender, 
#                                     default = label_parsed)) +
#   labs(y = "Density", x = "p-value")
# ggsave("sim/figs/oracle_lg_pval_dist.png", width = 7.5, height = 5)
# 
# #extract power at alpha = 0.05
# oracle_lg_power <- pval_large |> 
#   group_by(name) |> 
#   summarize(power = sum(pval < 0.05)/n())
# 
# write_csv(oracle_lg_power, "sim/tables/oracle_lg_power.csv")
