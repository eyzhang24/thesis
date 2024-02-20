library(tidyverse)
library(latex2exp)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

########
# before simulation
########

# correlation between exposures
target_first <- read_csv("madres_data/target_first.csv")

cor_mat <- cor(target_first[, 5:14], method = "spearman")
cor_mat[lower.tri(cor_mat)] <- NA
melt_cor <- reshape2::melt(cor_mat) |> 
  mutate(label = ifelse(value == 1, NA, round(value, 2))) |> 
  na.omit()
cor_orig <- melt_cor |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    limit = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
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

ggsave("index/figures/ch4_corr.png", width = 5, height = 5)

# log-transformed versus original
comb_small <- read_csv("madres_data/base_data.csv")

#log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric))

univ1 <- comb_small |> 
  select(10:19) |> 
  pivot_longer(cols = 1:10) |> 
  mutate(name = factor(name, levels = names(comb_small)[10:19])) |> 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free", nrow = 2) +
  labs(x = "Concentration (ng/mL)")
univ2 <- comb_log |> 
  select(10:19) |> 
  pivot_longer(cols = 1:10) |> 
  mutate(name = factor(name, levels = names(comb_small)[10:19])) |> 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free", nrow = 2) +
  labs(x = "Natural log concentration (log(ng/mL))")

cowplot::plot_grid(univ1, univ2, labels = "auto", nrow = 2)

ggsave("index/figures/ch4_univlog.png", width = 7.5, height = 5)

# covariates
cov_cont <- comb_log |> 
  select(age, bmi) |> 
  pivot_longer(cols = 1:2) |> 
  ggplot(aes(x = value)) +
  geom_density() + 
  facet_wrap(~name, scales = "free", ncol = 1) 
# cov_cont

df_forcovcat <- comb_log |> 
  select(smoke, race) |> 
  mutate(smoke = ifelse(smoke == 0, "Never-exposed", "Ever-exposed"), 
         race = case_when(
           race == 1 ~ "Non-Hisp. white", 
           race == 2 ~ "Non-Hisp. black", 
           race == 3 ~ "Non-Hisp. other", 
           race == 4 ~ "Hispanic born\nin US", 
           race == 5 ~ "Hispanic born\noutside US"
         )) |> 
  pivot_longer(cols = 1:2) |> 
  mutate(value = factor(
    value, levels = rev(c("Never-exposed", "Ever-exposed",
                      "Non-Hisp. white", "Non-Hisp. black", "Non-Hisp. other", 
                      "Hispanic born\nin US", "Hispanic born\noutside US"))
  )) 

cov_cat <- df_forcovcat |> 
  ggplot(aes(x = value)) +
  geom_bar(stat = "count", fill = "gray") +
  geom_text(aes(label = after_stat(count)), stat = "count", 
             size = 3, hjust = 1, nudge_y = -2) +
  facet_wrap(~name, scales = "free", ncol = 1) +
  coord_flip() +
  labs(x = NULL)
# cov_cat

cowplot::plot_grid(cov_cont, cov_cat, labels = "auto", nrow = 1, 
                   rel_widths = c(0.4, 0.6))
ggsave("index/figures/ch4_covdist.png", width = 6, height = 4)

#look at correlation between race and chemicals
name_order <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")
comb_log |> 
  select(c(6, 10:19)) |> 
  pivot_longer(cols = 2:11, names_to = "key", values_to = "value") |> 
  mutate(key = factor(key, levels = name_order)) |> 
  mutate(race = as.factor(race)) |> 
  ggplot(aes(x = race, y = value, color = race)) +
  geom_boxplot() +
  scale_color_discrete(
    name = "Race by ethnicity\nand birth place", 
    labels = c("Non-Hisp. white", "Non-Hisp. black", "Non-Hisp. other", 
               "Hispanic born\nin US", "Hispanic born\noutside US")) + 
  theme(legend.spacing.y = unit(0.25, 'cm')) +
  guides(color = guide_legend(byrow = TRUE)) +
  labs(x = "Race, coded", y = "Log concentration") +
  facet_wrap(~key, scales = "free_y")
ggsave("index/figures/ch4_race_exp.png", width = 6, height = 4)

########
# after simulation
########

#read small size back in
out1 <- read_rds("sim/sim_preds_sm.RDS")

comb_sim1 <- bind_rows(out1)

#density plots for exposures
name_order <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")
                # "age", "bmi")
comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(5:15) |> 
  pivot_longer(cols = 1:10) |>
  mutate(name = factor(name, levels = name_order)) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(10:19) |> pivot_longer(cols = 1:10) |> 
      mutate(name = factor(name, levels = name_order)),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free")
ggsave("index/figures/ch4_univ_exp_sim.png", width = 6, height = 4)

#density plot for covariates
cov_sim_p <- comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(age, bmi, sim) |> 
  pivot_longer(cols = 1:2) |>
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(age, bmi) |> pivot_longer(cols = 1:2),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free", ncol = 1)
#bar plot + violin plot for covariates
cov_sim_q <- comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(sim, smoke, race) |> 
  mutate(smoke = ifelse(smoke == 0, "Never-exposed", "Ever-exposed"), 
         race = case_when(
           race == 1 ~ "Non-Hisp. white", 
           race == 2 ~ "Non-Hisp. black", 
           race == 3 ~ "Non-Hisp. other", 
           race == 4 ~ "Hispanic born\nin US", 
           race == 5 ~ "Hispanic born\noutside US"
         )) |> 
  pivot_longer(cols = 2:3) |> 
  group_by(sim, name, value) |> 
  summarize(count = n()) |> 
  mutate(value = factor(
    value, levels = rev(c("Never-exposed", "Ever-exposed",
                          "Non-Hisp. white", "Non-Hisp. black", "Non-Hisp. other", 
                          "Hispanic born\nin US", "Hispanic born\noutside US"))
  )) |> 
  ggplot(aes(x = value, y = count)) +
  geom_bar(data = df_forcovcat, aes(x = value), inherit.aes = FALSE,
           stat = "count", fill = "skyblue") +
  geom_violin(color = "gray30", fill = "gray", alpha = 0.25) + 
  facet_wrap(~name, scales = "free", ncol = 1) +
  coord_flip() +
  labs(x = NULL)
# cov_sim_q

cowplot::plot_grid(cov_sim_p, cov_sim_q, labels = "auto", nrow = 1, 
                   rel_widths = c(0.4, 0.6))
  
ggsave("index/figures/ch4_univ_cov_sim.png", width = 6, height = 4)

#scatterplots???

#look at average correlation values 

cors <- out1 |> 
  purrr::map_df(\(x) {
    cor_mat <- cor(x[, 5:14], method = "spearman")
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
    limit = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
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

ggsave("index/figures/ch4_corr_sim.png", width = 10, height = 7)

#average value
cor_sim <- cors |> 
  group_by(Var1, Var2) |> 
  summarize(value = mean(value)) |> 
  mutate(label = round(value, 2)) |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    limit = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
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
cor_sim

ggsave("index/figures/ch4_corr_avg_sim.png", width = 5, height = 5)

#put original and simulated correlation together
top_row <- cowplot::plot_grid(cor_orig, cor_sim, labels = "auto", label_size = 16, 
                   nrow = 1, scale = 0.95)
top_row
ggsave("index/figures/ch4_corr_sim+orig.png", width = 10, height = 5)

#put original, simulated, and density correlation together
cowplot::plot_grid(top_row, cors_dens, labels = c("", "c"), label_size = 16, 
                   nrow = 2, rel_heights = c(5, 7), scale = c(1, 0.95))
ggsave("index/figures/ch4_corr_simorigdens.png", width = 10, height = 12)

#############
#match label with equation
############
#for mlr models
equations1 <-  c(TeX("No inter", output = "character"), 
                 TeX("0.3Hg$*$Ni"), TeX("0.1Hg$*($Ni$-1)^2$"), 
                TeX("0.3Cd$*$As"), TeX("0.1Cd$*($As$-1)^2$"), 
                TeX("0.3Hg$*$Co"), TeX("0.1Hg$*($Co$-1)^2$"), 
                TeX("0.3Hg$*$Ni$*$Tl"), TeX("0.1Hg$*($Ni$-1)^2*$Tl"), 
                TeX("0.6Hg$*$Ni"), TeX("0.2Hg$*($Ni$-1)^2$"), 
                TeX("0.6Cd$*$As"), TeX("0.2Cd$*($As$-1)^2$"), 
                TeX("0.6Hg$*$Co"), TeX("0.2Hg$*($Co$-1)^2$"), 
                TeX("0.6Hg$*$Ni$*$Tl"), TeX("0.2Hg$*($Ni$-1)^2*$Tl"))
names1 <- c("_base",
  "am1", "ap1", "bm1", "bp1", "cm1", "cp1", "dm1", "dp1",
  "am2", "ap2", "bm2", "bp2", "cm2", "cp2", "dm2", "dp2")
appender1 <- function(string) {
  return(equations1[match(string, names1)])
}

#for oracle models
equations <-  c(TeX("0.3Hg$*$Ni"), TeX("0.1Hg$*($Ni$-1)^2$"), 
                TeX("0.3Cd$*$As"), TeX("0.1Cd$*($As$-1)^2$"), 
                TeX("0.3Hg$*$Co"), TeX("0.1Hg$*($Co$-1)^2$"), 
                TeX("0.3Hg$*$Ni$*$Tl"), TeX("0.1Hg$*($Ni$-1)^2*$Tl"), 
                TeX("0.6Hg$*$Ni"), TeX("0.2Hg$*($Ni$-1)^2$"), 
                TeX("0.6Cd$*$As"), TeX("0.2Cd$*($As$-1)^2$"), 
                TeX("0.6Hg$*$Co"), TeX("0.2Hg$*($Co$-1)^2$"), 
                TeX("0.6Hg$*$Ni$*$Tl"), TeX("0.2Hg$*($Ni$-1)^2*$Tl"))
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
mlrmods <- read_rds("sim/mlr/mlr_mods_sm.RDS")
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
chem_mods <- read_rds("sim/mlr/chem_mods_sm.RDS")
rsq_chem <- chem_mods |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_mods))

rsq_chem |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/chem_sm_rsq.png", width = 7.5, height = 5)

#chem model output large
chem_modl <- read_rds("sim/mlr/chem_mods_lg.RDS")
rsq_cheml <- chem_modl |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_modl))

rsq_cheml |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/chem_lg_rsq.png", width = 7.5, height = 5)

###############
#oracle models
##############

#small models
oracle_mods <- read_rds("sim/oracle/oracle_mods_sm.RDS")

#extract r squared
rsquared2 <- oracle_mods |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(oracle_mods))

rsquared2 |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/oracle_sm_rsq.png", width = 7.5, height = 5)

#extract p values
keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')
pval <- oracle_mods |> 
  purrr::map_df(\(x) {
    x <- summary(x)$coefficients[,4]
    return(data.frame(pval = x[!(names(x) %in% keepnames)]))
  }) |> 
  mutate(name = names(oracle_mods)[101:1700])

pval |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender, 
                                    default = label_parsed)) +
  labs(y = "Density", x = "p-value")
ggsave("sim/figs/oracle_sm_pval_dist.png", width = 7.5, height = 5)

#extract power at alpha = 0.05
oracle_sm_power <- pval |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())


write_csv(oracle_sm_power, "sim/tables/oracle_sm_power.csv")

#large models
oracle_modl <- read_rds("sim/oracle/oracle_mods_lg.RDS")

#extract r squared
rsquared2l <- oracle_modl |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(oracle_modl))

rsquared2l |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender1, 
                                    default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
ggsave("sim/figs/oracle_sm_rsq.png", width = 7.5, height = 5)

#p-value
pval_large <- oracle_modl |> 
  purrr::map_df(\(x) {
    x <- summary(x)$coefficients[,4]
    return(data.frame(pval = x[!(names(x) %in% keepnames)]))
  }) |> 
  mutate(name = names(oracle_modl)[101:1700]) 

pval_large |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y", 
             labeller = as_labeller(appender, 
                                    default = label_parsed)) +
  labs(y = "Density", x = "p-value")
ggsave("sim/figs/oracle_lg_pval_dist.png", width = 7.5, height = 5)

#extract power at alpha = 0.05
oracle_lg_power <- pval_large |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())

write_csv(oracle_lg_power, "sim/tables/oracle_lg_power.csv")
