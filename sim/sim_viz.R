library(tidyverse)
library(latex2exp)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

#############
#predictors
#############

comb_small <- read_csv("madres_data/base_data.csv")

#log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric))

#look at predictor data distributions
out1 <- read_rds("sim/sim_preds_sm.RDS")

comb_sim1 <- bind_rows(out1)

#density plots
comb_sim1 |> 
  mutate(sim = as.factor(sim)) |> 
  select(3:15) |> 
  pivot_longer(cols = 1:12) |>
  mutate(name = factor(name, levels = names(comb_sim1)[3:14])) |> 
  ggplot(aes(x = value, group = sim)) +
  geom_line(stat = "density", color = "grey10", alpha = 0.01) + 
  geom_density(
    data = comb_log |> select(8:19) |> pivot_longer(cols = 1:12),
    mapping = aes(x = value), 
    color = "deepskyblue", linewidth = 0.75, inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free")

#ggpairs
comb_sim1_p <- comb_sim1 |> 
  filter(sim <= 200) |> 
  select(3:14) |>
  mutate(sim = 0.01) |> 
  bind_rows(mutate(select(comb_log, 8:19), 
                   sim = 1))

GGally::ggpairs(comb_sim1_p, aes(alpha = sim, color = sim), columns = 1:12, 
                upper = 'blank', diag = 'blank')

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
