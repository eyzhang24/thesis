library(tidyverse)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

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

#mlr model output
mlrmods <- read_rds("sim/mlr_mods_sm.RDS")
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
  facet_wrap(~name, scales = "free_y")

#chem model output
chem_oracle <- read_rds("sim/chem_oracle_sm.RDS")
rsq_chem <- chem_oracle |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_oracle))

rsq_chem |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y")

chem_oraclel <- read_rds("sim/chem_oracle_lg.RDS")
rsq_cheml <- chem_oraclel |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_oraclel))

rsq_cheml |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y")

oracle_mods <- read_rds("sim/oracle_mods_sm.RDS")
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
  facet_wrap(~name, scales = "free_y")

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
  facet_wrap(~name, scales = "free_y")

pval |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())

oracle_modl <- read_rds("sim/oracle_mods_lg.RDS")
pval_large <- oracle_modl |> 
  purrr::map_df(\(x) {
    x <- summary(x)$coefficients[,4]
    return(data.frame(pval = x[!(names(x) %in% keepnames)]))
  }) |> 
  mutate(name = names(oracle_modl)[101:1700])

pval_large |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y")

pval_large |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())

### test p-val for racexexp
pval <- oracle_mods |> 
  purrr::map_df(\(x) {
    x <- summary(x)$coefficients[,4]
    return(data.frame(pval = x[!(names(x) %in% keepnames)]))
  }) |> 
  mutate(name = names(oracle_mods))

pval |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y")

pval |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())

oracle_modl <- read_rds("sim/oracle_mods_lg_re.RDS")
pval_large <- oracle_modl |> 
  purrr::map_df(\(x) {
    x <- summary(x)$coefficients[,4]
    return(data.frame(pval = x[!(names(x) %in% keepnames)]))
  }) |> 
  mutate(name = names(oracle_modl))

pval_large |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~name, scales = "free_y")

pval_large |> 
  group_by(name) |> 
  summarize(power = sum(pval < 0.05)/n())
