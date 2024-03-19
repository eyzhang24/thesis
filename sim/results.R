library(tidyverse)

# plot pips bkmr small
ksm_pips <- read_csv("sim/bkmr_sm/pips.csv")

# boxplot of pip values
ksm_pips |> 
  ggplot(aes(x = variable, y = PIP)) +
  geom_boxplot() + 
  facet_wrap(~case)

# proportion > 0.5
ksm_pips |> 
  mutate(imp = PIP > 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)) |> 
  ggplot(aes(x = variable, y = sensitivity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~case)

# plot univariate relationships bkmr small


# plot bivariate relationships bkmr small
ksm_biv <- read_csv("sim/bkmr_sm/biv_expresp.csv")

# plot v1 vs. v2, 
ksm_biv |> 
  filter(case %in% 2:5) |> 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x") +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-6, 6) 

# int vs. rest significant visualization
ksm_ints <- read_csv("sim/bkmr_sm/ints.csv")

ksm_ints <- ksm_ints |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

summary(ksm_ints$signif)

ksm_ints |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = variable, y = sensitivity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~case, scales = "free_x")

# one vs. other significant visualization
ksm_intb <- read_csv("sim/bkmr_sm/int_bivar.csv")

ksm_intb <- ksm_intb |> 
  mutate(cond = paste0(z1, "+", z2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

table(ksm_intb$signif)

ksm_intb |> 
  group_by(case) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = "", y = sensitivity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~case)

# one vs. two others significant visualization
ksm_intt <- read_csv("sim/bkmr_sm/int_trivar.csv")
# glimpse(ksm_intt)

ksm_intt <- ksm_intt |> 
  mutate(cond = paste0(variable, " by ", fixedat1, "+", fixedat2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

ksm_intt |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |>
  ggplot(aes(x = cond, y = sensitivity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~case)

# bivar and trivar together
int_comb <- bind_rows(
  select(ksm_intb, cond, trial, case, signif), 
  select(ksm_intt, cond, trial, case, signif)
)

int_comb |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(cond, sensitivity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~case, scales = "free_x")

