library(tidyverse)
library(latex2exp)

col2hex <- function(x, alpha = FALSE) {
  args <- as.data.frame(t(col2rgb(x, alpha = alpha)))
  args <- c(args, list(names = x, maxColorValue = 255))
  do.call(rgb, args)
}

# set up ------------------------------------------------------------------

# set theme for plots
theme_set(theme_light())
theme_update(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

getsign_re <- function(chem) {
  return(chem %in% c("Hg", "Ni", "Sb", "Sn"))
}

which_race <- function(race) {
  case_when(
    race == 1 ~ "Non-Hisp. white", 
    race == 2 ~ "Non-Hisp. black", 
    race == 3 ~ "Non-Hisp. other", 
    race == 4 ~ "Hispanic born in US", 
    race == 5 ~ "Hispanic born outside US",
    race == 6 ~ "Collapsed non-Hisp."
  )
}

size_race <- function(race) {
  case_when(
    race == 1 ~ "16", 
    race == 2 ~ "27", 
    race == 3 ~ "13", 
    race == 4 ~ "87", 
    race == 5 ~ "109",
    race == 6 ~ "56"
  )
}

# naive -------------------------------------------------------------------

nsmre_pval <- read_csv("sim/_mlr/pvalsmre.csv")

nlgre_pval <- read_csv("sim/_mlr/pvallgre.csv")

# p-values small
nsmre_pvalc <- nsmre_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = getsign_re(var), 
         mod = "Naive MLR")

# p-values large
nlgre_pvalc <- nlgre_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = getsign_re(var), 
         mod = "Naive MLR")

# save sensitivity as table
naive_allre <- bind_rows(
  mutate(nsmre_pvalc, size = "Small", mod = "Naive MLR"),
  mutate(nlgre_pvalc, size = "Large", mod = "Naive MLR")
) |> 
  group_by(mod, case, size, var, sign) |> 
  summarize(sensitivity = sum(p < 0.05)/n())
write_csv(naive_allre, "index/data/naive_re_sens.csv")


# oracle ------------------------------------------------------------------

# small
osmre_pval <- read_csv("sim/_oracle/pvalsmre.csv")

# univariate p-values
osmre_pvalc <- osmre_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = TRUE)

# interaction terms
keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')

osmre_pvali <- osmre_pval |> 
  filter(!(var %in% keepnames)) |> 
  mutate(trial = case, 
         case = ceiling(case/100))

# combine univariate and interactions
osmre_comb <- bind_rows(
  mutate(osmre_pvali, variable = "Int"), 
  mutate(osmre_pvalc, variable = var)
)

# large
olgre_pval <- read_csv("sim/_oracle/pvallgre.csv")

# univariate p-values
olgre_pvalc <- olgre_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = TRUE)

# interaction terms
keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')

olgre_pvali <- olgre_pval |> 
  filter(!(var %in% keepnames)) |> 
  mutate(trial = case, 
         case = ceiling(case/100))

# combine univariate and interactions
olgre_comb <- bind_rows(
  mutate(olgre_pvali, variable = "Int"), 
  mutate(olgre_pvalc, variable = var)
)

# combine both oracle outputs
oracle_re <- bind_rows(
  mutate(osmre_comb, size = "Small"), 
  mutate(olgre_comb, size = "Large")
) |> 
  select(size, case, trial, variable, p) |> 
  arrange(case, trial, desc(size))

write_csv(oracle_re, "index/data/oracle_re_sens.csv")


# bkmr --------------------------------------------------------------------

library(intervals)
# critical value of 2.57583 for bonferroni corrected 95% CI

# data <- ksmre_ints |>
#   mutate(lower = est - 2.57583*sd,
#          upper = est + 2.57583*sd,
#          sign = (case <= 2 & race == 2) | (case >= 3 & race == 5)) |>
#   filter(trial == 301)

# function for checking overlap of confidence intervals
check_overlap <- function(data) {
  num_na <- sum(is.na(data$lower))
  sign <- data$sign 
  unsign <- !(data$sign)
  sign_intervals <- Intervals(na.omit(data[sign, c("lower", "upper")]))
  unsign_intervals <- Intervals(na.omit(data[unsign, c("lower", "upper")]))
  return(data.frame(
    num_na = num_na,
    sig_overlaps = length(interval_overlap(sign_intervals, unsign_intervals)[[1]]), 
    insig_overlaps = (length(unlist(interval_overlap(unsign_intervals, unsign_intervals))) == 
                        (length(unsign_intervals)/2)^2)
  ))
}

# check_overlap(data)

# detection of interactions 

# small
ksmre_ints <- read_csv("sim/re/ksmre_ints.csv")

# uncollapsed
ksmre_det1 <- ksmre_ints |> 
  filter(race != 6) |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 2) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# check number of models where at least one was unable to be fitted
table(ksmre_det1$num_na != 0)

# use collapsed
ksmre_det2 <- ksmre_ints |> 
  filter(race %in% c(4, 5, 6)) |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 6) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# large
klgre_ints <- read_csv("sim/re/klgre_ints.csv")

klgre_det <- klgre_ints |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 2) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# save the trials that produced significant results
klgre_det |> 
  filter(sig_overlaps != 4) |> 
  write_csv("sim/re/trials_klg_sig.csv")

# get sensitivity
ksmre_sens1 <- ksmre_det1 |> 
  group_by(case) |> 
  summarize(sens = sum((sig_overlaps/(4-num_na)) != 1 & insig_overlaps)/n()) |> 
  mutate(size = "Small uncollapsed")
ksmre_sens2 <- ksmre_det2 |> 
  group_by(case) |> 
  summarize(sens = sum((sig_overlaps/(2-num_na)) != 1 & insig_overlaps)/n()) |> 
  mutate(size = "Small collapsed")
# no trials that had a significant non-overlap AlSO had non-significant non-overlap
klgre_sens <- klgre_det |> 
  group_by(case) |> 
  summarize(sens = sum(sig_overlaps != 4 & insig_overlaps)/n()) |> 
  mutate(size = "Large")

# combine sensitivities together
bkmrre_senscomb <- bind_rows(ksmre_sens1, ksmre_sens2, klgre_sens) 
write_csv(bkmrre_senscomb, "index/data/bkmr_re_sens.csv")

# pip's
ksmre_pips <- read_csv("sim/re/ksmre_pips.csv")
klgre_pips <- read_csv("sim/re/klgre_pips.csv")

# summarize pip's
ksmre_pips2 <- ksmre_pips |> 
  filter(is.na(variable) | variable == "Hg") |> 
  group_by(case, race) |> 
  summarize(num_na = sum(is.na(variable))/100, 
            hg_detect = sum(PIP >= 0.5, na.rm = T)/100)

klgre_pips2 <- klgre_pips |> 
  filter(is.na(variable) | variable == "Hg") |> 
  group_by(case, race) |> 
  summarize(num_na = sum(is.na(variable))/n(), 
            hg_detect = sum(PIP >= 0.5, na.rm = T)/n())

bkmrre_pipscomb <- bind_rows(
  mutate(ksmre_pips2, size = "Small"), 
  mutate(klgre_pips2, size = "Large")
)  |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race))
write_csv(bkmrre_pipscomb, "index/data/bkmr_re_pips.csv")


# bsr ---------------------------------------------------------------------

# pip's

# pip's
ssmre_pips <- read_csv("sim/re/ssmre_pips.csv")
slgre_pips <- read_csv("sim/re/slgre_pips.csv")

# summarize pip's
ssmre_pips2 <- ssmre_pips |> 
  filter(is.na(variable) | variable == "Hg") |> 
  group_by(case, race) |> 
  summarize(num_na = sum(is.na(variable))/100, 
            hg_detect = sum(PIP >= 0.5, na.rm = T)/100)

slgre_pips2 <- slgre_pips |> 
  filter(is.na(variable) | variable == "Hg") |> 
  group_by(case, race) |> 
  summarize(num_na = sum(is.na(variable))/n(), 
            hg_detect = sum(PIP >= 0.5, na.rm = T)/n())

bsrre_pipscomb <- bind_rows(
  mutate(ssmre_pips2, size = "Small"), 
  mutate(slgre_pips2, size = "Large")
)  |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race))
write_csv(bsrre_pipscomb, "index/data/bsr_re_pips.csv")


# exposure response relationship ------------------------------------------


#### BKMR #####

# small
ksmre_expresp <- read_csv("sim/re/ksmre_expresp.csv")
# large
klgre_expresp <- read_csv("sim/re/klgre_expresp.csv")
# combine
bkmr_re_expresp <- bind_rows(
  mutate(ksmre_expresp, size = "Small"), 
  mutate(klgre_expresp, size = "Large")
) |> 
  filter(race != 6) 

# first two cases
kre12 <- bkmr_re_expresp |> 
  filter(case <= 2) |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race), 
         size = factor(size, levels = c("Small", "Large"))) 

k12 <- kre12 |>
  mutate(
    race2 = case_when(race == 2 ~ 5,
                      race == 5 ~ 2,
                      .default = race),
    race2 = factor(
      race2,
      levels = c(1, 2, 3, 4, 5),
      labels = c(1, 5, 3, 4, 2)
    ), 
    case = factor(ifelse(case == 1, "Lower", "Higher"), 
                  levels = c("Lower", "Higher"))
  ) |> 
  ggplot(aes(z1, est, color = race2)) +
  geom_line(aes(group = interaction(trial, race2)), alpha = 0.15) +
  ylim(-8, 8) + 
  facet_grid(size ~ case) +
  scale_color_manual(values = rev(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")), 
                     breaks = c(1, 2, 3, 4, 5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Hg", y = "Estimate", subtitle = "BKMR, interaction in Non-Hispanic black")

# second two cases
kre34 <- bkmr_re_expresp |> 
  filter(case > 2) |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race), 
         size = factor(size, levels = c("Small", "Large"))) 

k34 <- kre34 |>
  mutate(
    race2 = factor(
      race,
      levels = c(1, 2, 3, 4, 5),
    ), 
    case = factor(ifelse(case == 3, "Lower", "Higher"), 
                  levels = c("Lower", "Higher"))
  ) |> 
  ggplot(aes(z1, est, color = race2)) +
  geom_line(aes(group = interaction(trial, race2)), alpha = 0.15) +
  ylim(-8, 8) + 
  facet_grid(size ~ case) +
  scale_color_manual(values = rev(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")), 
                     breaks = c(1, 2, 3, 4, 5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Hg", y = "Estimate", subtitle = "BKMR, interaction in Hispanic born outside US")


#### BSR #####

# small
ssmre_expresp <- read_csv("sim/re/ssmre_expresp.csv")
# large
slgre_expresp <- read_csv("sim/re/slgre_expresp.csv")
# combine
bsr_re_expresp <- bind_rows(
  mutate(ssmre_expresp, size = "Small"), 
  mutate(slgre_expresp, size = "Large")
) |> 
  filter(race != 6) 

# first two cases
sre12 <- bsr_re_expresp |> 
  filter(case <= 2) |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race), 
         size = factor(size, levels = c("Small", "Large"))) 

s12 <- sre12 |>
  mutate(
    race2 = case_when(race == 2 ~ 5,
                      race == 5 ~ 2,
                      .default = race),
    race2 = factor(
      race2,
      levels = c(1, 2, 3, 4, 5),
      labels = c(1, 5, 3, 4, 2)
    ), 
    case = factor(ifelse(case == 1, "Lower", "Higher"), 
                  levels = c("Lower", "Higher"))
  ) |> 
  ggplot(aes(j1val, est, color = race2)) +
  geom_line(aes(group = interaction(trial, race2)), alpha = 0.15) +
  ylim(-6, 12) +
  facet_grid(size ~ case) +
  scale_color_manual(values = rev(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")), 
                     breaks = c(1, 2, 3, 4, 5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Hg", y = "Estimate", subtitle = "BSR, interaction in Non-Hispanic black")

# second two cases
sre34 <- bsr_re_expresp |> 
  filter(case > 2) |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race), 
         size = factor(size, levels = c("Small", "Large"))) 

s34 <- sre34 |>
  mutate(
    race2 = factor(
      race,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Non-Hisp. white", "Non-Hisp. other", "Non-Hisp. black", 
                 "Hispanic born in US", "Hispanic born outside US")
    ), 
    case = factor(ifelse(case == 3, "Lower", "Higher"), 
                  levels = c("Lower", "Higher"))
  ) |> 
  ggplot(aes(j1val, est, color = race2)) +
  geom_line(aes(group = interaction(trial, race2)), alpha = 0.15) +
  ylim(-6, 12) +
  facet_grid(size ~ case) +
  scale_color_manual(values = rev(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"))) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Hg", y = "Estimate", color = "Race", 
       subtitle = "BSR, interaction in Hispanic born outside US") +
  theme(legend.position = "bottom")

re_expresp <- cowplot::plot_grid(
  k12 + theme(legend.position = "none"),
  s12 + theme(legend.position = "none"), 
  k34 + theme(legend.position = "none"),
  s34 + theme(legend.position = "none"),
  nrow = 2
)

reerlegend <- cowplot::get_legend(s34)
cowplot::plot_grid(re_expresp, reerlegend, rel_heights = c(5, 0.35), ncol = 1)
ggsave("index/figures/ch4_re_expresp.png", width = 7.75, height = 6)

# run time ----------------------------------------------------------------

nsmre_times <- read_rds("sim/_mlr/mlr_mods_sm_re_times.RDS")
nlgre_times <- read_rds("sim/_mlr/mlr_mods_lg_re_times.RDS")
osmre_times <- read_rds("sim/_oracle/oracle_mods_sm_re_times.RDS")
olgre_times <- read_rds("sim/_oracle/oracle_mods_lg_re_times.RDS")
ksmre_times <- read_rds("sim/re/ksmre_times.RDS")
klgre_times <- read_rds("sim/re/klgre_times.RDS")

# change list to dataframe
nsmre_times2 <- unlist(nsmre_times) |> 
  as.data.frame(row.names = names(nsmre_times)) |> 
  rename(time = 1) |> 
  rownames_to_column(var = "case1") |> 
  mutate(case = match(case1, unique(case1)), 
         trial = 1:400, 
         time = as.difftime(time, units = "secs")) |> 
  select(-case1)
nlgre_times2 <- unlist(nlgre_times) |> 
  as.data.frame(row.names = names(nlgre_times)) |> 
  rename(time = 1) |> 
  rownames_to_column(var = "case1") |> 
  mutate(case = match(case1, unique(case1)), 
         trial = 1:400, 
         time = as.difftime(time, units = "secs")) |> 
  select(-case1)
osmre_times2 <- unlist(osmre_times) |> 
  as.data.frame(row.names = names(osmre_times)) |> 
  rename(time = 1) |> 
  rownames_to_column(var = "case1") |> 
  mutate(case = match(case1, unique(case1)), 
         trial = 1:400, 
         time = as.difftime(time, units = "secs")) |> 
  select(-case1)
olgre_times2 <- unlist(olgre_times) |> 
  as.data.frame(row.names = names(olgre_times)) |> 
  rename(time = 1) |> 
  rownames_to_column(var = "case1") |> 
  mutate(case = match(case1, unique(case1)), 
         trial = 1:400, 
         time = as.difftime(time, units = "secs")) |> 
  select(-case1)

re_times <- bind_rows(
  mutate(nsmre_times2, mod = "Naive", size = "Small", race = NA), 
  mutate(nlgre_times2, mod = "Naive", size = "Large", race = NA), 
  mutate(osmre_times2, mod = "Oracle", size = "Small", race = NA), 
  mutate(olgre_times2, mod = "Oracle", size = "Large", race = NA), 
  mutate(ksmre_times, mod = "BKMR", size = "Small", case = as.numeric(case)), 
  mutate(klgre_times, mod = "BKMR", size = "Large", case = as.numeric(case))
)

# function to format difftime objects with appropriate units
format_difftime <- function(difftime_obj) {
  value <- as.numeric(difftime_obj)
  # determine appropriate units
  if (abs(value) < 60) {
    duration <- as.numeric(difftime_obj, units = "secs")
    result <- sprintf("%.2g %s", duration, "s")
  } else if (abs(value) < 3600) {
    duration <- as.numeric(difftime_obj, units = "mins")
    result <- sprintf("%.2f %s", duration, "m")
  } else if (abs(value) < 86400) {
    duration <- as.numeric(difftime_obj, units = "hours")
    result <- sprintf("%.2f %s", duration, "h")
  } else {
    duration <- as.numeric(difftime_obj, units = "days")
    result <- sprintf("%.2f %s", duration, "d")
  }
  return(result)
}

# format into table of mean run-times
retime_table <- re_times |> 
  mutate(mod = factor(mod, levels = c("Naive", "Oracle", "BKMR")), 
         size = factor(size, levels = c("Small", "Large")), 
         case = as.numeric(case), 
         scenario = ifelse(case %in% c(1, 2), 
                           "Smaller race cat.", 
                           "Larger race cat.")) |> 
  group_by(mod, size, race) |> 
  summarize(mean_time = mean(time)) |> 
  ungroup() |> 
  rowwise() |> 
  mutate(mean_time = format_difftime(mean_time)) |> 
  arrange(size) |> 
  pivot_wider(names_from = size, values_from = mean_time, 
              names_sep = " ") |> 
  mutate(race = which_race(race))
write_csv(retime_table, "index/data/time_re.csv")

## OLD ##

# # first, label exposure-response with whether it was significant
# klgre_pipstojoin <- klgre_pips |> 
#   filter(is.na(variable) | variable == "Hg") |> 
#   mutate(sign = ifelse(is.na(PIP), FALSE, PIP >= 0.5)) |> 
#   select(case, trial, race, sign)
# klgre_expresp2 <- klgre_expresp |>
#   left_join(klgre_pipstojoin, by = c("case", "trial", "race")) |> 
#   mutate(which_race = which_race(race), 
#          size_race = size_race(race))
#"#A54657", #B084CC

# klgre_expresp |> 
#   mutate(race = factor(race),
#          which_race = which_race(race), 
#          size_race = size_race(race)) |> 
#   ggplot(aes(z1, est, color = race)) +
#   geom_line(aes(group = interaction(trial, race)), alpha = 0.2) +
#   ylim(-8, 8) +
#   facet_wrap(~case)  + 
#   guides(color = guide_legend(override.aes = list(alpha = 1)))

