library(tidyverse)
library(latex2exp)

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
ksmre_sens1 <- ksmre_ints |> 
  filter(race != 6) |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 2) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# use collapsed
ksmre_sens2 <- ksmre_ints |> 
  filter(race %in% c(4, 5, 6)) |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 6) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# large
klgre_ints <- read_csv("sim/re/klgre_ints.csv")

klgre_sens <- klgre_ints |> 
  mutate(lower = est - 2.57583*sd, 
         upper = est + 2.57583*sd, 
         sign = (case <= 2 & race == 2) | (case >= 3 & race == 5)) |> 
  group_by(case, trial) |> 
  group_modify(~ check_overlap(.))

# save the trials that produced significant results
klgre_sens |> 
  filter(sig_overlaps != 4) |> 
  write_csv("sim/re/trials_klg_sig.csv")


# pip's


# exposure-response relationship 
  


# bsr ---------------------------------------------------------------------



# run time ----------------------------------------------------------------


