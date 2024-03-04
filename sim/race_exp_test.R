library(tidyverse)

# test models for interaction b/t race and exp

# read in original data
comb_small <- read_csv("madres_data/base_data.csv")

#log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric))
# make sure to scale before modeling!!!

# race counts
# 1   2   3   4   5 
# 16  27  13  87 109 

# read in simulated data (small)
# not scaled yet
out1 <- read_rds("sim/sim_preds_sm.RDS")

# one <- out1[[1]]

# function for creating response
em2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5 + Hg, # double in group 2
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

ep2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5 + Hg) + # double in group 5
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

set.seed(1)
# test on original data
comb_resp_1 <- comb_log |> 
  mutate(across(10:19, scale)) |> 
  em2() |> 
  mutate(race = ifelse(race == 2, 1, 0)) |> # make baseline non-interaction
  mutate(across(race:smoke, as.factor))

comb_resp_2 <- comb_log |> 
  mutate(across(10:19, scale)) |> 
  ep2() |> 
  mutate(race = ifelse(race == 5, 1, 0)) |> # make baseline non-interaction
  mutate(across(race:smoke, as.factor))

# lm on original data
mod_o1 <- lm(y ~ Hg + Sb +
              I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
              race*Hg + 
              age + bmi + race + smoke, data = comb_resp_1)
mod_o2 <- lm(y ~ Hg + Sb +
               I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
               race*Hg + 
               age + bmi + race + smoke, data = comb_resp_2)

summary(mod_o1)
summary(mod_o2)

# now let's test it on all the simulated data
# numbers 1801-1900 will be larger eff size in smaller group
# numbers 2001-2100 will be larger eff size in larger group
out_small1 <- out1[1801:1900]
out_small2 <- out1[2001:2100]

# simulate response
set.seed(1)
out_resp1 <- out_small1 |> 
  map(\(x) {
    x |> 
      mutate(across(5:14, scale)) |> 
      em2() |> 
      mutate(race = ifelse(race == 2, 1, 0)) |> # make baseline non-interaction
      mutate(across(race:smoke, as.factor))
  })

set.seed(1)
out_resp2 <- out_small2 |> 
  map(\(x) {
    x |> 
      mutate(across(5:14, scale)) |> 
      ep2() |> 
      mutate(race = ifelse(race == 5, 1, 0)) |> # make baseline non-interaction
      mutate(across(race:smoke, as.factor))
  })

# look at linear models
mods1 <- out_resp1 |> 
  map(\(x) {
    lm(y ~ Hg + Sb +
         I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
         race*Hg + 
         age + bmi + race + smoke, data = x)
  })
mods2 <- out_resp2 |> 
  map(\(x) {
    lm(y ~ Hg + Sb +
         I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
         race*Hg + 
         age + bmi + race + smoke, data = x)
  })

# extract oracle model outputs
pval1 <- mods1 |> 
  map_dbl(\(x) summary(x)$coef["Hg:race1", 4])
pval2 <- mods2 |> 
  map_dbl(\(x) summary(x)$coef["Hg:race1", 4])

ggplot(NULL, aes(x = pval1)) +
  geom_histogram(binwidth = 0.05)
ggplot(NULL, aes(x = pval2)) +
  geom_histogram(binwidth = 0.05)

sum(pval1 < 0.05)/100
sum(pval2 < 0.05)/100
