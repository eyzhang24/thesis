# create functions for various response variables

# base case, no interactions
base_case <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

am1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*Ni + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

am2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Hg*Ni + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

ap1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.1*Hg*((Ni-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

ap2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.2*Hg*((Ni-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Cd*As + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Cd*As + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.1*Cd*((As-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.2*Cd*((As-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*Co + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Hg*Co + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.1*Hg*((Co-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.2*Hg*((Co-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*Ni*Tl + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Hg*Ni*Tl + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.1*Hg*((Ni-1)^2)*Tl +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.2*Hg*((Ni-1)^2)*Tl +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 2, 
                     race == 2 ~ -0.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ -0.25) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}
