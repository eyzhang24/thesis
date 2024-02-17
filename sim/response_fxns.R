# create functions for various response variables

# base case, no interactions
base_case <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

am1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.35*Hg*Ni + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

am2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.7*Hg*Ni + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

ap1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.13*Hg*((Ni-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

ap2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.26*Hg*((Ni-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.35*Cd*As + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.7*Cd*As + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.125*Cd*((As-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

bp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.25*Cd*((As-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*Co + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Hg*Co + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.15*Hg*((Co-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

cp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*((Co-1)^2) +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dm1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.3*Hg*Ni*Tl + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dm2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.6*Hg*Ni*Tl + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dp1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.09*Hg*((Ni-1)^2)*Tl +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

dp2 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           0.18*Hg*((Ni-1)^2)*Tl +
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}
