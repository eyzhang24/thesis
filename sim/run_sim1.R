library(tidyverse)
library(rslurm)
library(bkmr)
library(NLinteraction)
library(copula)

# read in data
comb_small <- read_csv("madres_data/base_data.csv")

# log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log))
comb_log_clip <- comb_log |> 
  select(5:19)

########
# simulate data
########
# read back in final copula
cfit_t <- read_rds("madres_data/tcop1.RDS")

# get rho and degrees of freedom
rho <- coef(cfit_t)[1:105]
df <- coef(cfit_t)[106]

# create function for simulation
simulate_data <- function(data, sampsize, rho, df = 1) {
  #' data = original observed data
  #' sampsize = size of simulated dataset
  #' rho = rho values from t-copula
  #' df = degrees of freedom from t-copula
  
  # simulate pseudo-observations from copula
  samp <- rCopula(sampsize, 
                  tCopula(rho, dim = ncol(data), dispstr = "un", df = df))
  # transform pseudo-observations to observed marginal distributions
  sampt <- 1:ncol(data) |> 
    purrr::map_dfc(
      \(x) {
        df <- data.frame(quantile(data[[x]], probs = samp[,x]), 
                         row.names = NULL)
        names(df) <- names(data)[x]
        return(df)
      }
    )
  return(sampt)
}

# create function to run on hpc
run_sim <- function() {
  set.seed(0)
  out <- 1:2500 |> 
    purrr::map(\(x) {
      mutate(simulate_data(comb_log_clip, sampsize = nrow(comb_log_clip), 
                           rho = rho, df = df), sim = x)
    })
  return(out)
}

# send job to hpc
sjob <- slurm_call(run_sim, 
                   global_objects = c('comb_log_clip', 'rho', 'df', 
                                      'simulate_data'),
                   jobname = 'sim_data')

# check on status of job
get_job_status(sjob)$completed
system('squeue --me')

# get output
out <- get_slurm_out(sjob)
write_rds(out, "sim_preds.RDS")

# read output back in
out <- read_rds("sim_preds.RDS")

#########
# test two way interactions
#########

# create response variable
set.seed(0)
out2 <- out |> 
  purrr::map(\(df) {
    df |> 
      mutate(across(age:Sn, ~c(scale(.)))) |> 
      mutate(across(mom_site:smoke, ~as.factor(round(.)))) |> 
      mutate(
        r00 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r01 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.6*Hg*Ni + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r02 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.2*Hg*((Ni-1)^2) +
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r03 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          1.5*sin(0.25*pi*Hg*Ni) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r11 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.6*Cd*As + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r12 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.2*Cd*((As-1)^2) +
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r13 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          1.5*sin(0.25*pi*Cd*As) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r21 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.6*Co*Ni + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r22 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.2*Co*((Ni-1)^2) +
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r23 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          1.5*sin(0.25*pi*Co*Ni) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
      )
  })

write_rds(out2, "sim_preds_resp.RDS")
out2 <- read_rds("sim_preds_resp.RDS")

# extract oracle model outputs
keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'smoke1', 'Cd', 'As', 'Co')
pval <- vector(mode='double', length=1000)
rsq <- vector(mode='double', length=1000)
for(i in 1:1000) {
  df2 <- out2[[i]]
  if(i <= 100) {
    m00 <- lm(r00 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m00)$r.squared
    x <- summary(m00)$coefficients[,4]
    pval[[i]] <- -99
  } else if(i <= 200) {
    m01 <- lm(r01 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                Hg*Ni + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m01)$r.squared
    x <- summary(m01)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 300) {
    m02 <- lm(r02 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Hg*((Ni-1)^2)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m02)$r.squared
    x <- summary(m02)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 400) {
    m03 <- lm(r03 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(sin(0.25*pi*Hg*Ni)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m03)$r.squared
    x <- summary(m03)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 500) {
    m11 <- lm(r11 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                Cd*As + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m11)$r.squared
    x <- summary(m11)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 600) {
    m12 <- lm(r12 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Cd*((As-1)^2)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m12)$r.squared
    x <- summary(m12)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 700) {
    m13 <- lm(r13 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(sin(0.25*pi*Cd*As)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m13)$r.squared
    x <- summary(m13)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 800) {
    m21 <- lm(r21 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                Co*Ni + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m21)$r.squared
    x <- summary(m21)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 900) {
    m22 <- lm(r22 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Co*((Ni-1)^2)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m22)$r.squared
    x <- summary(m22)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 1000) {
    m23 <- lm(r23 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(sin(0.25*pi*Co*Ni)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m23)$r.squared
    x <- summary(m23)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  }
}

#look at p-values
dfpval <- data.frame(pval)
dfpval$index <- 1:1000
dfpval$pval <- ifelse(dfpval$pval == -99, NA, dfpval$pval)
dfpval$group <- with(dfpval, factor(ceiling(index/100)))
dfpval |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
ggsave("pvalplot.png", width = 6, height = 4)
dfpval |> 
  group_by(group) |> 
  summarize(sum(pval < 0.05)/n())

#look at r-squared
dfrsq <- data.frame(rsq)
dfrsq$index <- 1:1000
dfrsq$group <- with(dfrsq, factor(ceiling(index/100)))
dfrsq |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
ggsave("rsqplot.png", width = 6, height = 4)

############
# test three-way and race
############

set.seed(0)
out3 <- out |> 
  purrr::map(\(df) {
    df |> 
      mutate(across(age:Sn, ~c(scale(.)))) |> 
      mutate(across(mom_site:smoke, ~as.factor(round(.)))) |> 
      mutate(
        r00 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r41 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.6*Hg*Ni*Tl + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r42 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          0.2*Hg*((Ni-1)^2)*Tl +
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r51 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          ifelse(race == 2, 2*Hg, 0) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r52 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          ifelse(race == 4, 2*Hg, 0) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5)
      )
  })

# extract oracle model outputs
keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'smoke1', 'Cd', 'As', 'Co')
keepnames2 <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
                'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
                'race2', 'race3', 'race4', 'smoke1', 'Cd', 'As', 'Co', 
                'Hg:race3', 'Hg:race4')
keepnames3 <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
                'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
                'race2', 'race3', 'race4', 'smoke1', 'Cd', 'As', 'Co', 
                'Hg:race2', 'Hg:race3')
pval <- vector(mode='double', length=500)
rsq <- vector(mode='double', length=500)
for(i in 1:500) {
  df2 <- out3[[i]]
  if(i <= 100) {
    m00 <- lm(r00 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m00)$r.squared
    x <- summary(m00)$coefficients[,4]
    pval[[i]] <- -99
  } else if(i <= 200) {
    m41 <- lm(r41 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                Hg:Ni:Tl + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m41)$r.squared
    x <- summary(m41)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 300) {
    m42 <- lm(r42 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Hg*((Ni-1)^2)*Tl) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m42)$r.squared
    x <- summary(m42)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 400) {
    m51 <- lm(r51 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                race*Hg + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m51)$r.squared
    x <- summary(m51)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames2)]
  } else if(i <= 500) {
    m52 <- lm(r52 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                race*Hg + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m52)$r.squared
    x <- summary(m52)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames3)]
  } 
}

#look at p-values
dfpval <- data.frame(pval)
dfpval$index <- 1:500
dfpval$pval <- ifelse(dfpval$pval == -99, NA, dfpval$pval)
dfpval$group <- with(dfpval, factor(ceiling(index/100)))
dfpval |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
dfpval |> 
  group_by(group) |> 
  summarize(sum(pval < 0.05)/n())

#look at r-squared
dfrsq <- data.frame(rsq)
dfrsq$index <- 1:500
dfrsq$group <- with(dfrsq, factor(ceiling(index/100)))
dfrsq |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
