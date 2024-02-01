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
simulate_data <- function(data, rho, df = 1) {
  #' data = original observed data
  #' rho = rho values from t-copula
  #' df = degrees of freedom from t-copula
  
  # simulate pseudo-observations from copula
  samp <- rCopula(nrow(data), 
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
      mutate(simulate_data(comb_log_clip, rho = rho, df = df), sim = x)
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
          2*sin(0.25*pi*Hg*Ni) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
        r11 = 
          Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
          Cd*As + 
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
          2*sin(0.25*pi*Cd*As) + 
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
          2*sin(0.25*pi*Co*Ni) + 
          age + 0.5*bmi + 
          case_when(race == 1 ~ 2, 
                    race == 2 ~ -0.5, 
                    race == 3 ~ 1, 
                    race == 4 ~ -0.25) +
          ifelse(smoke == 1, -1, 0.5) +
          rnorm(nrow(df), 0, 5), 
      )
  })


keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'smoke1')
pval <- vector(mode='list', length=1000)
rsq <- vector(mode='list', length=1000)
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
                Hg*Ni + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m11)$r.squared
    x <- summary(m11)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 600) {
    m12 <- lm(r12 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Hg*((Ni-1)^2)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m12)$r.squared
    x <- summary(m12)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 700) {
    m13 <- lm(r13 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(sin(0.25*pi*Hg*Ni)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m13)$r.squared
    x <- summary(m13)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 800) {
    m21 <- lm(r21 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                Hg*Ni + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m21)$r.squared
    x <- summary(m21)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 900) {
    m22 <- lm(r22 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(Hg*((Ni-1)^2)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m22)$r.squared
    x <- summary(m22)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  } else if(i <= 1000) {
    m23 <- lm(r23 ~ Hg + Sb +
                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                I(sin(0.25*pi*Hg*Ni)) + 
                age + bmi + race + smoke, data = df2)
    rsq[[i]] <- summary(m23)$r.squared
    x <- summary(m23)$coefficients[,4]
    pval[[i]] <- x[!(names(x) %in% keepnames)]
  }
}

#look at p-values
dfpval <- do.call(rbind.data.frame, pval)
names(dfpval) <- "pval"
dfpval$index <- 1:1000
dfpval$pval <- ifelse(dfpval$pval == -99, NA, dfpval$pval)
dfpval$group <- with(dfpval, factor(ceiling(index/100)))
dfpval |> 
  ggplot(aes(x = pval)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
ggsave("pvalplot.png", width = 6, height = 4)

#look at r-squared
dfrsq <- do.call(rbind.data.frame, rsq)
names(dfrsq) <- "rsq"
dfrsq$index <- 1:1000
dfrsq$group <- with(dfrsq, factor(ceiling(index/100)))
dfrsq |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~group, scales = "free_y")
ggsave("rsqplot.png", width = 6, height = 4)



