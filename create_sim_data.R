library(tidyverse)
library(copula)
library(latex2exp)

# read in data
comb <- read_rds("nhanes_data/processed_data.RDS")


# fit copula --------------------------------------------------------------

# log-transform target data
# all POP's, serum cotinine, and telomere
comb_log <- comb |> 
  mutate(across(c(LBX074LA:LBXF08LA, LBXCOT, TELOMEAN), log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric)) |> 
  #re-order variables
  select(SEQN:LBXF08LA, #ID and pop's
         BMXBMI:LBXBAPCT, #bio
         RIDAGEYR:DMDEDUC2, #demo, gender is binary, re and educ multi
         TELOMEAN) |> 
  #gender to 0's and 1's: 0 male, 1 female
  #should really be called sex
  mutate(RIAGENDR = RIAGENDR - 1)

# save a copy
write_csv(comb_log, "nhanes_data/log_processed_data.RDS")

# create pseudo observations
u <- pobs(comb_log[, 2:29])
# jitter gender uniformly 
prop_gender0 <- 1 - mean(comb_log$RIAGENDR)
set.seed(0)
u_gender <- comb_log$RIAGENDR |> 
  map_dbl(\(x) {
    ifelse(x == 0, runif(1, 0, prop_gender0), runif(1, prop_gender0, 1))
  })
u[, 28] <- u_gender

# check correlation matrix
# cor(comb_log[, c(2:28)], method = "spearman")
# cor(u, method = "spearman")

#fit copulas
cfit_gaus <- fitCopula(normalCopula(dim = 28, dispstr = "un"), u)
cfit_t <- fitCopula(tCopula(dim = 28, dispstr = "un", df.fixed = FALSE), u)
cfit_t2 <- fitCopula(tCopula(dim = 28, dispstr = "un", df = 4, df.fixed = TRUE), u)
cfit_t3 <- fitCopula(tCopula(dim = 28, dispstr = "un", df = 10, df.fixed = TRUE), u)
cfit_gum <- fitCopula(gumbelCopula(4, dim = 28), u)
cfit_frank <- fitCopula(frankCopula(4, dim = 28), u)
cfit_clay <- fitCopula(claytonCopula(4, dim = 28), u)
cfit_joe <- fitCopula(joeCopula(4, dim = 28), u)
cfit_gum2 <- fitCopula(gumbelCopula(2, dim = 28), u)
cfit_frank2 <- fitCopula(frankCopula(2, dim = 28), u)
cfit_clay2 <- fitCopula(claytonCopula(2, dim = 28), u)
cfit_joe2 <- fitCopula(joeCopula(2, dim = 28), u)

#evaluate fit
aic_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, cfit_clay, cfit_joe, 
                          cfit_gum2, cfit_frank2, cfit_clay2, cfit_joe2), AIC)
names(aic_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", "cfit_clay", "cfit_joe", 
                       "cfit_gum2", "cfit_frank2", "cfit_clay2", "cfit_joe2")
sort(aic_values)

lik_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, cfit_clay, cfit_joe, 
                          cfit_gum2, cfit_frank2, cfit_clay2, cfit_joe2), logLik)
names(lik_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", "cfit_clay", "cfit_joe", 
                       "cfit_gum2", "cfit_frank2", "cfit_clay2", "cfit_joe2")
sort(lik_values)

#guassian copula performs best, proceed with this
write_rds(cfit_t, "nhanes_sim/tcop.RDS")


# sim predictor data ------------------------------------------------------

# read back in
cfit_t <- read_rds("nhanes_sim/tcop.RDS")
comb_log <- read_csv("nhanes_data/log_processed_data.RDS")

# create smaller version of dataset, with only columns for simulation
comb_log_clip <- comb_log |> 
  select(LBX074LA:RIAGENDR)

# get rho and degrees of freedom
rho <- coef(cfit_t)[1:378]
df <- coef(cfit_t)[379]

# create function for simulation
simulate_data <- function(data, n, rho, df, prop_sex, prop_race, prop_edu) {
  #'data = original observed data
  #'n = sample size
  #'rho = rho values from t-copula
  #'df = degrees of freedom from t-copula
  #'prop_sex = proportion female from observed dataset
  #'prop_race = table with race/eth values
  #'prop_edu = table with education level values
  
  #simulate pseudo-observations from copula
  samp <- rCopula(n, tCopula(rho, dim = ncol(data), dispstr = "un", df = df))
  #transform pseudo-observations to observed marginal distributions
  sampt <- 1:ncol(data) |> 
    purrr::map_dfc(
      \(x) {
        if(names(data)[x] == "RIAGENDR") {
          df <- data.frame(ifelse(samp[,x] < prop_sex, 0, 1), 
                           row.names = NULL)
        } else {
          df <- data.frame(quantile(data[[x]], probs = samp[,x]), 
                           row.names = NULL)
        }
        names(df) <- names(data)[x]
        return(df)
      }
    ) |> 
    mutate(RIDRETH1 = sample(x = names(prop_race), prob = prop_race,
                             size = n, replace = T), 
           DMDEDUC2 = sample(x = names(prop_edu), prob = prop_edu,
                             size = n, replace = T)) 
  return(sampt)
}

# reproducibility
set.seed(0)

# create 2100 size 250 datasets
out <- 1:2100 |> 
  purrr::map(\(x) {
    mutate(simulate_data(comb_log_clip, 
                         n = 250, 
                         rho = rho, df = df, 
                         prop_sex = 1-mean(comb_log_clip$RIAGENDR), 
                         prop_race = table(comb_log$RIDRETH1), 
                         prop_edu = table(comb_log$DMDEDUC2)), 
           RIDRETH1 = as.numeric(RIDRETH1), 
           DMDEDUC2 = as.numeric(DMDEDUC2), 
           sim = x) 
  })

write_rds(out, "nhanes_sim/sim_preds_sm.RDS")

# reproducibility
set.seed(0)

# create 2100 size 1000 datasets
out2 <- 1:2100 |> 
  purrr::map(\(x) {
    mutate(simulate_data(comb_log_clip, 
                         n = 1000, 
                         rho = rho, df = df, 
                         prop_sex = 1-mean(comb_log_clip$RIAGENDR), 
                         prop_race = table(comb_log$RIDRETH1), 
                         prop_edu = table(comb_log$DMDEDUC2)), 
           RIDRETH1 = as.numeric(RIDRETH1), 
           DMDEDUC2 = as.numeric(DMDEDUC2), 
           sim = x) 
  })

write_rds(out2, "nhanes_sim/sim_preds_lg.RDS")


# sim response data -------------------------------------------------------

source("nhanes_sim/response_fxns.R")

## smaller size -----------------------------------------------------------

out1 <- read_rds("nhanes_sim/sim_preds_sm.RDS")

set.seed(0)
out1_resp1 <- out1 |> 
  purrr::map(\(x) {
    # get dataset number 
    no <- x$sim[1] 
    x_scale <- x |> 
      mutate(across(LBX074LA:RIDAGEYR, ~c(scale(.)))) 
    df <- case_when(
      no <= 100 ~ base_case(x_scale), 
      no <= 200 ~ am1(x_scale), 
      no <= 300 ~ am2(x_scale), 
      no <= 400 ~ ap1(x_scale), 
      no <= 500 ~ ap2(x_scale), 
      no <= 600 ~ bm1(x_scale), 
      no <= 700 ~ bm2(x_scale), 
      no <= 800 ~ bp1(x_scale), 
      no <= 900 ~ bp2(x_scale), 
      no <= 1000 ~ cm1(x_scale), 
      no <= 1100 ~ cm2(x_scale), 
      no <= 1200 ~ cp1(x_scale), 
      no <= 1300 ~ cp2(x_scale), 
      no <= 1400 ~ dm1(x_scale), 
      no <= 1500 ~ dm2(x_scale), 
      no <= 1600 ~ dp1(x_scale), 
      no <= 1700 ~ dp2(x_scale), 
      .default = x #note 1701 - 2100 not set yet
    ) 
  }) |> 
  purrr::set_names(nm = c(
    rep("_base", 100), 
    rep("am1", 100), 
    rep("am2", 100), 
    rep("ap1", 100), 
    rep("ap2", 100), 
    rep("bm1", 100), 
    rep("bm2", 100), 
    rep("bp1", 100), 
    rep("bp2", 100), 
    rep("cm1", 100), 
    rep("cm2", 100), 
    rep("cp1", 100), 
    rep("cp2", 100), 
    rep("dm1", 100), 
    rep("dm2", 100), 
    rep("dp1", 100), 
    rep("dp2", 100), 
    rep("unset", 400)
  )) 

# only export first 1700
write_rds(out1_resp1[1:1700], "nhanes_sim/sim_resp_sm_a.RDS")


## larger size -----------------------------------------------------------

out2 <- read_rds("nhanes_sim/sim_preds_lg.RDS")

set.seed(0)
out2_resp1 <- out2 |> 
  purrr::map(\(x) {
    # get dataset number 
    no <- x$sim[1] 
    x_scale <- x |> 
      mutate(across(LBX074LA:RIDAGEYR, ~c(scale(.)))) 
    df <- case_when(
      no <= 100 ~ base_case(x_scale), 
      no <= 200 ~ am1(x_scale), 
      no <= 300 ~ am2(x_scale), 
      no <= 400 ~ ap1(x_scale), 
      no <= 500 ~ ap2(x_scale), 
      no <= 600 ~ bm1(x_scale), 
      no <= 700 ~ bm2(x_scale), 
      no <= 800 ~ bp1(x_scale), 
      no <= 900 ~ bp2(x_scale), 
      no <= 1000 ~ cm1(x_scale), 
      no <= 1100 ~ cm2(x_scale), 
      no <= 1200 ~ cp1(x_scale), 
      no <= 1300 ~ cp2(x_scale), 
      no <= 1400 ~ dm1(x_scale), 
      no <= 1500 ~ dm2(x_scale), 
      no <= 1600 ~ dp1(x_scale), 
      no <= 1700 ~ dp2(x_scale), 
      .default = x #note 1701 - 2100 not set yet
    ) 
  }) |> 
  purrr::set_names(nm = c(
    rep("_base", 100), 
    rep("am1", 100), 
    rep("am2", 100), 
    rep("ap1", 100), 
    rep("ap2", 100), 
    rep("bm1", 100), 
    rep("bm2", 100), 
    rep("bp1", 100), 
    rep("bp2", 100), 
    rep("cm1", 100), 
    rep("cm2", 100), 
    rep("cp1", 100), 
    rep("cp2", 100), 
    rep("dm1", 100), 
    rep("dm2", 100), 
    rep("dp1", 100), 
    rep("dp2", 100), 
    rep("unset", 400)
  )) 

# only export first 1700
write_rds(out2_resp1[1:1700], "nhanes_sim/sim_resp_lg_a.RDS")

## race/ethnicity ---------------------------------------------------------

# read output back in, size 252
out1 <- read_rds("nhanes_sim/sim_preds_sm.RDS")

set.seed(0)
out1_resp1_re <- out1 |> 
  purrr::map(\(x) {
    # get dataset number 
    no <- x$sim[1] 
    x <- x |> 
      mutate(across(LBX074LA:RIDAGEYR, ~c(scale(.)))) 
    df <- case_when(
      no <= 1700 ~ x, #note 1 - 1700 are chemxchem
      no <= 1800 ~ em1(x), 
      no <= 1900 ~ em2(x), 
      no <= 2000 ~ ep1(x), 
      no <= 2100 ~ ep2(x), 
      .default = x 
    ) 
  }) |> 
  purrr::set_names(nm = c(
    rep("unset", 1700), 
    rep("em1", 100), 
    rep("em2", 100), 
    rep("ep1", 100), 
    rep("ep2", 100)
  )) 
  
write_rds(out1_resp1_re[1701:2100], "nhanes_sim/sim_resp_sm_re.RDS")

# read output back in, size 1000
out2 <- read_rds("nhanes_sim/sim_preds_lg.RDS")

set.seed(0)
out2_resp1_re <- out2 |> 
  purrr::map(\(x) {
    # get dataset number 
    no <- x$sim[1] 
    x <- x |> 
      mutate(across(LBX074LA:RIDAGEYR, ~c(scale(.)))) 
    df <- case_when(
      no <= 1700 ~ x, #note 1 - 1700 are chemxchem
      no <= 1800 ~ em1(x), 
      no <= 1900 ~ em2(x), 
      no <= 2000 ~ ep1(x), 
      no <= 2100 ~ ep2(x), 
      .default = x 
    ) 
  }) |> 
  purrr::set_names(nm = c(
    rep("unset", 1700), 
    rep("em1", 100), 
    rep("em2", 100), 
    rep("ep1", 100), 
    rep("ep2", 100)
  )) 

# get output
out2_resp1_re <- out2_resp1_re[1701:2100]
write_rds(out2_resp1_re, "nhanes_sim/sim_resp_lg_re.RDS")

# fit oracle/mlr ----------------------------------------------------------

# only use the first 1700 datasets for now
out2_resp1 <- read_rds("nhanes_sim/sim_resp_lg_a.RDS")

mlrs <- vector(mode='list', length = 1700)
names(mlrs) <- names(out1_resp1)
mlrtimes <- vector(mode = 'list', length = 1700)
names(mlrtimes) <- names(out1_resp1)

chems_only <- vector(mode = 'list', length = 1700)
names(chems_only) <- names(out1_resp1)

oracles <- vector(mode='list', length = 1700)
names(oracles) <- names(out1_resp1)
oracletimes <- vector(mode = 'list', length = 1700)
names(oracletimes) <- names(out1_resp1)


for(i in 1:1700) {
  df <- out2_resp1[[i]] |> 
    mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
    select(-sim)
  
  start.time <- Sys.time()
  mlrs[[i]] <- lm(y ~ ., data = df)
  end.time <- Sys.time()
  mlrtimes[[i]] <- end.time - start.time
  
  # chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
  #                         I(1/(1+exp(-4*LBX194LA))) + 
  #                         I(1/(1+exp(-4*LBXPCBLA))) + 
  #                         I(LBXF04LA^2) + LBXF04LA ,
  #                       data = df)
  
  if(i <= 100) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA,
                       data = df)
  } else if (i <= 300) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXD05LA*LBX194LA + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXD05LA*LBX194LA, 
                       data = df)
  } else if (i <= 500) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXD05LA*((LBX194LA-1)^2)) + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXD05LA*((LBX194LA-1)^2)), 
                       data = df)
  } else if (i <= 700) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXF08LA*LBXF03LA + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXF08LA*LBXF03LA, 
                       data = df)
  } else if (i <= 900) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXF08LA*((LBXF03LA-1)^2)) + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXF08LA*((LBXF03LA-1)^2)), 
                       data = df)
  } else if (i <= 1100) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBX074LA*LBX194LA + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBX074LA*LBX194LA, 
                       data = df)
  } else if (i <= 1300) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBX074LA*((LBX194LA-1)^2)) + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBX074LA*((LBX194LA-1)^2)), 
                       data = df)
  } else if (i <= 1500) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXD05LA:LBXPCBLA:LBX194LA + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         LBXD05LA:LBXPCBLA:LBX194LA, 
                       data = df)
  } else if (i <= 1700) {
    start.time <- Sys.time()
    oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA) + 
                         BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                         LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                         RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                       data = df)
    end.time <- Sys.time()
    oracletimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                         I(1/(1+exp(-4*LBX194LA))) + 
                         I(1/(1+exp(-4*LBXPCBLA))) + 
                         I(LBXF04LA^2) + LBXF04LA +
                         I(LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA), 
                       data = df)
  }
}

# write out 
write_rds(mlrs, "nhanes_sim/mods/mlr_mods_lg.RDS")
write_rds(mlrtimes, "nhanes_sim/times/mlr_mods_lg_times.RDS")
write_rds(oracles, "nhanes_sim/mods/oracle_mods_lg.RDS")
write_rds(oracletimes, "nhanes_sim/times/oracle_mods_lg_times.RDS")
write_rds(chems_only, "nhanes_sim/mods/chem_mods_lg.RDS")

## check p-values ---------------------------------------------------------


### chem --------------------------------------------------------------------

#chem model output small
chem_mods <- read_rds("nhanes_sim/mods/chem_mods_lg.RDS")
rsq_chem <- chem_mods |> 
  purrr::map_df(\(x) {
    data.frame(
      rsq = summary(x)$r.squared
    )
  }) |> 
  mutate(name = names(chem_mods))

rsqsmplot <- rsq_chem |> 
  ggplot(aes(x = rsq)) +
  geom_density() + 
  facet_wrap(~name, 
             # labeller = as_labeller(appender1, 
             #                        default = label_parsed), 
             ncol = 4) +
  labs(y = "Density", x = latex2exp::TeX("R$^2$"))
rsqsmplot
# ggsave("sim/figs/chem_sm_rsq.png", width = 7, height = 5)


### oracle -----------------------------------------------------------------



# readin
oracle_mods <- read_rds("nhanes_sim/mods/oracle_mods_lg.RDS")

# extract r squared
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
             # labeller = as_labeller(appender1,
             #                        default = label_parsed),
             ncol = 4) +
  labs(y = "Density", x = TeX("R$^2$"))
# ggsave("sim/figs/oracle_sm_rsq.png", width = 7.5, height = 5)

# extract p values
keepnames <- c('(Intercept)', 'LBXD05LA', 'LBX074LA', 
               'I(1/(1 + exp(-4 * LBX194LA)))', 'I(1/(1 + exp(-4 * LBXPCBLA)))', 
               'I(LBXF04LA^2)', 
               'LBXF04LA', 'LBX194LA', 'LBXF08LA', 'LBXF03LA',
               'BMXBMI', 'LBXCOT', 'LBXWBCSI', 'LBXLYPCT', 
               'LBXMOPCT', 'LBXNEPCT', 'LBXEOPCT', 'LBXBAPCT', 'RIDAGEYR', 
               'RIAGENDR1', 'RIDRETH12', 'RIDRETH13', 'RIDRETH14', 'RIDRETH15', 
               'DMDEDUC22', 'DMDEDUC23', 'DMDEDUC24', 'DMDEDUC25')
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
             # labeller = as_labeller(appender,
             #                        default = label_parsed)
             ) +
  labs(y = "Density", x = "p-value")
# ggsave("sim/figs/oracle_sm_pval_dist.png", width = 7.5, height = 5)

#extract power at alpha = 0.05
oracle_sm_power <- pval |>
  group_by(name) |>
  summarize(power = sum(pval < 0.05)/n())

oracle_sm_power

# write_csv(oracle_sm_power, "sim/tables/oracle_sm_power.csv")


