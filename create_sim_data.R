library(tidyverse)
library(copula)

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
  #gender to 0's and 1's
  mutate(RIAGENDR = RIAGENDR - 1)

# spearman rho


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
