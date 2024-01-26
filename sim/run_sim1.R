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



