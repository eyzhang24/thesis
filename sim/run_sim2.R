library(tidyverse)
library(rslurm)
library(bkmr)
library(NLinteraction)
# library(copula)

# read in data
comb_small <- read_csv("madres_data/base_data.csv")

# log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log))
comb_log_clip <- comb_log |> 
  select(7:19)

########
# simulate predictor data
########
cfit_gaus <- read_rds("sim/gauscop.RDS")

#get rho 
rho <- coef(cfit_gaus)

#create function for simulation
simulate_data <- function(data, n, rho, prop_smoke, prop_race) {
  #'data = original observed data
  #'n = sample size
  #'rho = rho values from normal copula
  #'prop_smoke = proportion smoke from observed dataset
  #'prop_race = table with race/eth values
  
  #simulate pseudo-observations from copula
  samp <- rCopula(n, normalCopula(rho, dim = ncol(data), dispstr = "un"))
  #transform pseudo-observations to observed marginal distributions
  sampt <- 1:ncol(data) |> 
    purrr::map_dfc(
      \(x) {
        if(names(data)[x] == "smoke") {
          df <- data.frame(ifelse(samp[,x] < prop_smoke, 0, 1), 
                           row.names = NULL)
        } else {
          df <- data.frame(quantile(data[[x]], probs = samp[,x]), 
                           row.names = NULL)
        }
        names(df) <- names(data)[x]
        return(df)
      }
    ) |> 
    mutate(race = sample(x = names(prop_race), prob = prop_race,
                         size = n, replace = T)) |> 
    relocate(race)
  return(sampt)
}

# create function to run size 252 samples on hpc
run_sim1 <- function() {
  set.seed(0)
  out <- 1:2100 |> 
    purrr::map(\(x) {
      mutate(simulate_data(comb_log_clip, n = nrow(comb_log_clip), rho = rho, 
                           prop_smoke = 1-mean(comb_log_clip$smoke), 
                           prop_race = table(comb_log$race)), 
             race = as.numeric(race), 
             sim = x) 
    })
  return(out)
}

# send job to hpc for size 252 samples
sjob1 <- slurm_call(run_sim1, 
                    global_objects = c('comb_log', 'comb_log_clip', 
                                       'rho', 'simulate_data'),
                    jobname = 'sim_data1')

# check on status of job
get_job_status(sjob1)$completed
system('squeue --me')

# get output
out1 <- get_slurm_out(sjob1)
write_rds(out1, "sim/sim_preds_sm.RDS")

# create function to run size 1000 samples on hpc
run_sim2 <- function() {
  set.seed(1)
  out <- 1:2100 |> 
    purrr::map(\(x) {
      mutate(simulate_data(comb_log_clip, n = 1000, rho = rho, 
                           prop_smoke = 1-mean(comb_log_clip$smoke), 
                           prop_race = table(comb_log$race)), 
             race = as.numeric(race), 
             sim = x)
    })
  return(out)
}

# send job to hpc for size 1000 samples
sjob2 <- slurm_call(run_sim2, 
                    global_objects = c('comb_log', 'comb_log_clip', 
                                       'rho', 'simulate_data'),
                    jobname = 'sim_data2')

# check on status of job
get_job_status(sjob2)$completed
system('squeue --me')

# get output
out2 <- get_slurm_out(sjob2)
write_rds(out2, "sim/sim_preds_lg.RDS")

#########
# simulate response, interxns b/t chemicals
#########
source("sim/response_fxns.R")

# read output back in, size 252
out1 <- read_rds("sim/sim_preds_sm.RDS")

run_resp1 <- function() {
  set.seed(0)
  out1_resp1 <- out1 |> 
    purrr::map(\(x) {
      # get dataset number 
      no <- x$sim[1] 
      x <- x |> 
        mutate(across(As:Sn, ~c(scale(.)))) 
      df <- case_when(
        no <= 100 ~ base_case(x), 
        no <= 200 ~ am1(x), 
        no <= 300 ~ am2(x), 
        no <= 400 ~ ap1(x), 
        no <= 500 ~ ap2(x), 
        no <= 600 ~ bm1(x), 
        no <= 700 ~ bm2(x), 
        no <= 800 ~ bp1(x), 
        no <= 900 ~ bp2(x), 
        no <= 1000 ~ cm1(x), 
        no <= 1100 ~ cm2(x), 
        no <= 1200 ~ cp1(x), 
        no <= 1300 ~ cp2(x), 
        no <= 1400 ~ dm1(x), 
        no <= 1500 ~ dm2(x), 
        no <= 1600 ~ dp1(x), 
        no <= 1700 ~ dp2(x), 
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
  return(out1_resp1)
}

runrespsm <- slurm_call(
  run_resp1, 
  global_objects = c('out1', 'base_case', 
                     'am1', 'am2', 'ap1', 'ap2', 
                     'bm1', 'bm2', 'bp1', 'bp2', 
                     'cm1', 'cm2', 'cp1', 'cp2', 
                     'dm1', 'dm2', 'dp1', 'dp2'),
  jobname = 'sim_resp1')

# check on status of job
get_job_status(runrespsm)$completed
system('squeue --me')

# get output
out1_resp1 <- get_slurm_out(runrespsm)
out1_resp1 <- read_rds("_rslurm_sim_resp1/results_0.RDS")
out1_resp1 <- out1_resp1[1:1700]
write_rds(out1_resp1, "sim/sim_resp_sm_a.RDS")

# read output back in, size 1000
out2 <- read_rds("sim/sim_preds_lg.RDS")

run_resp2 <- function() {
  set.seed(0)
  out2_resp1 <- out2 |> 
    purrr::map(\(x) {
      # get dataset number 
      no <- x$sim[1]
      x <- x |> 
        mutate(across(As:Sn, ~c(scale(.)))) 
      df <- case_when(
        no <= 100 ~ base_case(x), 
        no <= 200 ~ am1(x), 
        no <= 300 ~ am2(x), 
        no <= 400 ~ ap1(x), 
        no <= 500 ~ ap2(x), 
        no <= 600 ~ bm1(x), 
        no <= 700 ~ bm2(x), 
        no <= 800 ~ bp1(x), 
        no <= 900 ~ bp2(x), 
        no <= 1000 ~ cm1(x), 
        no <= 1100 ~ cm2(x), 
        no <= 1200 ~ cp1(x), 
        no <= 1300 ~ cp2(x), 
        no <= 1400 ~ dm1(x), 
        no <= 1500 ~ dm2(x), 
        no <= 1600 ~ dp1(x), 
        no <= 1700 ~ dp2(x), 
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
  return(out2_resp1)
}

sjob4 <- slurm_call(
  run_resp2, 
  global_objects = c('out2', 'base_case', 
                     'am1', 'am2', 'ap1', 'ap2', 
                     'bm1', 'bm2', 'bp1', 'bp2', 
                     'cm1', 'cm2', 'cp1', 'cp2', 
                     'dm1', 'dm2', 'dp1', 'dp2'),
  jobname = 'sim_resp2')

# check on status of job
get_job_status(sjob4)$completed
system('squeue --me')

# get output
out2_resp1 <- get_slurm_out(sjob4)
# out2_resp1 <- read_rds("_rslurm_sim_resp2/results_0.RDS")
out2_resp1 <- out2_resp1[1:1700]
write_rds(out2_resp1, "sim/sim_resp_lg_a.RDS")

# for(i in seq(50, 2050 ,50)) {
#   print(names(out2_resp1)[i])
#   print(summary(out2_resp1[[i]]$y))
# }

#########
# run mlr & oracle
#########

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")
run_mlr_sm <- function() {
  # initialize vectors
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
    df <- out1_resp1[[i]] |> 
      mutate(race = as.factor(race), smoke = as.factor(smoke)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrs[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ Hg + Sb + 
                            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))), 
                          data = df)
    
    if(i <= 100) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 300) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg*Ni + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 500) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Ni-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 700) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Cd*As + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 900) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Cd*((As-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 1100) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg*Co + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 1300) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Co-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 1500) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg:Ni:Tl + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    } else if (i <= 1700) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Ni-1)^2)*Tl) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
    }
  }
  return(list(mlrs, mlrtimes, oracles, oracletimes, chems_only))
}

sjob5 <- slurm_call(
  run_mlr_sm, 
  global_objects = c('out1_resp1'),
  jobname = 'mlr_sm')

# check on status of job
get_job_status(sjob5)$completed
system('squeue --me')

# get output
mlr_sm <- read_rds("_rslurm_mlr_sm/results_0.RDS")
# mlr_sm <- get_slurm_out(sjob5)
mlr_mods <- mlr_sm[[1]]
mlr_times <- mlr_sm[[2]]
oracle_mods <- mlr_sm[[3]]
oracle_times <- mlr_sm[[4]]
chem_mods <- mlr_sm[[5]]
write_rds(mlr_mods, "sim/mlr_mods_sm.RDS")
write_rds(mlr_times, "sim/mlr_mods_sm_times.RDS")
write_rds(oracle_mods, "sim/oracle_mods_sm.RDS")
write_rds(oracle_times, "sim/oracle_mods_sm_times.RDS")
write_rds(chem_mods, "sim/chem_mods_sm.RDS")

### larger sample size

out2_resp1 <- read_rds("sim/sim_resp_lg_a.RDS")
run_mlr_lg <- function() {
  # initialize vectors
  mlrl <- vector(mode='list', length = 1700)
  names(mlrl) <- names(out2_resp1)
  mlrtimel <- vector(mode = 'list', length = 1700)
  names(mlrtimel) <- names(out2_resp1)
  
  chems_onlyl <- vector(mode = 'list', length = 1700)
  names(chems_onlyl) <- names(out2_resp1)
  
  oraclel <- vector(mode='list', length = 1700)
  names(oraclel) <- names(out2_resp1)
  oracletimel <- vector(mode = 'list', length = 1700)
  names(oracletimel) <- names(out2_resp1)
  
  for(i in 1:1700) {
    df <- out2_resp1[[i]] |> 
      mutate(race = as.factor(race), smoke = as.factor(smoke)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrl[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimel[[i]] <- end.time - start.time
    
    chems_onlyl[[i]] <- lm(y ~ Hg + Sb + 
                             I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))), 
                           data = df)
    
    if(i <= 100) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 300) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg*Ni + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 500) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Ni-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 700) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Cd*As + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 900) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Cd*((As-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 1100) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg*Co + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 1300) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Co-1)^2)) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 1500) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           Hg:Ni:Tl + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    } else if (i <= 1700) {
      start.time <- Sys.time()
      oraclel[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*((Ni-1)^2)*Tl) + 
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimel[[i]] <- end.time - start.time
    }
  }
  return(list(mlrl, mlrtimel, oraclel, oracletimel, chems_onlyl))
}

sjob6 <- slurm_call(
  run_mlr_lg, 
  global_objects = c('out2_resp1'),
  jobname = 'mlr_lg')

# check on status of job
get_job_status(sjob6)$completed
system('squeue --me')

# get output
mlr_lg <- get_slurm_out(sjob6)
mlr_modl <- mlr_lg[[1]]
mlr_timel <- mlr_lg[[2]]
oracle_modl <- mlr_lg[[3]]
oracle_timel <- mlr_lg[[4]]
chem_modl <- mlr_lg[[5]]
write_rds(mlr_modl, "sim/mlr_mods_lg.RDS")
write_rds(mlr_timel, "sim/mlr_mods_lg_times.RDS")
write_rds(oracle_modl, "sim/oracle_mods_lg.RDS")
write_rds(oracle_timel, "sim/oracle_mods_lg_times.RDS")
write_rds(chem_modl, "sim/chem_mods_lg.RDS")

#########
# run bkmr
#########

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")

run_bkmr_sm <- function(vector) {
  bkmr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  
  set.seed(0)
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df <- out1_resp1[[i]]
    Z <- df |> 
      select(As:Sn)
    X <- df |>
      bind_cols(
        data.frame(model.matrix(~ race-1, data = 
                                  mutate(df, race = as.factor(race))))
      ) |> 
      select(race2:race5, smoke:bmi)
    y <- df$y
    
    # fit model and save time
    start.time <- Sys.time()
    mod <- kmbayes(y = y, Z = Z, X = X, 
                   iter = 50000, verbose = FALSE, varsel = TRUE)
    end.time <- Sys.time()
    bkmr_times[[i]] <- end.time - start.time
    
    # save model and remove from memory
    write_rds(mod, file = paste0("mods/bkmr_sm_", names(out1_resp1)[i], "_", i, ".RDS"))
    rm(mod)
  }
  write_rds(bkmr_times, file = "bkmr_sm_times.RDS")
  return(bkmr_times)
}

bksmjob01 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1:100),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm01',
  slurm_options = list(mem = '8G'))

bksmjob02 <- slurm_call(
  run_bkmr_sm, params = list(vector = 101:200),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm02',
  slurm_options = list(mem = '8G'))

bksmjob03 <- slurm_call(
  run_bkmr_sm, params = list(vector = 201:300),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm03',
  slurm_options = list(mem = '8G'))

bksmjob04 <- slurm_call(
  run_bkmr_sm, params = list(vector = 301:400),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm04',
  slurm_options = list(mem = '8G'))

bksmjob05 <- slurm_call(
  run_bkmr_sm, params = list(vector = 401:500),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm05',
  slurm_options = list(mem = '8G'))

bksmjob06 <- slurm_call(
  run_bkmr_sm, params = list(vector = 501:600),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm06',
  slurm_options = list(mem = '8G'))

bksmjob07 <- slurm_call(
  run_bkmr_sm, params = list(vector = 601:700),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm07',
  slurm_options = list(mem = '8G'))

bksmjob08 <- slurm_call(
  run_bkmr_sm, params = list(vector = 701:800),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm08',
  slurm_options = list(mem = '8G'))

bksmjob09 <- slurm_call(
  run_bkmr_sm, params = list(vector = 801:900),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm09',
  slurm_options = list(mem = '8G'))

bksmjob10 <- slurm_call(
  run_bkmr_sm, params = list(vector = 901:1000),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm10',
  slurm_options = list(mem = '8G'))

bksmjob11 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1001:1100),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm11',
  slurm_options = list(mem = '8G'))

bksmjob12 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1101:1200),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm12',
  slurm_options = list(mem = '8G'))

bksmjob13 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1201:1300),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm13',
  slurm_options = list(mem = '8G'))

bksmjob14 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1301:1400),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm14',
  slurm_options = list(mem = '8G'))

bksmjob15 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1401:1500),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm15',
  slurm_options = list(mem = '8G'))

bksmjob16 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1501:1600),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm16',
  slurm_options = list(mem = '8G'))

bksmjob17 <- slurm_call(
  run_bkmr_sm, params = list(vector = 1601:1700),
  global_objects = c('out1_resp1'),
  jobname = 'bkmr_sm17',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')

### larger sample size
out2_resp1 <- read_rds("sim/sim_resp_lg_a.RDS")

run_bkmr_lg <- function(vector) {
  bkmr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  
  set.seed(0)
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df <- out2_resp1[[i]]
    Z <- df |> 
      select(As:Sn)
    X <- df |>
      bind_cols(
        data.frame(model.matrix(~ race-1, data = 
                                  mutate(df, race = as.factor(race))))
      ) |> 
      select(race2:race5, smoke:bmi)
    y <- df$y
    knots <- fields::cover.design(Z, nd = 100)$design
    
    # fit model and save time
    start.time <- Sys.time()
    mod <- kmbayes(y = y, Z = Z, X = X, knots = knots,
                   iter = 50000, verbose = FALSE, varsel = TRUE)
    end.time <- Sys.time()
    bkmr_times[[i]] <- end.time - start.time
    
    # save model and remove from memory
    write_rds(mod, file = paste0("mods/bkmr_sm_", names(out2_resp1)[i], "_", i, ".RDS"))
    rm(mod)
  }
  write_rds(bkmr_times, file = "bkmr_lg_times.RDS")
  # return(bkmr_times)
}

# run for all samples
bklgjob01 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1:100),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg01',
  slurm_options = list(mem = '8G'))

bklgjob02 <- slurm_call(
  run_bkmr_lg, params = list(vector = 101:200),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg02',
  slurm_options = list(mem = '8G'))

bklgjob03 <- slurm_call(
  run_bkmr_lg, params = list(vector = 201:300),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg03',
  slurm_options = list(mem = '8G'))

bklgjob04 <- slurm_call(
  run_bkmr_lg, params = list(vector = 301:400),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg04',
  slurm_options = list(mem = '8G'))

bklgjob05 <- slurm_call(
  run_bkmr_lg, params = list(vector = 401:500),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg05',
  slurm_options = list(mem = '8G'))

bklgjob06 <- slurm_call(
  run_bkmr_lg, params = list(vector = 501:600),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg06',
  slurm_options = list(mem = '8G'))

bklgjob07 <- slurm_call(
  run_bkmr_lg, params = list(vector = 601:700),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg07',
  slurm_options = list(mem = '8G'))

bklgjob08 <- slurm_call(
  run_bkmr_lg, params = list(vector = 701:800),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg08',
  slurm_options = list(mem = '8G'))

bklgjob09 <- slurm_call(
  run_bkmr_lg, params = list(vector = 801:900),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg09',
  slurm_options = list(mem = '8G'))

bklgjob10 <- slurm_call(
  run_bkmr_lg, params = list(vector = 901:1000),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg10',
  slurm_options = list(mem = '8G'))

bklgjob11 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1001:1100),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg11',
  slurm_options = list(mem = '8G'))

bklgjob12 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1101:1200),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg12',
  slurm_options = list(mem = '8G'))

bklgjob13 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1201:1300),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg13',
  slurm_options = list(mem = '8G'))

bklgjob14 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1301:1400),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg14',
  slurm_options = list(mem = '8G'))

bklgjob15 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1401:1500),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg15',
  slurm_options = list(mem = '8G'))

bklgjob16 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1501:1600),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg16',
  slurm_options = list(mem = '8G'))

bklgjob17 <- slurm_call(
  run_bkmr_lg, params = list(vector = 1601:1700),
  global_objects = c('out2_resp1'),
  jobname = 'bkmr_lg17',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')

#########
# run bsr
#########

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")

run_bsr_sm <- function(vector) {
  bsr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  
  set.seed(0)
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df <- out1_resp1[[i]]
    X <- df |> 
      select(As:Sn) |> 
      as.matrix.data.frame()
    C <- df |>
      bind_cols(
        data.frame(model.matrix(~ race-1, data = 
                                  mutate(df, race = as.factor(race))))
      ) |> 
      select(race2:race5, smoke:bmi) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # fit model for d = {1, 2, 3, 4, 5} and save time
    list_times <- vector(mode = "list", length = 2)
    
    start.time <- Sys.time()
    mod1 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 1)
    mod2 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 2)
    mod3 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 3)
    mod4 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 4)
    end.time <- Sys.time()
    list_times[[1]] <- end.time - start.time
    
    ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
    
    print(paste0("-------- chose ", ind, " -----------"))
    
    start.time <- Sys.time()
    mod <- NLint(Y = Y, X = X, C = C, 
                 nIter = 50000, nBurn = 25000, ns = ind)
    end.time <- Sys.time()
    list_times[[2]] <- end.time - start.time
    
    bsr_times[[i]] <- list_times
    
    # save model and remove from memory
    write_rds(mod, file = 
                paste0("mods/bsr_sm_", names(out1_resp1)[i], "_", i, 
                       "df", ind, ".RDS"))
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  write_rds(bsr_times, file = "bsr_sm_times.RDS")
  return(bsr_times)
}

rjob01 <- slurm_call(
  run_bsr_sm, params = list(vector = 1:100),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm01',
  slurm_options = list(mem = '8G'))

rjob02 <- slurm_call(
  run_bsr_sm, params = list(vector = 101:200),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm02',
  slurm_options = list(mem = '8G'))

rjob03 <- slurm_call(
  run_bsr_sm, params = list(vector = 201:300),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm03',
  slurm_options = list(mem = '8G'))

rjob04 <- slurm_call(
  run_bsr_sm, params = list(vector = 301:400),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm04',
  slurm_options = list(mem = '8G'))

rjob05 <- slurm_call(
  run_bsr_sm, params = list(vector = 401:500),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm05',
  slurm_options = list(mem = '8G'))

rjob06 <- slurm_call(
  run_bsr_sm, params = list(vector = 501:600),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm06',
  slurm_options = list(mem = '8G'))

rjob07 <- slurm_call(
  run_bsr_sm, params = list(vector = 601:700),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm07',
  slurm_options = list(mem = '8G'))

rjob08 <- slurm_call(
  run_bsr_sm, params = list(vector = 701:800),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm08',
  slurm_options = list(mem = '8G'))

rjob09 <- slurm_call(
  run_bsr_sm, params = list(vector = 801:900),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm09',
  slurm_options = list(mem = '8G'))

rjob10 <- slurm_call(
  run_bsr_sm, params = list(vector = 901:1000),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm10',
  slurm_options = list(mem = '8G'))

rjob11 <- slurm_call(
  run_bsr_sm, params = list(vector = 1001:1100),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm11',
  slurm_options = list(mem = '8G'))

rjob12 <- slurm_call(
  run_bsr_sm, params = list(vector = 1101:1200),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm12',
  slurm_options = list(mem = '8G'))

rjob13 <- slurm_call(
  run_bsr_sm, params = list(vector = 1201:1300),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm13',
  slurm_options = list(mem = '8G'))

rjob14 <- slurm_call(
  run_bsr_sm, params = list(vector = 1301:1400),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm14',
  slurm_options = list(mem = '8G'))

rjob15 <- slurm_call(
  run_bsr_sm, params = list(vector = 1401:1500),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm15',
  slurm_options = list(mem = '8G'))

rjob16 <- slurm_call(
  run_bsr_sm, params = list(vector = 1501:1600),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm16',
  slurm_options = list(mem = '8G'))

rjob17 <- slurm_call(
  run_bsr_sm, params = list(vector = 1601:1700),
  global_objects = c('out1_resp1'),
  jobname = 'bsr_sm17',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')


### larger sample size

out2_resp1 <- read_rds("sim/sim_resp_lg_a.RDS")

run_bsr_lg <- function(vector) {
  bsr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  
  set.seed(0)
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df <- out2_resp1[[i]]
    X <- df |> 
      select(As:Sn) |> 
      as.matrix.data.frame()
    C <- df |>
      bind_cols(
        data.frame(model.matrix(~ race-1, data = 
                                  mutate(df, race = as.factor(race))))
      ) |> 
      select(race2:race5, smoke:bmi) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # fit model for d = {1, 2, 3, 4, 5} and save time
    list_times <- vector(mode = "list", length = 2)
    
    start.time <- Sys.time()
    mod1 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 1)
    mod2 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 2)
    mod3 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 3)
    mod4 <- NLint(Y = Y, X = X, C = C, 
                  nIter = 5000, nBurn = 2500, ns = 4)
    end.time <- Sys.time()
    list_times[[1]] <- end.time - start.time
    
    ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
    
    print(paste0("-------- chose ", ind, " -----------"))
    
    start.time <- Sys.time()
    mod <- NLint(Y = Y, X = X, C = C, 
                 nIter = 50000, nBurn = 25000, ns = ind)
    end.time <- Sys.time()
    list_times[[2]] <- end.time - start.time
    
    bsr_times[[i]] <- list_times
    
    # save model and remove from memory
    write_rds(mod, file = 
                paste0("mods/bsr_sm_", names(out2_resp1)[i], "_", i, 
                       "df", ind, ".RDS"))
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  write_rds(bsr_times, file = "bsr_sm_times.RDS")
  # return(bsr_times)
}

qjob01 <- slurm_call(
  run_bsr_lg, params = list(vector = 1:100),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg01',
  slurm_options = list(mem = '8G'))

qjob02 <- slurm_call(
  run_bsr_lg, params = list(vector = 101:200),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg02',
  slurm_options = list(mem = '8G'))

qjob03 <- slurm_call(
  run_bsr_lg, params = list(vector = 201:300),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg03',
  slurm_options = list(mem = '8G'))

qjob04 <- slurm_call(
  run_bsr_lg, params = list(vector = 301:400),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg04',
  slurm_options = list(mem = '8G'))

qjob05 <- slurm_call(
  run_bsr_lg, params = list(vector = 401:500),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg05',
  slurm_options = list(mem = '8G'))

qjob06 <- slurm_call(
  run_bsr_lg, params = list(vector = 501:600),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg06',
  slurm_options = list(mem = '8G'))

qjob07 <- slurm_call(
  run_bsr_lg, params = list(vector = 601:700),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg07',
  slurm_options = list(mem = '8G'))

qjob08 <- slurm_call(
  run_bsr_lg, params = list(vector = 701:800),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg08',
  slurm_options = list(mem = '8G'))

qjob09 <- slurm_call(
  run_bsr_lg, params = list(vector = 801:900),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg09',
  slurm_options = list(mem = '8G'))

qjob10 <- slurm_call(
  run_bsr_lg, params = list(vector = 901:1000),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg10',
  slurm_options = list(mem = '8G'))

qjob11 <- slurm_call(
  run_bsr_lg, params = list(vector = 1001:1100),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg11',
  slurm_options = list(mem = '8G'))

qjob12 <- slurm_call(
  run_bsr_lg, params = list(vector = 1101:1200),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg12',
  slurm_options = list(mem = '8G'))

qjob13 <- slurm_call(
  run_bsr_lg, params = list(vector = 1201:1300),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg13',
  slurm_options = list(mem = '8G'))

qjob14 <- slurm_call(
  run_bsr_lg, params = list(vector = 1301:1400),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg14',
  slurm_options = list(mem = '8G'))

qjob15 <- slurm_call(
  run_bsr_lg, params = list(vector = 1401:1500),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg15',
  slurm_options = list(mem = '8G'))

qjob16 <- slurm_call(
  run_bsr_lg, params = list(vector = 1501:1600),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg16',
  slurm_options = list(mem = '8G'))

qjob17 <- slurm_call(
  run_bsr_lg, params = list(vector = 1601:1700),
  global_objects = c('out2_resp1'),
  jobname = 'bsr_lg17',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')

