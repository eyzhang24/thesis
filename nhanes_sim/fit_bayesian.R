library(tidyverse)
library(bkmr)

#########
# run bkmr
#########

df <- read_rds("nhanes_data/processed_data.RDS")

# fit bkmr on the original dataset

# prep data
Z <- df |> select(LBX074LA:LBXF08LA)
X <- df |> 
  bind_cols(
    data.frame(model.matrix(
      ~ RIDRETH1 - 1, data = mutate(df, RIDRETH1 = as.factor(RIDRETH1)))), 
    data.frame(model.matrix(
      ~ DMDEDUC2 - 1, data = mutate(df, DMDEDUC2 = as.factor(DMDEDUC2))))
  ) |> 
  select(RIDAGEYR, RIAGENDR, RIDRETH12:RIDRETH15, DMDEDUC22:DMDEDUC25, 
         BMXBMI:LBXBAPCT)
y <- df$y <- df$TELOMEAN

# fit model
set.seed(0)
start.time <- Sys.time()
mod <- kmbayes(y = y, Z = Z, X = X, 
               iter = 500, varsel = TRUE)
end.time <- Sys.time() 


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