library(tidyverse)
library(rslurm)
library(bkmr)
library(NLinteraction)

# functions for creating response
# interxn in smaller group, smaller effect size
em1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5 + 0.5*Hg, # 1.5x in group 2
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5) +
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

# interxn in smaller group, larger effect size
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

# interxn in larger group, smaller effect size
ep1 <- function(df) {
  mutate(df, y = 
           Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
           age + 0.5*bmi + 
           case_when(race == 1 ~ 1, 
                     race == 2 ~ 1.5, 
                     race == 3 ~ 1, 
                     race == 4 ~ 1, 
                     race == 5 ~ 1.5 + 0.5*Hg) + # 1.5x in group 5
           ifelse(smoke == 1, -1, 0.5) +
           rnorm(nrow(df), 0, 5))
}

# interxn in larger group, larger effect size
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

# read output back in, size 252
out1 <- read_rds("sim/sim_preds_sm.RDS")

run_resp1_re <- function() {
  set.seed(0)
  out1_resp1 <- out1 |> 
    purrr::map(\(x) {
      # get dataset number 
      no <- x$sim[1] 
      x <- x |> 
        mutate(across(age:Sn, ~c(scale(.)))) 
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
  return(out1_resp1)
}

runrespsm_re <- slurm_call(
  run_resp1_re, 
  global_objects = c('out1', 
                     'em1', 'em2', 'ep1', 'ep2'),
  jobname = 'sim_resp1_re')

# check on status of job
get_job_status(runrespsm_re)$completed
system('squeue --me')

# get output
out1_resp1_re <- get_slurm_out(runrespsm_re)
# out1_resp1_re <- read_rds("_rslurm_sim_resp1_re/results_0.RDS")
out1_resp1_re <- out1_resp1_re[1701:2100]
write_rds(out1_resp1_re, "sim/sim_resp_sm_re.RDS")

# read output back in, size 1000
out2 <- read_rds("sim/sim_preds_lg.RDS")

run_resp2_re <- function() {
  set.seed(0)
  out2_resp1 <- out2 |> 
    purrr::map(\(x) {
      # get dataset number 
      no <- x$sim[1] 
      x <- x |> 
        mutate(across(age:Sn, ~c(scale(.)))) 
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
  return(out2_resp1)
}

runresplg_re <- slurm_call(
  run_resp2_re, 
  global_objects = c('out2', 
                     'em1', 'em2', 'ep1', 'ep2'),
  jobname = 'sim_resp2_re')

# check on status of job
get_job_status(runresplg_re)$completed
system('squeue --me')

# get output
out2_resp1_re <- get_slurm_out(runresplg_re)
# out2_resp1 <- read_rds("_rslurm_sim_resp2/results_0.RDS")
out2_resp1_re <- out2_resp1_re[1701:2100]
write_rds(out2_resp1_re, "sim/sim_resp_lg_re.RDS")

#########
# run mlr & oracle
#########

### smaller sample size

out1_resp1_re <- read_rds("sim/sim_resp_sm_re.RDS")
run_mlr_sm_re <- function() {
  # initialize vectors
  mlrs <- vector(mode='list', length = 400)
  names(mlrs) <- names(out1_resp1_re)
  mlrtimes <- vector(mode = 'list', length = 400)
  names(mlrtimes) <- names(out1_resp1_re)
  
  chems_only <- vector(mode = 'list', length = 400)
  names(chems_only) <- names(out1_resp1_re)
  chems_oracle <- vector(mode = 'list', length = 400) 
  names(mlrs) <- names(out1_resp1_re)
  
  oracles <- vector(mode='list', length = 400)
  names(oracles) <- names(out1_resp1_re)
  oracletimes <- vector(mode = 'list', length = 400)
  names(oracletimes) <- names(out1_resp1_re)
  
  for(i in 1:400) {
    df <- out1_resp1_re[[i]] |> 
      mutate(race = as.factor(race), smoke = as.factor(smoke)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrs[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ Hg + Sb + 
                            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))), 
                          data = df)
    
    if(i <= 200) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*(race == 2)) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      chems_oracle[[i]] <- lm(y ~ Hg + Sb +
                                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                                I(Hg*(race == 2)), data = df)
    } else if (i <= 400) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*(race == 5)) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      chems_oracle[[i]] <- lm(y ~ Hg + Sb +
                                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                                I(Hg*(race == 5)), data = df)
    } 
  }
  return(list(mlrs, mlrtimes, oracles, oracletimes, chems_only, chems_oracle))
}

mlrsm_job <- slurm_call(
  run_mlr_sm_re, 
  global_objects = c('out1_resp1_re'),
  jobname = 'mlr_sm_re')

# check on status of job
get_job_status(mlrsm_job)$completed
system('squeue --me')

# get output
mlr_sm <- read_rds("_rslurm_mlr_sm_re/results_0.RDS")
# mlr_sm <- get_slurm_out(sjob5)
mlr_mods <- mlr_sm[[1]]
mlr_times <- mlr_sm[[2]]
oracle_mods <- mlr_sm[[3]]
oracle_times <- mlr_sm[[4]]
chem_mods <- mlr_sm[[5]]
chems_oracle <- mlr_sm[[6]]
write_rds(mlr_mods, "sim/mlr_mods_sm_re.RDS")
write_rds(mlr_times, "sim/mlr_mods_sm_re_times.RDS")
write_rds(oracle_mods, "sim/oracle_mods_sm_re.RDS")
write_rds(oracle_times, "sim/oracle_mods_sm_re_times.RDS")
write_rds(chem_mods, "sim/chem_mods_re_sm.RDS")
write_rds(chems_oracle, "sim/chem_oracle_re_sm.RDS")

### larger sample size

out2_resp1_re <- read_rds("sim/sim_resp_lg_re.RDS")
run_mlr_lg_re <- function() {
  # initialize vectors
  mlrs <- vector(mode='list', length = 400)
  names(mlrs) <- names(out2_resp1_re)
  mlrtimes <- vector(mode = 'list', length = 400)
  names(mlrtimes) <- names(out2_resp1_re)
  
  chems_only <- vector(mode = 'list', length = 400)
  names(chems_only) <- names(out2_resp1_re)
  chems_oracle <- vector(mode = 'list', length = 400) 
  names(mlrs) <- names(out2_resp1_re)
  
  oracles <- vector(mode='list', length = 400)
  names(oracles) <- names(out2_resp1_re)
  oracletimes <- vector(mode = 'list', length = 400)
  names(oracletimes) <- names(out2_resp1_re)
  
  for(i in 1:400) {
    df <- out2_resp1_re[[i]] |> 
      mutate(race = as.factor(race), smoke = as.factor(smoke)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrs[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ Hg + Sb + 
                            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))), 
                          data = df)
    
    if(i <= 200) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*(race == 2)) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      chems_oracle[[i]] <- lm(y ~ Hg + Sb +
                                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                                I(Hg*(race == 2)), data = df)
    } else if (i <= 400) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ Hg + Sb +
                           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                           I(Hg*(race == 5)) +
                           age + bmi + race + smoke, data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      chems_oracle[[i]] <- lm(y ~ Hg + Sb +
                                I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
                                I(Hg*(race == 5)), data = df)
    } 
  }
  return(list(mlrs, mlrtimes, oracles, oracletimes, chems_only, chems_oracle))
}

mlrlg_job <- slurm_call(
  run_mlr_lg_re, 
  global_objects = c('out2_resp1_re'),
  jobname = 'mlr_lg_re')

# check on status of job
get_job_status(mlrlg_job)$completed
system('squeue --me')

# get output
mlr_lg <- read_rds("_rslurm_mlr_lg_re/results_0.RDS")
mlr_modl <- mlr_lg[[1]]
mlr_timel <- mlr_lg[[2]]
oracle_modl <- mlr_lg[[3]]
oracle_timel <- mlr_lg[[4]]
chem_modl <- mlr_lg[[5]]
chem_oraclel <- mlr_lg[[6]]
write_rds(mlr_modl, "sim/mlr_mods_lg_re.RDS")
write_rds(mlr_timel, "sim/mlr_mods_lg_re_times.RDS")
write_rds(oracle_modl, "sim/oracle_mods_lg_re.RDS")
write_rds(oracle_timel, "sim/oracle_mods_lg_re_times.RDS")
write_rds(chem_modl, "sim/chem_mods_re_lg.RDS")
write_rds(chem_oraclel, "sim/chem_oracle_re_sm.RDS")

###########
# run stratified bkmr
###########

out1_resp1_re <- read_rds("sim/sim_resp_sm_re.RDS")

run_bkmr_sm_re <- function(vector) {
  bkmr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  if(!dir.exists("times")) {
    dir.create("times")
  }
  
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out1_resp1_re[[i]]
    
    # for each race level, run bmkr
    list_times <- vector(mod = "list", length = 6)
    list_mods <- vector(mod = "list", length = 6)
    
    set.seed(0)
    for(j in 1:5) {
      df <- df_full[df_full$race == j, ]
      Z <- df |> 
        select(As:Sn)
      X <- df |>
        select(smoke:bmi)
      y <- df$y
      
      # fit model and save time
      start.time <- Sys.time()
      mod <- tryCatch({
        kmbayes(y = y, Z = Z, X = X, 
                iter = 50000, verbose = FALSE, varsel = TRUE)
      }, error = function(e) {
        NA
      })
      end.time <- Sys.time()
      list_times[[j]] <- end.time - start.time
      list_mods[[j]] <- mod
    }
    
    # combine 1-3 r/e
    df <- df_full[df_full$race %in% c(1, 2, 3), ]
    Z <- df |> 
      select(As:Sn)
    X <- df |>
      select(smoke:bmi)
    y <- df$y
    
    # fit model and save time
    start.time <- Sys.time()
    mod <- tryCatch({
      kmbayes(y = y, Z = Z, X = X, 
              iter = 50000, verbose = FALSE, varsel = TRUE)
    }, error = function(e) {
      NA
    })
    end.time <- Sys.time()
    list_times[[6]] <- end.time - start.time
    list_mods[[6]] <- mod
    
    
    bkmr_times[[i]] <- list_times
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bkmr_sm_", names(out1_resp1_re)[i], "_", i, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bkmr_smf_", names(out1_resp1_re)[i], "_", i, ".RDS"))
    rm(mod)
  }
  # write_rds(bkmr_times, file = "bkmr_sm_times.RDS")
  return(bkmr_times)
}

ujob01 <- slurm_call(
  run_bkmr_sm_re, params = list(vector = 1:100),
  global_objects = c('out1_resp1_re'),
  jobname = 'ksmre01',
  slurm_options = list(mem = '8G'))

ujob02 <- slurm_call(
  run_bkmr_sm_re, params = list(vector = 101:200),
  global_objects = c('out1_resp1_re'),
  jobname = 'ksmre02',
  slurm_options = list(mem = '8G'))

ujob03 <- slurm_call(
  run_bkmr_sm_re, params = list(vector = 201:300),
  global_objects = c('out1_resp1_re'),
  jobname = 'ksmre03',
  slurm_options = list(mem = '8G'))

ujob04 <- slurm_call(
  run_bkmr_sm_re, params = list(vector = 301:400),
  global_objects = c('out1_resp1_re'),
  jobname = 'ksmre04',
  slurm_options = list(mem = '8G'))

# larger size
out2_resp1_re <- read_rds("sim/sim_resp_lg_re.RDS")
run_bkmr_lg_re <- function(vector) {
  bkmr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  if(!dir.exists("times")) {
    dir.create("times")
  }
  
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out2_resp1_re[[i]]
    
    # for each race level, run bmkr
    list_times <- vector(mod = "list", length = 5)
    list_mods <- vector(mod = "list", length = 5)
    
    set.seed(0)
    for(j in 1:5) {
      df <- df_full[df_full$race == j, ]
      Z <- df |> 
        select(As:Sn)
      X <- df |>
        select(smoke:bmi)
      y <- df$y
      
      # fit model and save time
      start.time <- Sys.time()
      mod <- kmbayes(y = y, Z = Z, X = X, 
                     iter = 50000, verbose = FALSE, varsel = TRUE)
      end.time <- Sys.time()
      list_times[[j]] <- end.time - start.time
      list_mods[[j]] <- mod
    }
    
    bkmr_times[[i]] <- list_times
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bkmr_lg_", names(out2_resp1_re)[i], "_", i, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bkmr_lgf_", names(out2_resp1_re)[i], "_", i, ".RDS"))
    rm(mod)
  }
  # write_rds(bkmr_times, file = "bkmr_lg_times_re.RDS")
  return(bkmr_times)
}

tjob01 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 1:100),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre01',
  slurm_options = list(mem = '8G'))

tjob02 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 101:200),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre02',
  slurm_options = list(mem = '8G'))

tjob03 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 201:300),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre03',
  slurm_options = list(mem = '8G'))

tjob04 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 301:400),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre04',
  slurm_options = list(mem = '8G'))

tjob01_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(1:100)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre01b',
  slurm_options = list(mem = '8G'))

tjob02_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(101:200)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre02b',
  slurm_options = list(mem = '8G'))

tjob03_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(201:300)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre03b',
  slurm_options = list(mem = '8G'))

tjob04_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(301:400)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre04b',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')

###########
# run stratified bsr
###########

out1_resp1_re <- read_rds("sim/sim_resp_sm_re.RDS")

# one <- out1_resp1_re[[1]]
# df <- one[one$race == 1, ]
# X <- df |> 
#   select(As:Sn) |> 
#   as.matrix.data.frame()
# C <- df |> 
#   select(smoke:bmi) |> 
#   as.matrix.data.frame()
# Y <- df$y
# set.seed(0)
# mod1 <- NLint(Y = Y, X = X, C = C, 
#               nIter = 500, nBurn = 250, ns = 1)


run_bsr_sm_re <- function(vector) {
  bsr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  if(!dir.exists("times")) {
    dir.create("times")
  }
  
  # extract names of files
  list_files <- list.files("mods", full.names = TRUE)
  nums <- as.numeric(sub(".+_(.+)d.+", "\\1", list_files))
  if(length(nums) == 0) {
    final <- 0
    starting <- min(vector)
  } else {
    final <- max(nums)
    starting <- final + 1
  }
  if(file.exists("final.txt")) {
    write(paste0("final ran = ", final, ", starting at ", starting), "final.txt", append = T)
  } else {
    writeLines(paste0("final ran = ", final, ", starting at ", starting), "final.txt")
  }
  
  
  for(i in starting:max(vector)) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out1_resp1_re[[i]]
    
    # for each race level, run bsr
    list_times <- vector(mod = "list", length = 6) 
    list_mods <- vector(mod = "list", length = 6)
    
    set.seed(0)
    for(j in 1:5) {
      df <- df_full[df_full$race == j, ]
      X <- df |> 
        select(As:Sn) |> 
        as.matrix.data.frame()
      C <- df |> 
        select(smoke:bmi) |> 
        as.matrix.data.frame()
      Y <- df$y
      
      # waic for choosing df
      list_times_small <- vector(mode = "list", length = 2)
      
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
      list_times_small[[1]] <- end.time - start.time
      
      ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
      
      # fit model and save time
      start.time <- Sys.time()
      mod <- tryCatch({
        NLint(Y = Y, X = X, C = C, 
              nIter = 50000, nBurn = 25000, ns = ind)
      }, error = function(e) {
        NA
      })
      end.time <- Sys.time()
      list_times_small[[2]] <- end.time - start.time
      
      bsr_times[[j]] <- list_times_small
      list_mods[[j]] <- mod
    }
    
    # combine 1-3 r/e
    df <- df_full[df_full$race %in% c(1, 2, 3), ]
    X <- df |> 
      select(As:Sn) |> 
      as.matrix.data.frame()
    C <- df |> 
      select(smoke:bmi) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # waic for choosing df
    list_times_small <- vector(mode = "list", length = 2)
    
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
    list_times_small[[1]] <- end.time - start.time
    
    ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
    
    # fit model and save time
    start.time <- Sys.time()
    mod <- tryCatch({
      NLint(Y = Y, X = X, C = C, 
            nIter = 50000, nBurn = 25000, ns = ind)
    }, error = function(e) {
      NA
    })
    end.time <- Sys.time()
    list_times_small[[2]] <- end.time - start.time
    
    list_times[[6]] <- list_times_small
    list_mods[[6]] <- mod
    
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bsr_sm_", names(out1_resp1_re)[i], "_", i, 
                       "df", ind, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bsr_sm_", names(out1_resp1_re)[i], "_", i, 
                       "df", ind, ".RDS"))
    
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  # write_rds(bsr_times, file = "bsr_smf_times.RDS")
  return(bsr_times)
}

vjob01 <- slurm_call(
  run_bsr_sm_re, params = list(vector = 1:100),
  global_objects = c('out1_resp1_re'),
  jobname = 'ssmre01',
  slurm_options = list(mem = '8G'))

vjob02 <- slurm_call(
  run_bsr_sm_re, params = list(vector = 101:200),
  global_objects = c('out1_resp1_re'),
  jobname = 'ssmre02',
  slurm_options = list(mem = '8G'))

vjob03 <- slurm_call(
  run_bsr_sm_re, params = list(vector = 201:300),
  global_objects = c('out1_resp1_re'),
  jobname = 'ssmre03',
  slurm_options = list(mem = '8G'))

vjob04 <- slurm_call(
  run_bsr_sm_re, params = list(vector = 301:400),
  global_objects = c('out1_resp1_re'),
  jobname = 'ssmre04',
  slurm_options = list(mem = '8G'))

# start sending them in individually
for (i in 1:100) {
  slurm_call(run_bsr_sm_re, params = list(vector = i), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre01_', i))
}

for (i in 101:200) {
  slurm_call(run_bsr_sm_re, params = list(vector = i), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre02_', i))
}

for (i in 201:300) {
  slurm_call(run_bsr_sm_re, params = list(vector = i), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre03_', i))
}

for (i in 301:400) {
  slurm_call(run_bsr_sm_re, params = list(vector = i), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre04_', i))
}

# larger size
out2_resp1_re <- read_rds("sim/sim_resp_lg_re.RDS")

run_bsr_lg_re <- function(vector) {
  bsr_times <- vector(mode = "list", length = length(vector))
  
  # create folder for model output
  if(!dir.exists("mods")) {
    dir.create("mods")
  }
  if(!dir.exists("times")) {
    dir.create("times")
  }
  
  # extract names of files
  list_files <- list.files("mods", full.names = TRUE)
  nums <- as.numeric(sub(".+_(.+)d.+", "\\1", list_files))
  if(length(nums) == 0) {
    final <- 0
    starting <- min(vector)
  } else {
    final <- max(nums)
    starting <- final + 1
  }
  if(file.exists("final.txt")) {
    write(paste0("final ran = ", final, ", starting at ", starting), "final.txt", append = T)
  } else {
    writeLines(paste0("final ran = ", final, ", starting at ", starting), "final.txt")
  }
  
  
  for(i in starting:max(vector)) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out2_resp1_re[[i]]
    
    # for each race level, run bsr
    list_times <- vector(mod = "list", length = 6) 
    list_mods <- vector(mod = "list", length = 6)
    
    set.seed(0)
    for(j in 1:5) {
      df <- df_full[df_full$race == j, ]
      X <- df |> 
        select(As:Sn) |> 
        as.matrix.data.frame()
      C <- df |> 
        select(smoke:bmi) |> 
        as.matrix.data.frame()
      Y <- df$y
      
      # waic for choosing df
      list_times_small <- vector(mode = "list", length = 2)
      
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
      list_times_small[[1]] <- end.time - start.time
      
      ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
      
      # fit model and save time
      start.time <- Sys.time()
      mod <- tryCatch({
        NLint(Y = Y, X = X, C = C, 
              nIter = 50000, nBurn = 25000, ns = ind)
      }, error = function(e) {
        NA
      })
      end.time <- Sys.time()
      list_times_small[[2]] <- end.time - start.time
      
      bsr_times[[j]] <- list_times_small
      list_mods[[j]] <- mod
    }
    
    # combine 1-3 r/e
    df <- df_full[df_full$race %in% c(1, 2, 3), ]
    X <- df |> 
      select(As:Sn) |> 
      as.matrix.data.frame()
    C <- df |> 
      select(smoke:bmi) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # waic for choosing df
    list_times_small <- vector(mode = "list", length = 2)
    
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
    list_times_small[[1]] <- end.time - start.time
    
    ind <- which.min(c(mod1$waic, mod2$waic, mod3$waic, mod4$waic))
    
    # fit model and save time
    start.time <- Sys.time()
    mod <- tryCatch({
      NLint(Y = Y, X = X, C = C, 
            nIter = 50000, nBurn = 25000, ns = ind)
    }, error = function(e) {
      NA
    })
    end.time <- Sys.time()
    list_times_small[[2]] <- end.time - start.time
    
    list_times[[6]] <- list_times_small
    list_mods[[6]] <- mod
    
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bsr_lg_", names(out2_resp1_re)[i], "_", i, 
                       "df", ind, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bsr_lg_", names(out2_resp1_re)[i], "_", i, 
                       "df", ind, ".RDS"))
    
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  # write_rds(bsr_times, file = "bsr_lgf_times.RDS")
  return(bsr_times)
}

wjob01 <- slurm_call(
  run_bsr_lg_re, params = list(vector = 1:100),
  global_objects = c('out2_resp1_re'),
  jobname = 'slgre01',
  slurm_options = list(mem = '8G'))

wjob02 <- slurm_call(
  run_bsr_lg_re, params = list(vector = 101:200),
  global_objects = c('out2_resp1_re'),
  jobname = 'slgre02',
  slurm_options = list(mem = '8G'))

wjob03 <- slurm_call(
  run_bsr_lg_re, params = list(vector = 201:300),
  global_objects = c('out2_resp1_re'),
  jobname = 'slgre03',
  slurm_options = list(mem = '8G'))

wjob04 <- slurm_call(
  run_bsr_lg_re, params = list(vector = 301:400),
  global_objects = c('out2_resp1_re'),
  jobname = 'slgre04',
  slurm_options = list(mem = '8G'))

# start sending them in individually
for (i in 1:100) {
  slurm_call(run_bsr_lg_re, params = list(vector = i), 
             global_objects = c('out2_resp1_re'), 
             jobname = paste0('slgre01_', i))
}

for (i in 101:200) {
  slurm_call(run_bsr_lg_re, params = list(vector = i), 
             global_objects = c('out2_resp1_re'), 
             jobname = paste0('slgre02_', i))
}

for (i in 201:300) {
  slurm_call(run_bsr_lg_re, params = list(vector = i), 
             global_objects = c('out2_resp1_re'), 
             jobname = paste0('slgre03_', i))
}

for (i in 301:400) {
  slurm_call(run_bsr_lg_re, params = list(vector = i), 
             global_objects = c('out2_resp1_re'), 
             jobname = paste0('slgre04_', i))
}

# check on status of job
system('squeue --me')

system('squeue --format="%.6i %.6P %.12j %.11M %.11l %.4D %.8R" --me')
