library(tidyverse)
library(rslurm)
library(bkmr)
library(NLinteraction)


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
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 == 3 ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA
      )) |> 
      mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrs[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimes[[i]] <- end.time - start.time

    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA +
                            I(1/(1+exp(-4*LBX194LA))) +
                            I(1/(1+exp(-4*LBXPCBLA))) +
                            I(LBXF04LA^2) + LBXF04LA ,
                          data = df)
    
    if(i <= 200) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                           I(1/(1+exp(-4*LBX194LA))) + 
                           I(1/(1+exp(-4*LBXPCBLA))) + 
                           I(LBXF04LA^2) + LBXF04LA +
                           I(LBXD05LA*(RIDRETH1 == 4)) +
                           BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                           LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                           RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                         data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      
      chems_oracle[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                              I(1/(1+exp(-4*LBX194LA))) + 
                              I(1/(1+exp(-4*LBXPCBLA))) + 
                              I(LBXF04LA^2) + LBXF04LA +
                                I(LBXD05LA*(RIDRETH1 == 4)),
                            data = df)
    } else if (i <= 400) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                           I(1/(1+exp(-4*LBX194LA))) + 
                           I(1/(1+exp(-4*LBXPCBLA))) + 
                           I(LBXF04LA^2) + LBXF04LA +
                           I(LBXD05LA*(RIDRETH1 == 3)) +
                           BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                           LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                           RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                         data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      
      chems_oracle[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                                I(1/(1+exp(-4*LBX194LA))) + 
                                I(1/(1+exp(-4*LBXPCBLA))) + 
                                I(LBXF04LA^2) + LBXF04LA +
                                I(LBXD05LA*(RIDRETH1 == 3)),
                              data = df)
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
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 %in% c(3) ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA
      )) |> 
      mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
      select(-sim)
    
    start.time <- Sys.time()
    mlrs[[i]] <- lm(y ~ ., data = df)
    end.time <- Sys.time()
    mlrtimes[[i]] <- end.time - start.time
    
    chems_only[[i]] <- lm(y ~ LBXD05LA + LBX074LA +
                            I(1/(1+exp(-4*LBX194LA))) +
                            I(1/(1+exp(-4*LBXPCBLA))) +
                            I(LBXF04LA^2) + LBXF04LA ,
                          data = df)
    
    if(i <= 200) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                           I(1/(1+exp(-4*LBX194LA))) + 
                           I(1/(1+exp(-4*LBXPCBLA))) + 
                           I(LBXF04LA^2) + LBXF04LA +
                           I(LBXD05LA*(RIDRETH1 == 4)) +
                           BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                           LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                           RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                         data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      
      chems_oracle[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                                I(1/(1+exp(-4*LBX194LA))) + 
                                I(1/(1+exp(-4*LBXPCBLA))) + 
                                I(LBXF04LA^2) + LBXF04LA +
                                I(LBXD05LA*(RIDRETH1 == 4)),
                              data = df)
    } else if (i <= 400) {
      start.time <- Sys.time()
      oracles[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                           I(1/(1+exp(-4*LBX194LA))) + 
                           I(1/(1+exp(-4*LBXPCBLA))) + 
                           I(LBXF04LA^2) + LBXF04LA +
                           I(LBXD05LA*(RIDRETH1 == 3)) +
                           BMXBMI + LBXCOT + LBXWBCSI + LBXLYPCT + LBXMOPCT +
                           LBXNEPCT + LBXEOPCT + LBXBAPCT + 
                           RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2, 
                         data = df)
      end.time <- Sys.time()
      oracletimes[[i]] <- end.time - start.time
      
      chems_oracle[[i]] <- lm(y ~ LBXD05LA + LBX074LA + 
                                I(1/(1+exp(-4*LBX194LA))) + 
                                I(1/(1+exp(-4*LBXPCBLA))) + 
                                I(LBXF04LA^2) + LBXF04LA +
                                I(LBXD05LA*(RIDRETH1 == 3)),
                              data = df)
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
    df_full <- out1_resp1_re[[i]] |> 
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 %in% c(3) ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA
      ))
    
    # for each race level, run bmkr
    list_times <- vector(mod = "list", length = 3)
    list_mods <- vector(mod = "list", length = 3)
    
    set.seed(0)
    for(j in 1:3) {
      df <- df_full[df_full$RIDRETH1 == j, ]
      Z <- df |> select(LBX074LA:LBXF08LA)
      X <- df |> 
        mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
               DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
               DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
               DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
        select(-DMDEDUC2) |> 
        select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT)
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
  
  # list_files <- list.files("mods", full.names = TRUE)
  # nums <- as.numeric(sub(".*_(\\d+)\\.RDS$", "\\1", list_files))
  # if(length(nums) == 0) {
  #   final <- 0
  #   starting <- min(vector)
  # } else {
  #   final <- max(nums)
  #   starting <- final + 1
  # }
  
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out2_resp1_re[[i]] |> 
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 %in% c(3) ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA
      ))
    
    # for each race level, run bmkr
    list_times <- vector(mod = "list", length = 3)
    list_mods <- vector(mod = "list", length = 3)
    
    set.seed(0)
    for(j in 1:3) {
      df <- df_full[df_full$RIDRETH1 == j, ]
      Z <- df |> select(LBX074LA:LBXF08LA)
      X <- df |> 
        mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
               DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
               DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
               DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
        select(-DMDEDUC2) |> 
        select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT)
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
    
    # # combine 2, 5 r/e
    # df <- df_full[df_full$RIDRETH1 %in% c(2, 5), ]
    # Z <- df |> select(LBX074LA:LBXF08LA)
    # X <- df |> 
    #   mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
    #          DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
    #          DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
    #          DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
    #   select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT)
    # y <- df$y 
    # 
    # # fit model and save time
    # start.time <- Sys.time()
    # mod <- tryCatch({
    #   kmbayes(y = y, Z = Z, X = X, 
    #           iter = 50000, verbose = FALSE, varsel = TRUE)
    # }, error = function(e) {
    #   NA
    # })
    # end.time <- Sys.time()
    # list_times[[6]] <- end.time - start.time
    # list_mods[[6]] <- mod
    
    bkmr_times[[i]] <- list_times
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bkmr_lg_", names(out2_resp1_re)[i], "_", i, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bkmr_lg_", names(out2_resp1_re)[i], "_", i, ".RDS"))
    rm(mod)
  }
  # write_rds(bkmr_times, file = "bkmr_lg_times_re.RDS")
  return(bkmr_times)
}

tjob01 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 1:50),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre01',
  slurm_options = list(mem = '8G'))

tjob02 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 101:150),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre02',
  slurm_options = list(mem = '8G'))

tjob03 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 201:250),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre03',
  slurm_options = list(mem = '8G'))

tjob04 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = 301:350),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre04',
  slurm_options = list(mem = '8G'))

tjob01_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(51:100)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre01b',
  slurm_options = list(mem = '8G'))

tjob02_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(151:200)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre02b',
  slurm_options = list(mem = '8G'))

tjob03_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(251:300)),
  global_objects = c('out2_resp1_re'),
  jobname = 'klgre03b',
  slurm_options = list(mem = '8G'))

tjob04_2 <- slurm_call(
  run_bkmr_lg_re, params = list(vector = rev(351:400)),
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


run_bsr_sm_re <- function(vector, re) {
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
  check <- any(!grepl("_re", list_files)) && length(list_files) > 0
  if (check) {
    return("Exiting function: '_re' not found in at least one of the files.")
  }
  
  # get "re" files
  filtered_files <- list_files[grepl("re", list_files)]
  
  # if filtered_files is not empty
  if (length(filtered_files) > 0) {
    # extract numbers between "re" and ".RDS"
    subtract <- as.numeric(gsub(".*re(\\d+)\\.RDS.*", "\\1", filtered_files))
    
    # return final vector to cycle through
    vector <- setdiff(vector, subtract)
    
  } 
  
  # if(grepl("pone", list_files)) {
  for(i in vector) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df_full <- out1_resp1_re[[i]] |>
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1,
        RIDRETH1 %in% c(3) ~ 2,
        RIDRETH1 %in% c(4, 5) ~ 3,
        .default = NA
      ))
    
    # for each race level, run bsr
    list_times <- vector(mod = "list", length = 3)
    list_mods <- vector(mod = "list", length = 3)
    
    set.seed(0)
    
    for (j in re) {
      df <- df_full[df_full$RIDRETH1 == j, ]
      X <- df |> select(LBX074LA:LBXF08LA) |> as.matrix.data.frame()
      C <- df |>
        mutate(
          DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0),
          DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0),
          DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0),
          DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)
        ) |>
        select(-DMDEDUC2) |>
        select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT) |>
        as.matrix.data.frame()
      Y <- df$y
      
      # waic for choosing df
      list_times_small <- vector(mode = "list", length = 2)
      
      start.time <- Sys.time()
      mod1 <- tryCatch({
        NLint(Y = Y, X = X, C = C,
          nIter = 5000, nBurn = 2500, ns = 1
        )
      }, error = function(e) {
        NA
      })
      mod2 <- tryCatch({
        NLint(Y = Y, X = X, C = C,
          nIter = 5000, nBurn = 2500, ns = 2
        )
      }, error = function(e) {
        NA
      })
      mod3 <- tryCatch({
        NLint(Y = Y, X = X, C = C,
          nIter = 5000, nBurn = 2500, ns = 3
        )
      }, error = function(e) {
        NA
      })
      mod4 <- tryCatch({
        NLint(Y = Y, X = X, C = C,
          nIter = 5000, nBurn = 2500, ns = 4
        )
      }, error = function(e) {
        NA
      })
      end.time <- Sys.time()
      list_times_small[[1]] <- end.time - start.time
      
      ind <- which.min(c(
        ifelse(is.na(mod1), NA, mod1$waic),
        ifelse(is.na(mod2), NA, mod2$waic),
        ifelse(is.na(mod3), NA, mod3$waic),
        ifelse(is.na(mod4), NA, mod4$waic)
      ))
      
      # fit model and save time
      start.time <- Sys.time()
      mod <- tryCatch({
        NLint(Y = Y, X = X, C = C,
          nIter = 50000, nBurn = 25000, ns = ind
        )
      }, error = function(e) {
        NA
      })
      end.time <- Sys.time()
      list_times_small[[2]] <- end.time - start.time
      
      list_times[[j]] <- list_times_small
      list_mods[[j]] <- mod
      write_rds(mod, file =
                  paste0("mods/bsr_sm_", names(out1_resp1_re)[i], "_",
                    i, "_re", j, ".RDS"
                  ))
      write_rds(list_times_small,
                file =
                  paste0("times/bsr_sm_", names(out1_resp1_re)[i], "_", i, ".RDS"))
    }
    
    
    # save model and remove from memory
    write_rds(list_mods, file =
                paste0("mods/bsr_sm_", names(out1_resp1_re)[i], "_", i, ".RDS"))
    write_rds(list_times, file =
                paste0("times/bsr_sm_", names(out1_resp1_re)[i], "_", i, ".RDS"))
    bsr_times[i] <- list_times
  }
  
  write_rds(bsr_times, file = "bsr_smf_times.RDS")
  return(bsr_times)
}

# start sending them in individually
for (i in 1:100) {
  slurm_call(run_bsr_sm_re, params = list(vector = i, re = 1:3), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre01_', i), 
             slurm_options = list(mem = '16G'))
}

for (i in 101:200) {
  slurm_call(run_bsr_sm_re, params = list(vector = i, re = 1:3), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre02_', i), 
             slurm_options = list(mem = '16G'))
}

for (i in 201:300) {
  slurm_call(run_bsr_sm_re, params = list(vector = i, re = 1:3), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre03_', i), 
             slurm_options = list(mem = '16G'))
}

for (i in 301:400) {
  slurm_call(run_bsr_sm_re, params = list(vector = i, re = 1:3), 
             global_objects = c('out1_resp1_re'), 
             jobname = paste0('ssmre04_', i), 
             slurm_options = list(mem = '16G'))
}

# # still running 
# sr <- c(371, 352, 313, 314, 288, 262, 270, 256, 231, 236, 224, 136, 117, 108, 107, 35, 36, 2)
# 
# # for still running, fit 5 and 2/5
# for (i in sr) {
#   header <- ceiling(i/100)
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0("ssmre0", header, "_", i))
# }
# 
# # for those that fixed after first re-try, fit 2/5 re
# for (i in setdiff(1:100, sr)) {
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0('ssmre01_', i))
# }
# 
# for (i in setdiff(101:200, sr)) {
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0('ssmre02_', i))
# }
# 
# for (i in setdiff(201:300, sr)) {
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0('ssmre03_', i))
# }
# 
# for (i in setdiff(301:400, sr)) {
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0('ssmre04_', i))
# }
# 
# # for those which failed b/c of category 2
# wh <- c(98, 108, 137, 174)
# for (i in wh) {
#   header <- ceiling(i/100)
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0("ssmre0", header, "_", i))
# }
# 
# # for 35 and 256, grp 5 won't coverage
# g5 <- c(35, 256)
# for (i in g5) {
#   header <- ceiling(i/100)
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0("ssmre0", header, "_", i))
# }
# 
# # for those which still don't have first four
# f4 <- c(108, 136, 224, 288)
# for (i in f4) {
#   header <- ceiling(i/100)
#   slurm_call(run_bsr_sm_re, params = list(vector = i), 
#              global_objects = c('out1_resp1_re'), 
#              jobname = paste0("ssmre0", header, "_", i))
# }
# 
# job1 <- slurm_call(run_bsr_sm_re, params = list(vector = 108), 
#            global_objects = c('out1_resp1_re'), 
#            jobname = paste0("ssmre0", 1, "_", i))

##############
# need to combine all the split up ones into one after finishes running
##############

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
    df_full <- out2_resp1_re[[i]] |> 
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 %in% c(3) ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA
      ))
    
    # for each race level, run bsr
    list_times <- vector(mod = "list", length = 3) 
    list_mods <- vector(mod = "list", length = 3)
    
    set.seed(0)
    for(j in 1:3) {
      df <- df_full[df_full$RIDRETH1 == j, ]
      X <- df |> select(LBX074LA:LBXF08LA) |> as.matrix.data.frame()
      C <- df |> 
        mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
               DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
               DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
               DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
        select(-DMDEDUC2) |> 
        select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT) |> 
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
    
    
    # save model and remove from memory
    write_rds(list_mods, file = 
                paste0("mods/bsr_lg_", names(out2_resp1_re)[i], "_", i, 
                       ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bsr_lg_", names(out2_resp1_re)[i], "_", i, 
                       ".RDS"))
    
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  # write_rds(bsr_times, file = "bsr_lgf_times.RDS")
  return(bsr_times)
}

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


# qgcompint ---------------------------------------------------------------

library(qgcompint)

# smaller size

out1_resp1_re <- read_rds("sim/sim_resp_sm_re.RDS")

run_qgc_sm_re <- function(vector) {
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
    df_full <- out1_resp1_re[[i]] |> 
      mutate(
        RIDRETH1 = case_when(
          RIDRETH1 %in% c(1, 2) ~ 1, 
          RIDRETH1 %in% c(3) ~ 2, 
          RIDRETH1 %in% c(4, 5) ~ 3, 
          .default = NA)
      )
    
    df <- df_full |> 
      mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
      as.data.frame() |> 
      rename_all(tolower)
    
    # times
    list_times <- vector(mode = "list", length = 2)
    
    # fit bootstrap
    start.time <- Sys.time()
    mod1 <- tryCatch({qgcomp.emm.boot(
      f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
        lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
        lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
        bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
        lbxbapct + ridageyr + riagendr + dmdeduc2, 
      data = df, 
      expnms = names(df)[1:18], 
      emmvar = "ridreth1", 
      family = gaussian(), 
      q = 4, 
      alpha = 0.05, 
      B = 10000, 
      seed = 0
    )}, error = function(e) {
      NA
    })
    end.time <- Sys.time()
    list_times[[1]] <- end.time - start.time
    
    # fit no boostrap
    start.time <- Sys.time()
    mod2 <- tryCatch({
      qgcomp.emm.noboot(
      f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
        lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
        lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
        bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
        lbxbapct + ridageyr + riagendr + dmdeduc2, 
      data = df, 
      expnms = names(df)[1:18], 
      emmvar = "ridreth1", 
      family = gaussian(), 
      q = 4, 
      alpha = 0.05
      )}, error = function(e) {
      NA
    })
    end.time <- Sys.time()
    list_times[[2]] <- end.time - start.time
    
    # # collapse 2/5
    # dfsmall <- df |> 
    #   mutate(ridreth1 = ifelse(ridreth1 == 5, 2, ridreth1), 
    #          ridreth1 = as.factor(ridreth1))
    # 
    # # fit bootstrap
    # start.time <- Sys.time()
    # mod3 <- tryCatch({qgcomp.emm.boot(
    #   f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
    #     lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
    #     lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
    #     bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
    #     lbxbapct + ridageyr + riagendr + dmdeduc2, 
    #   data = df, 
    #   expnms = names(df)[1:18], 
    #   emmvar = "ridreth1", 
    #   family = gaussian(), 
    #   q = 4, 
    #   alpha = 0.05, 
    #   B = 10000, 
    #   seed = 0
    # )}, error = function(e) {
    #   NA
    # })
    # end.time <- Sys.time()
    # list_times[[3]] <- end.time - start.time
    # 
    # # fit no boostrap
    # start.time <- Sys.time()
    # mod4 <- tryCatch({
    #   qgcomp.emm.noboot(
    #     f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
    #       lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
    #       lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
    #       bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
    #       lbxbapct + ridageyr + riagendr + dmdeduc2, 
    #     data = df, 
    #     expnms = names(df)[1:18], 
    #     emmvar = "ridreth1", 
    #     family = gaussian(), 
    #     q = 4, 
    #     alpha = 0.05
    #   )}, error = function(e) {
    #     NA
    #   })
    # end.time <- Sys.time()
    # list_times[[4]] <- end.time - start.time
    # 
    # save models
    list_mods <- list(mod1, mod2)
    
    # write out
    write_rds(list_mods, file = 
                paste0("mods/qgc_sm_", names(out1_resp1_re)[i], "_", i, 
                       "df.RDS"))
    write_rds(list_times, file = 
                paste0("times/qgc_sm_", names(out1_resp1_re)[i], "_", i, 
                       "df.RDS"))
    
  }
  
}

# run for all 
qgjob01 <- slurm_call(
  run_qgc_sm_re, params = list(vector = 1:100),
  global_objects = c('out1_resp1_re'),
  jobname = 'qsmre01',
  slurm_options = list(mem = '8G'))

qgjob02 <- slurm_call(
  run_qgc_sm_re, params = list(vector = 101:200),
  global_objects = c('out1_resp1_re'),
  jobname = 'qsmre02',
  slurm_options = list(mem = '8G'))

qgjob03 <- slurm_call(
  run_qgc_sm_re, params = list(vector = 201:300),
  global_objects = c('out1_resp1_re'),
  jobname = 'qsmre03',
  slurm_options = list(mem = '8G'))

qgjob04 <- slurm_call(
  run_qgc_sm_re, params = list(vector = 301:400),
  global_objects = c('out1_resp1_re'),
  jobname = 'qsmre04',
  slurm_options = list(mem = '8G'))

# larger size

out2_resp1_re <- read_rds("sim/sim_resp_lg_re.RDS")

library(qgcompint)

run_qgc_lg_re <- function(vector) {
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
    df_full <- out2_resp1_re[[i]] |> 
      mutate(RIDRETH1 = case_when(
        RIDRETH1 %in% c(1, 2) ~ 1, 
        RIDRETH1 %in% c(3) ~ 2, 
        RIDRETH1 %in% c(4, 5) ~ 3, 
        .default = NA))
    
    df <- df_full |> 
      mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
      as.data.frame() |> 
      rename_all(tolower)
    
    # times
    list_times <- vector(mode = "list", length = 2)
    
    # fit bootstrap
    start.time <- Sys.time()
    mod1 <- tryCatch({
      qgcomp.emm.boot(
        f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
          lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
          lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
          bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
          lbxbapct + ridageyr + riagendr + dmdeduc2, 
        data = df, 
        expnms = names(df)[1:18], 
        emmvar = "ridreth1", 
        family = gaussian(), 
        q = 4, 
        alpha = 0.05, 
        B = 10000, 
        seed = 0
      )
    }, error = function(e) {
      NA
    })
    
    end.time <- Sys.time()
    list_times[[1]] <- end.time - start.time
    
    # fit no boostrap
    start.time <- Sys.time()
    mod2 <- tryCatch({
      qgcomp.emm.noboot(
        f = y ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
          lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
          lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
          bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
          lbxbapct + ridageyr + riagendr + dmdeduc2, 
        data = df, 
        expnms = names(df)[1:18], 
        emmvar = "ridreth1", 
        family = gaussian(), 
        q = 4, 
        alpha = 0.05
      )
    }, error = function(e) {
      NA
    })
    
    end.time <- Sys.time()
    list_times[[2]] <- end.time - start.time
    
    # save models
    list_mods <- list(mod1, mod2)
    
    # write out
    write_rds(list_mods, file = 
                paste0("mods/qgc_lg_", names(out2_resp1_re)[i], "_", i, 
                       "df.RDS"))
    write_rds(list_times, file = 
                paste0("times/qgc_lg_", names(out2_resp1_re)[i], "_", i, 
                       "df.RDS"))
    
  }
  
}

# run for all 
qgjob05 <- slurm_call(
  run_qgc_lg_re, params = list(vector = 1:100),
  global_objects = c('out2_resp1_re'),
  jobname = 'qlgre01',
  slurm_options = list(mem = '8G'))

qgjob06 <- slurm_call(
  run_qgc_lg_re, params = list(vector = 101:200),
  global_objects = c('out2_resp1_re'),
  jobname = 'qlgre02',
  slurm_options = list(mem = '8G'))

qgjob07 <- slurm_call(
  run_qgc_lg_re, params = list(vector = 201:300),
  global_objects = c('out2_resp1_re'),
  jobname = 'qlgre03',
  slurm_options = list(mem = '8G'))

qgjob08 <- slurm_call(
  run_qgc_lg_re, params = list(vector = 301:400),
  global_objects = c('out2_resp1_re'),
  jobname = 'qlgre04',
  slurm_options = list(mem = '8G'))

# check on status of job
system('squeue --me')

system('squeue --format="%.6i %.6P %.12j %.11M %.11l %.4D %.8R" --me')
