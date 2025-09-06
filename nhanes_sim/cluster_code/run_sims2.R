library(tidyverse)
library(rslurm)
library(bkmr)
library(NLinteraction)

#########
# run mlr and oracle ------------------------------------------------------
#########

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")
run_mlr_sm <- function() {
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
  mlrs <- vector(mode='list', length = 1700)
  names(mlrs) <- names(out2_resp1)
  mlrtimes <- vector(mode = 'list', length = 1700)
  names(mlrtimes) <- names(out2_resp1)
  
  chems_only <- vector(mode = 'list', length = 1700)
  names(chems_only) <- names(out2_resp1)
  
  oracles <- vector(mode='list', length = 1700)
  names(oracles) <- names(out2_resp1)
  oracletimes <- vector(mode = 'list', length = 1700)
  names(oracletimes) <- names(out2_resp1)
  
  
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
  
  return(list(mlrs, mlrtimes, oracles, oracletimes, chems_only))
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

############
# bkmr --------------------------------------------------------------------
###########

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")

run_bkmr_sm <- function(vector) {
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
    df <- out1_resp1[[i]]
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
    y <- df$y 
    
    # fit model and save time
    set.seed(0)
    start.time <- Sys.time()
    mod <- kmbayes(y = y, Z = Z, X = X, 
                          iter = 50000, verbose = FALSE, varsel = TRUE)
    end.time <- Sys.time()
    bkmr_times[[i]] <- end.time - start.time
    
    # save model and remove from memory
    write_rds(mod, file = paste0("mods/bkmr_sm_", names(out1_resp1)[i], "_", i, ".RDS"))
    write_rds(bkmr_times[[i]], file = 
                paste0("times/bkmr_sm_", names(out1_resp1)[i], "_", i, ".RDS"))
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
  if(!dir.exists("times")) {
    dir.create("times")
  }
  
  # extract names of files
  list_files <- list.files("mods", full.names = TRUE)
  nums <- as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", list_files))
  print(nums)
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
  
  print(starting)
  print(max(vector))
  for(i in starting:max(vector)) {
    print(paste0("------------- run ", i, "--------------"))
    
    # prepare data
    df <- out2_resp1[[i]]
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
    y <- df$y 
    knots <- fields::cover.design(Z, nd = 100)$design
    
    # fit model and save time
    set.seed(0)
    start.time <- Sys.time()
    mod <- kmbayes(y = y, Z = Z, X = X, knots = knots,
                   iter = 50000, verbose = FALSE, varsel = TRUE)
    end.time <- Sys.time()
    bkmr_times[[i]] <- end.time - start.time
    
    # save model and remove from memory
    write_rds(mod, file = paste0("mods/bkmr_lg_", names(out2_resp1)[i], "_", i, ".RDS"))
    write_rds(bkmr_times[[i]], file = 
                paste0("times/bkmr_lg_", names(out2_resp1)[i], "_", i, ".RDS"))
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


###############
# bsr  --------------------------------------------------------------------
###############

### smaller sample size

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")

run_bsr_sm <- function(vector) {
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
    df <- out1_resp1[[i]]
    X <- df |> 
      select(LBX074LA:LBXF08LA) |> 
      as.matrix.data.frame()
    C <- df |>
      bind_cols(
        data.frame(model.matrix(
          ~ RIDRETH1 - 1, data = mutate(df, RIDRETH1 = as.factor(RIDRETH1)))), 
        data.frame(model.matrix(
          ~ DMDEDUC2 - 1, data = mutate(df, DMDEDUC2 = as.factor(DMDEDUC2))))
      ) |> 
      select(RIDAGEYR, RIAGENDR, RIDRETH12:RIDRETH15, DMDEDUC22:DMDEDUC25, 
             BMXBMI:LBXBAPCT) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # fit model for d = {1, 2, 3, 4, 5} and save time
    list_times <- vector(mode = "list", length = 2)
    
    set.seed(0)
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
    write_rds(list_times, file = 
                 paste0("times/bsr_sm_", names(out1_resp1)[i], "_", i, 
                        "df", ind, ".RDS"))
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  write_rds(bsr_times, file = "bsr_smf_times.RDS")
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
    df <- out2_resp1[[i]]
    X <- df |> 
      select(LBX074LA:LBXF08LA) |> 
      as.matrix.data.frame()
    C <- df |>
      bind_cols(
        data.frame(model.matrix(
          ~ RIDRETH1 - 1, data = mutate(df, RIDRETH1 = as.factor(RIDRETH1)))), 
        data.frame(model.matrix(
          ~ DMDEDUC2 - 1, data = mutate(df, DMDEDUC2 = as.factor(DMDEDUC2))))
      ) |> 
      select(RIDAGEYR, RIAGENDR, RIDRETH12:RIDRETH15, DMDEDUC22:DMDEDUC25, 
             BMXBMI:LBXBAPCT) |> 
      as.matrix.data.frame()
    Y <- df$y
    
    # fit model for d = {1, 2, 3, 4, 5} and save time
    list_times <- vector(mode = "list", length = 2)
    
    set.seed(0)
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
                paste0("mods/bsr_lgf_", names(out2_resp1)[i], "_", i, 
                       "df", ind, ".RDS"))
    write_rds(list_times, file = 
                paste0("times/bsr_lgf_", names(out2_resp1)[i], "_", i, 
                       "df", ind, ".RDS"))
    
    rm(mod1, mod2, mod3, mod4, mod)
  }
  
  write_rds(bsr_times, file = "bsr_lgf_times.RDS")
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
  slurm_options = list(mem = '8G')) # 20 day time limit

# check on status of job
system('squeue --me')

