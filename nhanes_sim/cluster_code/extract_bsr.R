library(NLinteraction)
library(tidyverse)

cnames <- c('LBX074LA', 'LBX099LA', 'LBX138LA', 'LBX153LA', 'LBX170LA', 
  'LBX180LA', 'LBX187LA', 'LBX194LA', 'LBXPCBLA', 'LBXHXCLA', 'LBX118LA', 
  'LBXD03LA', 'LBXD05LA', 'LBXD07LA', 'LBXF03LA', 'LBXF04LA', 'LBXF05LA', 
  'LBXF08LA')

######
# extract bsr smalls!
######

# extracting file names

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_ssm <- list_files[grepl("_rslurm_bsr_sm", list_files)]

ssm_subf <- list.dirs(list_ssm, full.names = TRUE, recursive = TRUE)
ssm_mod <- ssm_subf[grepl("mods", ssm_subf)]

ssm_labels <- gsub("\\D", "", ssm_mod)
ssm_labels <- ifelse(ssm_labels == "", 1, as.numeric(ssm_labels))

# get paths
ssm_paths <- ssm_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ssm_labels)

# extract PIPs
ssm_pips <- names(ssm_paths) |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- data.frame(variable = cnames, PIP = bsr$MainPIP) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns)
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_pips, "sim/bsr_sm/pips.csv")

# extract bivariate pip's
ssm_pip_biv <- names(ssm_paths) |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- reshape2::melt(bsr$InteractionPIP, 
                                 na.rm = TRUE, 
                                 value.name = "PIP") |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns)
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_pip_biv, "sim/bsr_sm/pip_biv.csv")

# extract trivariate pip's
ssm_pip_triv <- names(ssm_paths)[14:17] |> # trivariate only
  purrr::map_df(\(x) {
    print(paste0("starting at ", x)) 
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y) 
        result <- data.frame(
          PIP = InteractionProb(NLmod = bsr, Xsub = c(13, 9, 8)), 
          trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
          df = bsr$ns
        )
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_pip_triv, "sim/bsr_sm/pip_triv.csv")

# trivariate fdr's
combs <- combn(1:18, 3)
lcombs <- split(combs, col(combs))
fcombs <- lcombs[!sapply(lcombs, function(x) all(x == c(8, 9, 13)))]

ssm_pip_triv_fdr <- names(ssm_paths)[14:17] |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x)) 
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y) 
        yres <- fcombs |> 
          purrr::map_df(\(z) {
          zres <- data.frame(
            PIP = InteractionProb(NLmod = bsr, Xsub = z), 
            comb = paste(z, collapse = ", "),
            trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
            df = bsr$ns
          )
          return(zres)
        })
        rm(bsr)
        return(yres)
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_pip_triv_fdr, "sim/bsr_sm/pip_triv_fdr.csv")

# read back in data
out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")
source("extract_fxns.R")

# extract bivariate relationships
ssm_biv <- names(ssm_paths)[2:13] |> # bivariate only 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    indices <- case_when(
      x %in% 2:5 ~ c(13, 8), # D05 and 194
      x %in% 6:9 ~ c(18, 15), # F08 and F03
      x %in% 10:13 ~ c(1, 8) # 074 & 194
    )
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        trial <- as.numeric(sub(".+_(.+)d.+", "\\1", y))
        if(trial %% 5 == 0) print(paste0("index ", y))
        df <- out1_resp1[[trial]]
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
        
        bsr <- read_rds(y)
        
        result1 <- bivarsurf_bsr(bsr, X = X, C = C, 
                                 j1 = indices[1], j2 = indices[2], 
                                gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[indices[1]], 
                 j2 = cnames[indices[2]])
        result2 <- bivarsurf_bsr(bsr, X = X, C = C, j1 = indices[2], j2 = indices[1], 
                                 gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                 quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[indices[2]], 
                 j2 = cnames[indices[1]])
        
        rm(bsr)
        return(bind_rows(result1, result2))
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_biv, "sim/bsr_sm/biv_expresp.csv")

# extract trivariate relationships
ssm_triv <- names(ssm_paths)[14:17] |> # trivariate only 
  purrr::map_df(\(x) {
    message("starting ", x)
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        trial <- as.numeric(sub(".+_(.+)d.+", "\\1", y))
        if(trial %% 5 == 0) message("   index ", trial)
        df <- out1_resp1[[trial]]
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
        
        bsr <- read_rds(y)
        
        result1 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 13, j2 = 9, j3 = 8,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        result2 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 9, j2 = 8, j3 = 13,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        result3 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 8, j2 = 13, j3 = 9,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        
        rm(bsr)
        return(bind_rows(result1, result2, result3))
      }) |> 
      mutate(case = x)
  })
write_csv(ssm_triv, "sim/bsr_sm/triv_expresp.csv")

# can we get confidence intervals from bsr at all? 

# extract run times
ssm_subf <- list.dirs(list_ssm, full.names = TRUE, recursive = TRUE)
ssm_times <- ssm_subf[grepl("times", ssm_subf)]

ssm_labelt <- gsub("\\D", "", ssm_times)
ssm_labelt <- ifelse(ssm_labelt == "", 1, as.numeric(ssm_labelt))

# get paths
ssm_patht <- ssm_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ssm_labelt)

t <- read_rds(ssm_patht[[1]][1])

times <- names(ssm_patht) |> 
  map_df(\(x) {
    ssm_patht[[x]] |> 
      map_df(\(y) {
        time = read_rds(y)
        # print(time)
        return(data.frame(
          time_selection = time[[1]],
          time = time[[2]], 
          trial = as.numeric(sub(".+_(.+)d.+", "\\1", y))))
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/bsr_sm/times.RDS")

######
# extract bsr larges!
######

# extracting file names

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_slg <- list_files[grepl("_rslurm_bsr_lg", list_files)]

slg_subf <- list.dirs(list_slg, full.names = TRUE, recursive = TRUE)
slg_mod <- slg_subf[grepl("mods", slg_subf)]

slg_labels <- gsub("\\D", "", slg_mod)
slg_labels <- ifelse(slg_labels == "", 1, as.numeric(slg_labels))

# get paths
slg_paths <- slg_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = slg_labels)

# extract PIPs
slg_pips <- names(slg_paths) |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    slg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- data.frame(variable = cnames, PIP = bsr$MainPIP) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns)
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(slg_pips, "sim/bsr_lg/pips.csv")

# extract bivariate pip's
slg_pip_biv <- names(slg_paths) |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    slg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- reshape2::melt(bsr$InteractionPIP, 
                                 na.rm = TRUE, 
                                 value.name = "PIP") |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns)
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(slg_pip_biv, "sim/bsr_lg/pip_biv.csv")

# extract trivariate pip's
slg_pip_triv <- names(slg_paths)[14:17] |> # trivariate only
  purrr::map_df(\(x) {
    print(paste0("starting at ", x)) 
    slg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y) 
        result <- data.frame(
          PIP = InteractionProb(NLmod = bsr, Xsub = c(13, 9, 8)), 
          trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
          df = bsr$ns
        )
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(slg_pip_triv, "sim/bsr_lg/pip_triv.csv")

# trivariate fdr's

# read back in data
out2_resp1 <- read_rds("sim/sim_resp_lg_a.RDS")
source("extract_fxns.R")

# extract bivariate relationships
slg_biv <- names(slg_paths)[2:13] |> # bivariate only 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    indices <- case_when(
      x %in% 2:5 ~ c(13, 8), # D05 and 194
      x %in% 6:9 ~ c(18, 15), # F08 and F03
      x %in% 10:13 ~ c(1, 8) # 074 & 194
    )
    slg_paths[[x]] |> 
      purrr::map_df(\(y) {
        trial <- as.numeric(sub(".+_(.+)d.+", "\\1", y))
        if(trial %% 5 == 0) print(paste0("index ", y))
        df <- out2_resp1[[trial]]
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
        
        bsr <- read_rds(y)
        
        result1 <- bivarsurf_bsr(bsr, X = X, C = C, j1 = indices[1], j2 = indices[2], 
                                 gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                 quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[indices[1]], 
                 j2 = cnames[indices[2]])
        result2 <- bivarsurf_bsr(bsr, X = X, C = C, j1 = indices[2], j2 = indices[1], 
                                 gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                 quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[indices[2]], 
                 j2 = cnames[indices[1]])
        
        rm(bsr)
        return(bind_rows(result1, result2))
      }) |> 
      mutate(case = x)
  })
write_csv(slg_biv, "sim/bsr_lg/biv_expresp.csv")

# extract trivariate relationships
slg_triv <- names(slg_paths)[14:17] |> # trivariate only 
  purrr::map_df(\(x) {
    message("starting ", x)
    slg_paths[[x]] |> 
      purrr::map_df(\(y) {
        trial <- as.numeric(sub(".+_(.+)d.+", "\\1", y))
        if(trial %% 5 == 0) message("   index ", trial)
        df <- out2_resp1[[trial]]
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
        
        bsr <- read_rds(y)
        
        result1 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 13, j2 = 9, j3 = 8,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        result2 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 9, j2 = 8, j3 = 13,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        result3 <- trivarsurf_bsr(bsr, X = X, C = C, j1 = 8, j2 = 13, j3 = 9,
                                  gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                                  quantile_rest = 0.5) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
                 df = bsr$ns, 
                 j1 = cnames[13], 
                 j2 = cnames[9], 
                 j3 = cnames[8])
        
        rm(bsr)
        return(bind_rows(result1, result2, result3))
      }) |> 
      mutate(case = x)
  })
write_csv(slg_triv, "sim/bsr_lg/triv_expresp.csv")

# trivariate fdr's
combs <- combn(1:18, 3)
lcombs <- split(combs, col(combs))
fcombs <- lcombs[!sapply(lcombs, function(x) all(x == c(8, 9, 13)))]

slg_pip_triv_fdr <- names(slg_paths)[14:17] |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x)) 
    ssm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y) 
        yres <- fcombs |> 
          purrr::map_df(\(z) {
            zres <- data.frame(
              PIP = InteractionProb(NLmod = bsr, Xsub = z), 
              comb = paste(z, collapse = ", "),
              trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)), 
              df = bsr$ns
            )
            return(zres)
          })
        rm(bsr)
        return(yres)
      }) |> 
      mutate(case = x)
  })
write_csv(slg_pip_triv_fdr, "sim/bsr_lg/pip_triv_fdr.csv")

# extract run times
slg_subf <- list.dirs(list_slg, full.names = TRUE, recursive = TRUE)
slg_times <- slg_subf[grepl("times", slg_subf)]

slg_labelt <- gsub("\\D", "", slg_times)
slg_labelt <- ifelse(slg_labelt == "", 1, as.numeric(slg_labelt))

# get paths
slg_patht <- slg_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = slg_labelt)

t <- read_rds(slg_patht[[1]][1])

times <- names(slg_patht) |> 
  map_df(\(x) {
    slg_patht[[x]] |> 
      map_df(\(y) {
        time = read_rds(y)
        # print(time)
        return(data.frame(
          time_selection = time[[1]],
          time = time[[2]], 
          trial = as.numeric(sub(".+_(.+)d.+", "\\1", y))))
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/bsr_lg/times.RDS")

# extract bsr re's -------------------------------------------------------

source("extract_fxns.R")

# read back in data
out1_resp1_re <- read_rds("sim/sim_resp_sm_re.RDS")

# small

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_ssmre <- list_files[grepl("_rslurm_ssmre", list_files)]

ssmre_subf <- list.dirs(list_ssmre, full.names = TRUE, recursive = TRUE)
ssmre_mod <- ssmre_subf[grepl("mods", ssmre_subf)]

ssmre_labels <- gsub("\\D", "", ssmre_mod)
ssmre_labels <- ifelse(ssmre_labels == "", 1, as.numeric(ssmre_labels))

# get paths
ssmre_paths <- ssmre_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ssmre_labels)

one <- read_rds(ssmre_paths[[1]][[1]])
onex <- one[[1]]

# extract PIPs
ssmre_pips <- names(ssmre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    ssmre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        
        1:6 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod) == "logical") {
              result <- data.frame(variable = NA, PIP = NA, 
                                   df = NA, race = z)
            } else {
              result <- data.frame(variable = cnames, PIP = mod$MainPIP) |>
                mutate(df = mod$ns, race = z)
            }
            return(result)
          }) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(ssmre_pips, "sim/re/ssmre_pips.csv")
write.csv(ssmre_pips, "sim/re/ssmre_pips.csv")

# get exposure response relationships
ssmre_univ <- names(ssmre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    ssmre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        trial <- as.numeric(sub(".+_(.+)d.+", "\\1", y))
        if(trial %% 10 == 0) message("index ", y)
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
        
        list_mod <- read_rds(y)
        1:6 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(j1val = NA, est = NA, 
                               lower = NA, upper = NA, race = z)
            } else { 
              # compute estimated Hg-response relationship
              df <- univarsurf_bsr(mod, X, C, j1 = 13) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(ssmre_univ, "sim/re/ssmre_expresp.csv")  

# extract run times
ssmre_subf <- list.dirs(list_ssmre, full.names = TRUE, recursive = TRUE)
ssmre_times <- ssmre_subf[grepl("times", ssmre_subf)]

ssmre_labelt <- gsub("\\D", "", ssmre_times)
ssmre_labelt <- ifelse(ssmre_labelt == "", 1, as.numeric(ssmre_labelt))

# get paths
ssmre_patht <- ssmre_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ssmre_labelt)

t <- read_rds(ssmre_patht[[1]][[1]])

times <- names(ssmre_patht) |> 
  map_df(\(x) {
    ssmre_patht[[x]] |> 
      map_df(\(y) {
        list_time = read_rds(y)
        df <- 1:6 |> 
          map_df(\(z) {
            data.frame(time_selection = list_time[[z]][[1]], 
                       time = list_time[[z]][[2]], 
                       race = z)
          }) |> 
          mutate(trial = as.numeric(sub(".+_(.+)d.+", "\\1", y)))
        return(df)
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/re/ssmre_times.RDS")


# large

# read back in data
out2_resp1_re <- read_rds("sim/sim_resp_lg_re.RDS")

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_slgre <- list_files[grepl("_rslurm_slgre", list_files)]

slgre_subf <- list.dirs(list_slgre, full.names = TRUE, recursive = TRUE)
slgre_mod <- slgre_subf[grepl("mods", slgre_subf)]

slgre_labels <- gsub("\\D", "", slgre_mod)
slgre_labels <- ifelse(slgre_labels == "", 1, as.numeric(slgre_labels))

# get paths
slgre_paths <- slgre_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = slgre_labels)

one <- read_rds(slgre_paths[[1]][[1]])
oney <- SingVarRiskSummaries(one[[2]], q.fixed = 0.5)
onex <- SingVarRiskSummaries(one[[2]], which.z = 4, q.fixed = 0.5)
onez <- PredictorResponseUnivar(one[[2]], which.z = 4)

# extract PIPs
extract_sre_pip <- function(case) {
  output <- slgre_paths[[case]] |> 
    purrr::map_df(\(y) { # for each trial
      list_mod <- read_rds(y)
      1:3 |> # for each stratified model
        purrr::map_df(\(z) {
          mod <- list_mod[[z]]
          if(class(mod) == "logical") {
            result <- data.frame(variable = NA, PIP = NA, 
                                 df = NA, race = z)
          } else {
            result <- data.frame(variable = cnames, PIP = mod$MainPIP) |>
              mutate(df = mod$ns, race = z)
          }
          return(result)
        }) |> 
        mutate(trial = as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y)))
    }) |> 
    mutate(case = case)
}
slgre_pips1 <- extract_sre_pip(case = 1)
slgre_pips2 <- extract_sre_pip(case = 2)
slgre_pips3 <- extract_sre_pip(case = 3)
slgre_pips4 <- extract_sre_pip(case = 4)
slgre_pips <- bind_rows(slgre_pips1, slgre_pips2, slgre_pips3, slgre_pips4)

write_csv(slgre_pips, "sim/re/slgre_pips.csv")
# write.csv(slgre_pips, "sim/re/slgre_pips.csv")

# get exposure response relationships
extract_slgre_univ <- function(x) {
  output <- slgre_paths[[x]] |> 
    purrr::map_df(\(y) { # for each trial
      trial <- as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y))
      if(trial %% 10 == 0) message("index ", y)
      df <- out2_resp1_re[[trial]] |> 
        mutate(RIDRETH1 = case_when(
          RIDRETH1 %in% c(1, 2) ~ 1, 
          RIDRETH1 %in% c(3) ~ 2, 
          RIDRETH1 %in% c(4, 5) ~ 3, 
          .default = NA
        ))
      X <- df |> 
        select(LBX074LA:LBXF08LA) |> 
        as.matrix.data.frame()
      C <- df |> 
        mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
               DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
               DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
               DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
        select(-DMDEDUC2) |> 
        select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT) |> 
        as.matrix.data.frame()
      
      list_mod <- read_rds(y)
      1:3 |> # for each stratified model
        purrr::map_df(\(z) {
          mod <- list_mod[[z]]
          Xt <- X[df$RIDRETH1 == z, ]
          Ct <- C[df$RIDRETH1 == z, ]
          if(class(mod)[1] == "logical") { # if output failed
            df <- data.frame(j1val = NA, est = NA, 
                             lower = NA, upper = NA, race = z)
          } else { 
            # compute estimated Hg-response relationship
            df <- univarsurf_bsr(mod, Xt, Ct, j1 = 13) |> 
              mutate(race = z)
          }
          # print(df)
          return(df)
        }) |> 
        mutate(trial = as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y)))
    }) |> 
    mutate(case = x)
}
slgre_univ <- names(slgre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    slgre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        trial <- as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y))
        if(trial %% 10 == 0) message("index ", y)
        df <- out2_resp1_re[[trial]] |> 
          mutate(RIDRETH1 = case_when(
            RIDRETH1 %in% c(1, 2) ~ 1, 
            RIDRETH1 %in% c(3) ~ 2, 
            RIDRETH1 %in% c(4, 5) ~ 3, 
            .default = NA
          ))
        X <- df |> 
          select(LBX074LA:LBXF08LA) |> 
          as.matrix.data.frame()
        C <- df |> 
          mutate(DMDEDUC22 = ifelse(DMDEDUC2 == 2, 1, 0), 
                 DMDEDUC23 = ifelse(DMDEDUC2 == 3, 1, 0), 
                 DMDEDUC24 = ifelse(DMDEDUC2 == 4, 1, 0), 
                 DMDEDUC25 = ifelse(DMDEDUC2 == 5, 1, 0)) |> 
          select(-DMDEDUC2) |> 
          select(RIDAGEYR, RIAGENDR, starts_with("DMDEDUC2"), BMXBMI:LBXBAPCT) |> 
          as.matrix.data.frame()
        
        list_mod <- read_rds(y)
        1:3 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            Xt <- X[df$RIDRETH1 == z, ]
            Ct <- C[df$RIDRETH1 == z, ]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(j1val = NA, est = NA, 
                               lower = NA, upper = NA, race = z)
            } else { 
              # compute estimated Hg-response relationship
              df <- univarsurf_bsr(mod, Xt, Ct, j1 = 13) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
slre_univ1 <- extract_slgre_univ(x = 1)
slre_univ2 <- extract_slgre_univ(x = 2)
slre_univ3 <- extract_slgre_univ(x = 3)
slre_univ4 <- extract_slgre_univ(x = 4)
slgre_univ <- bind_rows(slre_univ1, slre_univ2, slre_univ3, slre_univ4)
write_csv(slgre_univ, "sim/re/slgre_expresp.csv")  

# extract run times
slgre_subf <- list.dirs(list_slgre, full.names = TRUE, recursive = TRUE)
slgre_times <- slgre_subf[grepl("times", slgre_subf)]

slgre_labelt <- gsub("\\D", "", slgre_times)
slgre_labelt <- ifelse(slgre_labelt == "", 1, as.numeric(slgre_labelt))

# get paths
slgre_patht <- slgre_times |> 
  purrr::map(\(x) {
    temp <- list.files(x, full.names = TRUE)
    temp <- temp[which(!grepl("df", temp))]
    return(temp)
  }) |> 
  setNames(nm = slgre_labelt)

# slgre_patht <- slgre_patht0[which(!grepl("df", slgre_patht0))]

t <- read_rds(slgre_patht[[1]][1])

times <- names(slgre_patht) |> 
  map_df(\(x) {
    slgre_patht[[x]] |> 
      map_df(\(y) {
        list_time = read_rds(y)
        df <- 1:3 |> 
          map_df(\(z) {
            data.frame(time_selection = list_time[[z]][[1]], 
                       time = list_time[[z]][[2]], 
                       race = z)
          }) |> 
          mutate(trial = as.numeric(sub(".*_([0-9]+)\\.RDS", "\\1", y)))
        return(df)
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/re/slgre_times.RDS")
