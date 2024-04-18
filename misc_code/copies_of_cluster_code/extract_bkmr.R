library(tidyverse)
library(bkmr)

######
# extract bkmr smalls!
######

# extracting file names
# want all _rslurm_bkmr_sm folders
# go into mods, extract all RDS's 
# create list of lists

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_ksm <- list_files[grepl("_rslurm_bkmr_sm", list_files)]

ksm_subf <- list.dirs(list_ksm, full.names = TRUE, recursive = TRUE)
ksm_mod <- ksm_subf[grepl("mods", ksm_subf)]

ksm_labels <- gsub("\\D", "", ksm_mod)
ksm_labels <- ifelse(ksm_labels == "", 1, as.numeric(ksm_labels))

# get paths
ksm_paths <- ksm_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ksm_labels)

# extract PIP's
ksm_pips <- names(ksm_paths) |> 
  purrr::map_df(\(x) {
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- data.frame(
          ExtractPIPs(bkmr)
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_pips, "sim/bkmr_sm/pips.csv")

# extract univariate relationships
ksm_univ <- names(ksm_paths) |> 
  purrr::map_df(\(x) {
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- data.frame(
          PredictorResponseUnivar(bkmr)
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_univ, "sim/bkmr_sm/univ_expresp.csv")

# extract bivariate relationships
ksm_biv <- names(ksm_paths)[2:13] |> # only for 2-way interaction
  purrr::map_df(\(x) {
    indices <- case_when(
      x %in% 2:5 ~ c(4, 5), # Hg and Ni
      x %in% 6:9 ~ c(1, 2), # Cd and As
      x %in% 10:13 ~ c(3, 5) # Co and Ni
    )
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        bivar <- PredictorResponseBivar(bkmr, 
                               z.pairs = rbind(indices), verbose = FALSE)
        result <- data.frame(
          PredictorResponseBivarLevels(
            pred.resp.df = bivar, 
            Z = bkmr$Z, qs = c(0.1, 0.5, 0.9))
          
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
# write_csv(ksm_biv, "sim/bkmr_sm/biv_expresp.csv")
# write_rds(ksm_biv, "sim/bkmr_sm/biv_expresp.rds")
ksm_biv_nona <- na.omit(ksm_biv)
write_csv(ksm_biv_nona, "sim/bkmr_sm/biv_expresp.csv")

# load in fxn
source("extract_fxns.R")

# extract trivariate relationships
ksm_triv <- names(ksm_paths)[14:17] |> # only for 3-way
  purrr::map_df(\(x) {
    message("starting ", x)
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- trivarsurf_bkmr(bkmr, 4, 5, 6, 
                                  qs.diff = c(0.1, 0.5, 0.9), 
                                  q.fixed = 0.5, ngrid = 50) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_triv, "sim/bkmr_sm/triv_expresp.csv")

# extract one vs. rest bivariate interactions
ksm_ints <- names(ksm_paths)[2:17] |> 
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    indices <- case_when(
      x %in% 2:5 ~ list(c(4, 5)), # Hg and Ni
      x %in% 6:9 ~ list(c(1, 2)), # Cd and As
      x %in% 10:13 ~ list(c(3, 5)), # Co and Ni
      x %in% 14:17 ~ list(c(4, 5, 6)) # Hg, Ni, Tl
    )
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- SingVarIntSummaries(bkmr, 
                                    which.z = indices[[1]], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    method = "approx")
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_ints, "sim/bkmr_sm/ints.csv")

# load in fxn
source("extract_fxns.R")

# extract one vs. other bivariate interactions
ksm_intb <- names(ksm_paths)[2:13] |> # only for two-way
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    indices <- case_when(
      x %in% 2:5 ~ list(c(4, 5)), # Hg and Ni
      x %in% 6:9 ~ list(c(1, 2)), # Cd and As
      x %in% 10:13 ~ list(c(3, 5)) # Co and Ni
    )
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- bivarinter_bkmr(bkmr, 
                                z1 = indices[[1]][1], 
                                z2 = indices[[1]][2], 
                                qs.diff = c(0.25, 0.75), 
                                qs.fixed = c(0.25, 0.75),
                                q.rest = 0.5)
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_intb, "sim/bkmr_sm/int_bivar.csv")

# FDR for Hg-Ni
indices <- t(combn(1:10, 2))
ksm_intb_hgni <- names(ksm_paths)[2:5] |> # only for Hg-Ni
  purrr::map_df(\(x) {
    message("starting ", x)
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- 1:45 |> 
          purrr::map_df(\(z) {
            ints <- bivarinter_bkmr(bkmr, 
                                    z1 = indices[z, ][1], 
                                    z2 = indices[z, ][2], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    q.rest = 0.5)
            data.frame(ints)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_intb_hgni, "sim/bkmr_sm/int_bivar_fullhgni.csv")

# FDR for rest
indices <- t(combn(1:10, 2))
ksm_intb_rest <- names(ksm_paths)[6:13] |> # only for Hg-Ni
  purrr::map_df(\(x) {
    message("starting ", x)
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- 1:45 |> 
          purrr::map_df(\(z) {
            ints <- bivarinter_bkmr(bkmr, 
                                    z1 = indices[z, ][1], 
                                    z2 = indices[z, ][2], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    q.rest = 0.5)
            data.frame(ints)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_intb_rest, "sim/bkmr_sm/int_bivar_fullrest.csv")

# extract one vs. 2 others trivariate interactions
ksm_intt <- names(ksm_paths)[14:17] |> # only for three-way
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    ksm_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- trivarinter_bkmr(bkmr, 
                                 z1 = 4, z2 = 5, z3 = 6,
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 q.rest = 0.5)
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(ksm_intt, "sim/bkmr_sm/int_trivar.csv")

# extract run times
ksm_subf <- list.dirs(list_ksm, full.names = TRUE, recursive = TRUE)
ksm_times <- ksm_subf[grepl("times", ksm_subf)]

ksm_labelt <- gsub("\\D", "", ksm_times)
ksm_labelt <- ifelse(ksm_labelt == "", 1, as.numeric(ksm_labelt))

# get paths
ksm_patht <- ksm_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ksm_labelt)

t <- read_rds(ksm_patht[[1]][1])

times <- names(ksm_patht) |> 
  map_df(\(x) {
    ksm_patht[[x]] |> 
      map_df(\(y) {
        # time = read_rds(y)
        # print(time)
        return(data.frame(
          time = read_rds(y),
          trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y))))
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/bkmr_sm/times.RDS")

######
# extract bkmr larges!
######
  
list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_klg <- list_files[grepl("_rslurm_bkmr_lg", list_files)]

klg_subf <- list.dirs(list_klg, full.names = TRUE, recursive = TRUE)
klg_mod <- klg_subf[grepl("mods", klg_subf)]

klg_labels <- gsub("\\D", "", klg_mod)
klg_labels <- ifelse(klg_labels == "", 1, as.numeric(klg_labels))

# get paths
klg_paths <- klg_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = klg_labels)

# extract PIP's
klg_pips <- names(klg_paths) |> 
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- data.frame(
          ExtractPIPs(bkmr)
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_pips, "sim/bkmr_lg/pips.csv")

# extract univariate relationships
klg_univ <- names(klg_paths) |> 
  purrr::map_df(\(x) {
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- data.frame(
          PredictorResponseUnivar(bkmr)
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_univ, "sim/bkmr_lg/univ_expresp.csv")

# extract bivariate relationships
klg_biv <- names(klg_paths)[2:13] |> # only for 2-way interaction
  purrr::map_df(\(x) {
    print(paste0("starting at ", x))
    indices <- case_when(
      x %in% 2:5 ~ c(4, 5), # Hg and Ni
      x %in% 6:9 ~ c(1, 2), # Cd and As
      x %in% 10:13 ~ c(3, 5) # Co and Ni
    )
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        bivar <- PredictorResponseBivar(bkmr, 
                                        z.pairs = rbind(indices), verbose = FALSE)
        result <- data.frame(
          PredictorResponseBivarLevels(
            pred.resp.df = bivar, 
            Z = bkmr$Z, qs = c(0.1, 0.5, 0.9))
          
        ) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
# write_csv(klg_biv, "sim/bkmr_lg/biv_expresp.csv")
# write_rds(klg_biv, "sim/bkmr_lg/biv_expresp.rds")
klg_biv_nona <- na.omit(klg_biv)
write_csv(klg_biv_nona, "sim/bkmr_lg/biv_expresp.csv")

# load in fxn
source("extract_fxns.R")

# extract trivariate relationships
klg_triv <- names(klg_paths)[14:17] |> # only for 3-way
  purrr::map_df(\(x) {
    message("starting ", x)
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- trivarsurf_bkmr(bkmr, 4, 5, 6, 
                                  qs.diff = c(0.1, 0.5, 0.9), 
                                  q.fixed = 0.5, ngrid = 50) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_triv, "sim/bkmr_lg/triv_expresp.csv")

# extract one vs. rest bivariate interactions
klg_ints <- names(klg_paths)[2:17] |> 
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    indices <- case_when(
      x %in% 2:5 ~ list(c(4, 5)), # Hg and Ni
      x %in% 6:9 ~ list(c(1, 2)), # Cd and As
      x %in% 10:13 ~ list(c(3, 5)), # Co and Ni
      x %in% 14:17 ~ list(c(4, 5, 6)) # Hg, Ni, Tl
    )
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- SingVarIntSummaries(bkmr, 
                                    which.z = indices[[1]], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    method = "approx")
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_ints, "sim/bkmr_lg/ints.csv")

# load in fxn
source("extract_fxns.R")

# extract one vs. other bivariate interactions
klg_intb <- names(klg_paths)[2:13] |> # only for two-way
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    indices <- case_when(
      x %in% 2:5 ~ list(c(4, 5)), # Hg and Ni
      x %in% 6:9 ~ list(c(1, 2)), # Cd and As
      x %in% 10:13 ~ list(c(3, 5)) # Co and Ni
    )
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- bivarinter_bkmr(bkmr, 
                                z1 = indices[[1]][1], 
                                z2 = indices[[1]][2], 
                                qs.diff = c(0.25, 0.75), 
                                qs.fixed = c(0.25, 0.75),
                                q.rest = 0.5)
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_intb, "sim/bkmr_lg/int_bivar.csv")

# FDR for Hg-Ni
tictoc::tic()
indices <- t(combn(1:10, 2))
klg_intb_hgni <- names(klg_paths)[2:5] |> # only for Hg-Ni
  purrr::map_df(\(x) {
    message("starting ", x)
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- 1:45 |> 
          purrr::map_df(\(z) {
            ints <- bivarinter_bkmr(bkmr, 
                                    z1 = indices[z, ][1], 
                                    z2 = indices[z, ][2], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    q.rest = 0.5)
            data.frame(ints)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_intb_hgni, "sim/bkmr_lg/int_bivar_fullhgni.csv")
tictoc::toc()

# FDR for rest
tictoc::tic()
indices <- t(combn(1:10, 2))
klg_intb_rest <- names(klg_paths)[6:13] |> # only for Hg-Ni
  purrr::map_df(\(x) {
    message("starting ", x)
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        result <- 1:45 |> 
          purrr::map_df(\(z) {
            ints <- bivarinter_bkmr(bkmr, 
                                    z1 = indices[z, ][1], 
                                    z2 = indices[z, ][2], 
                                    qs.diff = c(0.25, 0.75), 
                                    qs.fixed = c(0.25, 0.75),
                                    q.rest = 0.5)
            data.frame(ints)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_intb_rest, "sim/bkmr_lg/int_bivar_fullrest.csv")
tictoc::toc()

# extract one vs. 2 others trivariate interactions
klg_intt <- names(klg_paths)[14:17] |> # only for three-way
  purrr::map_df(\(x) {
    print(paste0("starting ", x))
    klg_paths[[x]] |> 
      purrr::map_df(\(y) {
        bkmr <- read_rds(y)
        ints <- trivarinter_bkmr(bkmr, 
                                 z1 = 4, z2 = 5, z3 = 6,
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 q.rest = 0.5)
        result <- data.frame(ints) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
        rm(bkmr)
        return(result)
      }) |> 
      mutate(case = x)
  })
write_csv(klg_intt, "sim/bkmr_lg/int_trivar.csv")


# extract bkmr re's -------------------------------------------------------

# small

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_ksmre <- list_files[grepl("_rslurm_ksmre", list_files)]

ksmre_subf <- list.dirs(list_ksmre, full.names = TRUE, recursive = TRUE)
ksmre_mod <- ksmre_subf[grepl("mods", ksmre_subf)]

ksmre_labels <- gsub("\\D", "", ksmre_mod)
ksmre_labels <- ifelse(ksmre_labels == "", 1, as.numeric(ksmre_labels))

# get paths
ksmre_paths <- ksmre_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ksmre_labels)

one <- read_rds(ksmre_paths[[1]][[1]])
oney <- SingVarRiskSummaries(one[[2]], q.fixed = 0.5)
onex <- SingVarRiskSummaries(one[[2]], which.z = 4, q.fixed = 0.5)
onez <- PredictorResponseUnivar(one[[2]], which.z = 4)

# extract PIPs
ksmre_pips <- names(ksmre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    ksmre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:6 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(variable = NA, PIP = NA, race = z)
            } else { 
              # get pips
              df <- ExtractPIPs(mod) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(ksmre_pips, "sim/re/ksmre_pips.csv")

# get confidence intervals
ksmre_ints <- names(ksmre_paths) |> 
  purrr::map_df(\(x) { # for each case
    ksmre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:6 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(est = NA, sd = NA, race = z)
            } else { 
              # compute estimate of 0.75 vs. 0.25 quantiles
              df <- SingVarRiskSummaries(mod, which.z = 4, q.fixed = 0.5) |> 
                select(est, sd) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(ksmre_ints, "sim/re/ksmre_ints.csv")

# get exposure response relationships
ksmre_univ <- names(ksmre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    ksmre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:6 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(z1 = NA, est = NA, se = NA, race = z)
            } else { 
              # compute estimated Hg-response relationship
              df <- PredictorResponseUnivar(mod, which.z = 4) |> 
                select(z1 = z, est, se) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(ksmre_univ, "sim/re/ksmre_expresp.csv")  

# extract run times
ksmre_subf <- list.dirs(list_ksmre, full.names = TRUE, recursive = TRUE)
ksmre_times <- ksmre_subf[grepl("times", ksmre_subf)]

ksmre_labelt <- gsub("\\D", "", ksmre_times)
ksmre_labelt <- ifelse(ksmre_labelt == "", 1, as.numeric(ksmre_labelt))

# get paths
ksmre_patht <- ksmre_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = ksmre_labelt)

t <- read_rds(ksmre_patht[[1]][1])

times <- names(ksmre_patht) |> 
  map_df(\(x) {
    ksmre_patht[[x]] |> 
      map_df(\(y) {
        # time = read_rds(y)
        # print(time)
        list_time <- read_rds(y)
        df <- 1:6 |> 
          map_df(\(z) {
            data.frame(time = list_time[[z]], 
                       race = z)
          })
        return(mutate(df, trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y))))
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/re/ksmre_times.RDS")

# large

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_klgre <- list_files[grepl("_rslurm_klgre", list_files)]

klgre_subf <- list.dirs(list_klgre, full.names = TRUE, recursive = TRUE)
klgre_mod <- klgre_subf[grepl("mods", klgre_subf)]

klgre_labels <- gsub("\\D", "", klgre_mod)
klgre_labels <- ifelse(klgre_labels == "", 1, as.numeric(klgre_labels))

# get paths
klgre_paths <- klgre_mod |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = klgre_labels)

one <- read_rds(klgre_paths[[1]][[1]])
oney <- SingVarRiskSummaries(one[[2]], q.fixed = 0.5)
onex <- SingVarRiskSummaries(one[[2]], which.z = 4, q.fixed = 0.5)
onez <- PredictorResponseUnivar(one[[2]], which.z = 4)

# get pips
klgre_pips <- names(klgre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    klgre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:5 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(variable = NA, PIP = NA, race = z)
            } else { 
              # get pips
              df <- ExtractPIPs(mod) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(klgre_pips, "sim/re/klgre_pips.csv")

# get confidence intervals
klgre_ints <- names(klgre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    klgre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:5 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(est = NA, sd = NA, race = z)
            } else { 
              # compute estimate of 0.75 vs. 0.25 quantiles
              df <- SingVarRiskSummaries(mod, which.z = 4, q.fixed = 0.5) |> 
                select(est, sd) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(klgre_ints, "sim/re/klgre_ints.csv")

# get exposure response relationships
klgre_univ <- names(klgre_paths) |> 
  purrr::map_df(\(x) { # for each case
    message("starting ", x)
    klgre_paths[[x]] |> 
      purrr::map_df(\(y) { # for each trial
        list_mod <- read_rds(y)
        1:5 |> # for each stratified model
          purrr::map_df(\(z) {
            mod <- list_mod[[z]]
            if(class(mod)[1] == "logical") { # if output failed
              df <- data.frame(z1 = NA, est = NA, se = NA, race = z)
            } else { 
              # compute estimated Hg-response relationship
              df <- PredictorResponseUnivar(mod, which.z = 4) |> 
                select(z1 = z, est, se) |> 
                mutate(race = z)
            }
            # print(df)
            return(df)
          }) |> 
          mutate(trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y)))
      }) |> 
      mutate(case = x)
  })
write_csv(klgre_univ, "sim/re/klgre_expresp.csv")  

# extract run times
klgre_subf <- list.dirs(list_klgre, full.names = TRUE, recursive = TRUE)
klgre_times <- klgre_subf[grepl("times", klgre_subf)]

klgre_labelt <- gsub("\\D", "", klgre_times)
klgre_labelt <- ifelse(klgre_labelt == "", 1, as.numeric(klgre_labelt))

# get paths
klgre_patht <- klgre_times |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = klgre_labelt)

t <- read_rds(klgre_patht[[1]][1])

times <- names(klgre_patht) |> 
  map_df(\(x) {
    klgre_patht[[x]] |> 
      map_df(\(y) {
        # time = read_rds(y)
        # print(time)
        list_time <- read_rds(y)
        df <- 1:5 |> 
          map_df(\(z) {
            data.frame(time = list_time[[z]], 
                       race = z)
          })
        return(mutate(df, trial = as.numeric(sub(".*_(\\d+)\\.RDS", "\\1", y))))
      }) |> 
      mutate(case = x)
  })

write_rds(times, "sim/re/klgre_times.RDS")
