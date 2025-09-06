library(tidyverse)

# extract p-values from naive and oracle

# naive small
mlr_sm <- read_rds("sim/mlr_mods_sm.RDS")

one <- mlr_sm[[1]]

mlrsm_pval <- 1:1700 |> 
  map_df(\(x) {
    mod <- mlr_sm[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(mlrsm_pval, "sim/_mlr/pvalsm.csv")

# naive large
mlr_lg <- read_rds("sim/mlr_mods_lg.RDS")

mlrlg_pval <- 1:1700 |> 
  map_df(\(x) {
    mod <- mlr_lg[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(mlrlg_pval, "sim/_mlr/pvallg.csv")

# oracle small
orac_sm <- read_rds("sim/oracle_mods_sm.RDS")

oracsm_pval <- 1:1700 |> 
  map_df(\(x) {
    mod <- orac_sm[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(oracsm_pval, "sim/_oracle/pvalsm.csv")

# oracle large
orac_lg <- read_rds("sim/oracle_mods_lg.RDS")

oraclg_pval <- 1:1700 |> 
  map_df(\(x) {
    mod <- orac_lg[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(oraclg_pval, "sim/_oracle/pvallg.csv")

# naive small re
mlr_sm_re <- read_rds("sim/mlr_mods_sm_re.RDS")

mlrsmre_pval <- 1:400 |> 
  map_df(\(x) {
    mod <- mlr_sm_re[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(mlrsmre_pval, "sim/_mlr/pvalsmre.csv")

# naive large re
mlr_lg_re <- read_rds("sim/mlr_mods_lg_re.RDS")

mlrlgre_pval <- 1:400 |> 
  map_df(\(x) {
    mod <- mlr_lg_re[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(mlrlgre_pval, "sim/_mlr/pvallgre.csv")

# oracle small re
orac_sm_re <- read_rds("sim/oracle_mods_sm_re.RDS")

oracsmre_pval <- 1:400 |> 
  map_df(\(x) {
    mod <- orac_sm_re[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(oracsmre_pval, "sim/_oracle/pvalsmre.csv")

# oracle large re
orac_lg_re <- read_rds("sim/oracle_mods_lg_re.RDS")

oraclgre_pval <- 1:400 |> 
  map_df(\(x) {
    mod <- orac_lg_re[[x]]
    data.frame(summary(mod)$coefficients) |> 
      rownames_to_column(var = "var") |> 
      select(var, p = 5) |> 
      mutate(case = x)
  })
write_csv(oraclgre_pval, "sim/_oracle/pvallgre.csv")
