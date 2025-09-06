library(tidyverse)
library(qgcompint)

# read in data
comb <- read_rds("nhanes_data/processed_data.RDS")

# log-transform target data
# all POP's, serum cotinine, and telomere
comb_log <- comb |> 
  haven::zap_label() |> 
  mutate(across(c(LBX074LA:LBXF08LA, LBXCOT, TELOMEAN), log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric)) |> 
  #re-order variables
  select(SEQN:LBXF08LA, #ID and pop's
         BMXBMI:LBXBAPCT, #bio
         RIDAGEYR:DMDEDUC2, #demo, gender is binary, re and educ multi
         TELOMEAN) |> 
  #gender to 0's and 1's: 0 male, 1 female
  #should really be called sex
  mutate(RIAGENDR = RIAGENDR - 1) |> 
  mutate(across(RIAGENDR:DMDEDUC2, as.factor)) |> 
  as.data.frame() 
names(comb_log) <- tolower(names(comb_log))

mod <- qgcomp.emm.boot(
  f = telomean ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
    lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
    lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
    bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
    lbxbapct + ridageyr + riagendr + dmdeduc2, 
  data = comb_log, 
  expnms = names(comb_log)[2:19], 
  emmvar = "ridreth1", 
  family = gaussian(), 
  q = 4, 
  alpha = 0.05, 
  B = 100, 
  seed = 0
)
mod
plot(mod, emmval = 1)
plot(mod, emmval = 2)
plot(mod, emmval = 3)
plot(mod, emmval = 4)
plot(mod, emmval = 5)

# combine poc
comb_log2 <- comb_log |> 
  mutate(ridreth1 = ifelse(ridreth1 == 3, 1, 0), 
         ridreth1 = as.factor(ridreth1))
mod2 <- qgcomp.emm.boot(
  f = telomean ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
    lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
    lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
    bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
    lbxbapct + ridageyr + riagendr + dmdeduc2, 
  data = comb_log2, 
  expnms = names(comb_log)[2:19], 
  emmvar = "ridreth1", 
  family = gaussian(), 
  q = 4, 
  alpha = 0.05, 
  B = 100, 
  seed = 0
)
mod2
plot(mod2, emmval = 0)
plot(mod2, emmval = 1)

mod2.n <- qgcomp.emm.noboot(
  f = telomean ~ lbx074la + lbx099la + lbx138la + lbx153la + lbx170la + lbx180la +
    lbx187la + lbx194la + lbxpcbla + lbxhxcla + lbx118la + lbxd03la +
    lbxd05la + lbxd07la + lbxf03la + lbxf04la + lbxf05la + lbxf08la +
    bmxbmi + lbxcot + lbxwbcsi + lbxlypct + lbxmopct + lbxnepct + lbxeopct +
    lbxbapct + ridageyr + riagendr + dmdeduc2, 
  data = comb_log2, 
  expnms = names(comb_log)[2:19], 
  emmvar = "ridreth1", 
  family = gaussian(), 
  q = 4, 
  alpha = 0.05
)
mod2.n
plot(mod2.n, emmval = 0)
plot(mod2.n, emmval = 1)
