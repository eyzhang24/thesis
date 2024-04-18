# this code splits up the large "sim_preds_lg.RDS" file 
# so that we can push it to git

library(readr)

## starting from the big file

# read in big file
sim_preds_lg <- read_rds("sim/sim_preds_lg.RDS")

# split up
sim_preds_lg1 <- sim_preds_lg[1:700]
sim_preds_lg2 <- sim_preds_lg[701:1400]
sim_preds_lg3 <- sim_preds_lg[1401:2100]

# write back to rds
write_rds(sim_preds_lg1, "sim/sim_preds_lg1.RDS")
write_rds(sim_preds_lg2, "sim/sim_preds_lg2.RDS")
write_rds(sim_preds_lg3, "sim/sim_preds_lg3.RDS")

## starting from separated files

# read in separated files
sim_preds_lg1 <- read_rds("sim/sim_preds_lg1.RDS")
sim_preds_lg2 <- read_rds("sim/sim_preds_lg2.RDS")
sim_preds_lg3 <- read_rds("sim/sim_preds_lg3.RDS")

# merge
sim_preds_lg <- c(sim_preds_lg1, sim_preds_lg2, sim_preds_lg3)

# write back to big rds for other use in code
write_rds(sim_preds_lg, "sim/sim_preds_lg.RDS")