library(tidyverse)

# for POP's, keep list used by Gibson et al

# demo data
demo <- haven::read_xpt("nhanes_data/raw_data/DEMO_B.XPT") |> 
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC2)
  # riagendr: 1 male, 2 female
  # ridreth1: 1 mex amer, 2 oth hisp, 3 nonhisp white, 4 nonhisp black, 5 other
  # dmdeduc2: 1 <9th, 2 9-11, 3 hs, 4 aa, 5 college
# dioxins, furans, coplanar pcb's
dcp <- haven::read_xpt("nhanes_data/raw_data/L28POC_B.XPT") |> 
  select(SEQN, 
         # lipid adjusted non-dioxin pcb's
         LBX074LA, #pcb74
         LBX099LA, #pcb99
         LBX138LA, #pcb138
         LBX153LA, #pcb153
         LBX170LA, #pcb170
         LBX180LA, #pcb180
         LBX187LA, #pcb187
         LBX194LA, #pcb194
         # lipid adjusted non-orth pcb's
         LBXPCBLA, #pcb126 
         LBXHXCLA, #pcb194 
         # lipid adjusted mono ortho pcb
         LBX118LA, #pcb118
         # lipid adjusted dioxins
         LBXD03LA, #1,2,3,6,7,8-hxcdd
         LBXD05LA, #1,2,3,4,6,7,8-hpcdd
         LBXD07LA, #1,2,3,4,6,7,8,9-ocdd
         # lipid adjusted furans
         LBXF03LA, #2,3,4,7,8-pncdf
         LBXF04LA, #1,2,3,4,7,8-hxcdf
         LBXF05LA, #1,2,3,6,7,8-hxcdf
         LBXF08LA, #1,2,3,4,6,7,8-hxcdf
         )
# bmi
bmi <- haven::read_xpt("nhanes_data/raw_data/BMX_B.XPT") |> 
  select(SEQN, BMXBMI)
# cotinine
cot <- haven::read_xpt("nhanes_data/raw_data/L06_B.XPT") |> 
  select(SEQN, LBXCOT)
# blood count
blood <- haven::read_xpt("nhanes_data/raw_data/L25_B.XPT") |> 
  select(SEQN, 
         LBXWBCSI, #white blood cell count
         LBXLYPCT, #lymphocyte percent
         LBXMOPCT, #monocyte percent
         LBXNEPCT, #neutrophil percent
         LBXEOPCT, #eosinophil percent
         LBXBAPCT, #basophil percent
         )
# telomeres
telo <- haven::read_xpt("nhanes_data/raw_data/TELO_B.XPT") |> 
  select(SEQN, TELOMEAN)

# join together
# dcp has smallest sample size, so start there
joined <- dcp |> 
  inner_join(demo, by = "SEQN") |> 
  inner_join(bmi, by = "SEQN") |> 
  inner_join(cot, by = "SEQN") |> 
  inner_join(blood, by = "SEQN") |> 
  inner_join(telo, by = "SEQN")

# now, omit NA's
joined_small <- joined |> 
  na.omit()
# sample size of 1,003 here, identical to gibson et al

# re-code factor variables
joined_recode <- joined_small |> 
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")), 
    RIDRETH1 = factor(RIDRETH1, levels = 1:5, 
                      labels = c("Mex. Ameri.", "Oth. Hisp.", "Non-Hisp. White", 
                                 "Non-Hisp. Black", "Other")), 
    DMDEDUC2 = factor(
      ifelse(DMDEDUC2 == 1, DMDEDUC2, DMDEDUC2 - 1), levels = 1:4, 
      labels = c("Less than HS", "HS", "Some College", "College")
      )
  )
Hmisc::label(joined_recode$RIAGENDR) <- "Gender"
Hmisc::label(joined_recode$RIDRETH1) <- "Race/Ethnicity - Recode"
Hmisc::label(joined_recode$DMDEDUC2) <- "Education Level - Adults 20+"

# save to rds
write_rds(joined_small, "nhanes_data/processed_data.RDS")

# # non-dioxin pcb's
# pcb <- haven::read_xpt("nhanes_data/raw_data/SSPCB_B.XPT")
