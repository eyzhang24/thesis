# create functions for various response variables

# base case, no interactions
base_case <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

am1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.475*LBXD05LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

am2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.95*LBXD05LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

ap1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.24*LBXD05LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

ap2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.48*LBXD05LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

# univariately insignificant
bm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA +
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.5*LBXF08LA*LBXF03LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           
           rnorm(nrow(df), 0, 7))
}

bm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           1*LBXF08LA*LBXF03LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

bp1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.17*LBXF08LA*((LBXF03LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

bp2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.34*LBXF08LA*((LBXF03LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

cm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.6*LBX074LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

cm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           1.2*LBX074LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

cp1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBX074LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

cp2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.5*LBX074LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

dm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.36*LBXD05LA*LBXPCBLA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

dm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.72*LBXD05LA*LBXPCBLA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

dp1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.125*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

dp2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 7))
}

# interxn in smaller group, smaller effect size
em1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1.5, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5 + 0.5*LBXD05LA, # 1.5x in group 4
                     RIDRETH1 == 5 ~ 1.5 + 0.5*LBXD05LA) +
           rnorm(nrow(df), 0, 7))
}

# interxn in smaller group, larger effect size
em2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1.5, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5 + LBXD05LA, # double in group 4
                     RIDRETH1 == 5 ~ 1.5 + LBXD05LA) +
           rnorm(nrow(df), 0, 7))
}

# interxn in larger group, smaller effect size
ep1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1.5, 
                     RIDRETH1 == 3 ~ 1 + 0.5*LBXD05LA, # 1.5x in group 3
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1.5) +
           rnorm(nrow(df), 0, 7))

}

# interxn in larger group, larger effect size
ep2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2)*LBXPCBLA +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1.5,
                     RIDRETH1 == 3 ~ 1 + LBXD05LA, # double in group 3
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1.5) +
           rnorm(nrow(df), 0, 7))
}
