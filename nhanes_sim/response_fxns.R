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
           rnorm(nrow(df), 0, 1))
}

am1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.35*LBXD05LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

am2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.7*LBXD05LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

ap1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.125*LBXD05LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

ap2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXD05LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

# univariately insignificant
bm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA +
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.35*LBXF08LA*LBXF03LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           
           rnorm(nrow(df), 0, 1))
}

bm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.7*LBXF08LA*LBXF03LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

bp1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.125*LBXF08LA*((LBXF03LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

bp2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.25*LBXF08LA*((LBXF03LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

cm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.35*LBX074LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

cm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.7*LBX074LA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

cp1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.125*LBX074LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

cp2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.125*LBX074LA*((LBX194LA-1)^2) +
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

dm1 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.35*LBXD05LA*LBXPCBLA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
}

dm2 <- function(df) {
  mutate(df, y = 
           LBXD05LA + LBX074LA + 
           3/(1+exp(-4*LBX194LA)) + 
           1.5/(1+exp(-4*LBXPCBLA)) - 
           0.75*(LBXF04LA^2) + 0.5*LBXF04LA + 
           0.7*LBXD05LA*LBXPCBLA*LBX194LA + 
           RIDAGEYR + 0.5*LBXLYPCT + 0.5*LBXMOPCT + 
           case_when(RIDRETH1 == 1 ~ 1.5, 
                     RIDRETH1 == 2 ~ 1, 
                     RIDRETH1 == 3 ~ 1, 
                     RIDRETH1 == 4 ~ 1.5, 
                     RIDRETH1 == 5 ~ 1) +
           rnorm(nrow(df), 0, 1))
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
           rnorm(nrow(df), 0, 1))
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
           rnorm(nrow(df), 0, 1))
}
