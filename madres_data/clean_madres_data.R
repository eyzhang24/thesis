library(tidyverse)
library(copula)
library(copulaedas)
library(GGally)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

target <- read_csv("madres_data/1945_TARGETED_DATA.csv")
epi <- read_csv("madres_data/1945_EPI_DATA.csv")

########
#clean target data
########

target_small <- target |> 
  #if below LOD, use LOD / sqrt(2)
  mutate(conc_mod = ifelse(Comment_code == 37, 
                           LOD / sqrt(2), 
                           Concentration)) |> 
  #adjust for urine specific gravity: Ac = A × [(SGmean –1)/(SG–1)]
  mutate(conc_mod = conc_mod * ((mean(target$SG)-1)/(SG-1))) |> 
  select(Project_ID, SID, PID, child_PID, Analyte_Code, conc_mod) |> 
  group_by(SID) |> 
  mutate(Project_ID = min(Project_ID)) |> 
  ungroup() |> 
  pivot_wider(names_from = Analyte_Code, values_from = conc_mod) |> 
  #howe kept As, Cd, Co, Hg, Ni, Tl, and Pb in main, Mo, Sb, and Sn in supp
  #don't have modified version of As used in their paper
  select(Project_ID, SID, PID, child_PID, As, Cd, Co, Hg, Ni, Tl, Pb, Mo, Sb, Sn)

#save
write_csv(target_small, "madres_data/target_small.csv")

#only keep data from first trimester
target_first <- target_small |> 
  group_by(child_PID) |> 
  filter(Project_ID == min(Project_ID)) |> 
  ungroup()

#correlation structure
cor(target_first[, 5:14])
# ggpairs(target_small, columns = 5:11)
target_first |> 
  select(5:14) |> 
  gather() |> 
  ggplot(aes(x = value)) +
  geom_density() + 
  facet_wrap(~key, scales = "free")

#save
write_csv(target_first, "madres_data/target_first.csv")

########
#experiment with copula
########

#spearman rho
cor(target_first[, 5:14], method = "spearman")

#create pseudo observations
u <- pobs(target_first[, 5:14])

#fit copulas
cfit_gaus <- fitCopula(normalCopula(dim = 7, dispstr = "un"), u)
cfit_t <- fitCopula(tCopula(dim = 7, dispstr = "un", df.fixed = TRUE), u)
cfit_gum <- fitCopula(gumbelCopula(dim = 7), u)
cfit_frank <- fitCopula(frankCopula(dim = 7), u)
cfit_clay <- fitCopula(claytonCopula(dim = 7), u)
cfit_joe <- fitCopula(joeCopula(dim = 7), u)

#evaluate fit
aic_values <- sapply(list(cfit_gaus, cfit_t, cfit_gum, cfit_frank, 
                          cfit_clay, cfit_joe), AIC)
names(aic_values) <- c("cfit_gaus", "cfit_t", "cfit_gum", "cfit_frank", 
                       "cfit_clay", "cfit_joe")
sort(aic_values)

#t copula performs best, proceed with this

########
#fit t copula and simulate data
########
cfit_t <- fitCopula(tCopula(dim = 10, dispstr = "un"), u)
#should df.fixed = TRUE (current default), or specified beforehand? 

#get rho and degrees of freedom
rho <- coef(cfit_t)[1:45]
df <- coef(cfit_t)[26]

##TEST FOR ONE ITERATION##

#simulate pseudo observations (quantile values)
set.seed(0)
samp <- rCopula(nrow(target_first),tCopula(rho,dim=10,df=df,dispstr = "un"))

#transform back to marginal distributions
sampt <- apply(samp, 2, 
               function(col) quantile(target_first[,5:14][[1]], probs = col))
sampt <- as.data.frame(sampt)
names(sampt) <- names(target_first[,5:14])

#check correlation structure
cor(samp, method = "spearman")
cor(sampt, method = "spearman")
cor(target_first[,5:14], method = "spearman")
#they all appear similar!

##MOVE ON TO MULTIPLE ITERATIONS##

#create function for simulation
simulate_data <- function(data, rho, df = 1) {
  #'data = original observed data
  #'rho = rho values from t-copula
  #'df = degrees of freedom from t-copula
  
  #simulate pseudo-observations from copula
  #currently using default df from original tcopula fit
  samp <- rCopula(nrow(data), 
                  tCopula(rho, dim = ncol(data), dispstr = "un", df = df))
  #, df.fixed = TRUE
  #transform pseudo-observations to observed marginal distributions
  sampt <- apply(samp, 2, 
                 function(col) quantile(data[[1]], probs = col)) |> 
    as.data.frame()
  names(sampt) <- names(data)
  return(sampt)
}

#create multiple simulated datasets
set.seed(0)
simulated <- 1:10 |> 
  purrr::map(\(x) simulate_data(target_first[,5:14], rho = rho, df = df))
# s1 <- simulated[[2]]

########
#clean epi data
########

#'Final models were adjusted for 
#'recruitment site, 
#'maternal age, 
#'pre- pregnancy BMI, 
#'race by ethnicity and birth place, 
#'any maternal anemia in pregnancy, 
#'any smoke exposure — maternal or other in pregnancy, 
#'and urinary AsB — an objective biomarker of fish and seafood consumption

#select relevant variables
epi_small <- epi |> 
  select(child_pid, mom_site, 
         t1_mat_age, #age, trimester 1 is most complete #t3_mat_age, birth_mat_age
         t1_pre_BMI, #bmi
         # baby_hispanic, baby_race, #baby r/e (less complete)
         t1_demo_race, t1_demo_hispanic, t1_demo_usa, #maternal m/e
         #ever-exposure to smoke
         t1_smoke, t1_smoke_preg, t2_smoke, t2_smoke_preg, t3_smoke, t3_smoke_preg, 
         #seafood intake
         t2_seafood_can, t2_seafood_fresh, t2_seafood_fried, t2_seafood_fsticks, 
         t2_seafood_oily, t2_seafood_shell, 
         t3_seafood_can, t3_seafood_fresh, t3_seafood_fried, t3_seafood_fsticks, 
         t3_seafood_oily, t3_seafood_shell, 
         birthweight #birthweight
         #can't find anemia measure or AsB
         ) |> 
  #replace -99 with NA
  mutate(across(where(is.numeric), ~ifelse(.==-99, NA, .))) 

#count number of NA values
apply(is.na(epi_small), 2, sum)

#handle NA values
epi_imp <- epi_small |> 
  #na's for smoke during preg, set to 0
  mutate(across(c(t1_smoke_preg, t2_smoke_preg, t3_smoke_preg), 
                ~ifelse(is.na(.), 0, .))) |> 
  #otherwise, impute median
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(.,na.rm = TRUE), .)))

########
#combine
########
comb <- epi |> 
  left_join(target_first, by = c("child_pid" = "child_PID"))
  
