library(tidyverse)
library(copula)
library(GGally)
library(bkmr)

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
  #make new variables
  mutate(mom_site = as.factor(mom_site), 
         race = as.factor(case_when(
           t1_demo_hispanic == 0 & t1_demo_race == 2 ~ 1, #white
           t1_demo_hispanic == 0 ~ 2, #other, non-hispanic
           t1_demo_hispanic == 1 & t1_demo_usa == 1 ~ 3, #hispanic in US
           t1_demo_hispanic == 1 & t1_demo_usa == 0 ~ 4, #hispanic NOT in US
           .default = NA
         )), 
         smoke = as.factor(ifelse(
           t1_smoke_preg == 1 | t2_smoke_preg == 1 | t3_smoke_preg == 1 |
             t1_smoke == 1 | t2_smoke == 1 | t3_smoke == 1, 1, 0
         ))) |> 
  #replace -99 with NA
  mutate(across(where(is.numeric), ~ifelse(. == -99, NA, .)))  |> 
  dplyr::select(child_pid, mom_site, 
                age = t1_mat_age, #age, trimester 1
                bmi = t1_pre_BMI, #bmi
                race, #maternal r/e
                smoke, #ever-exposure to smoke
                gender, birthweight, GA #birthweight and gestational age
                #can't find anemia measure or AsB
  ) 

#count number of NA values
apply(is.na(epi_small), 2, sum)

#handle NA values
epi_imp <- epi_small |> 
  #don't include birthweight
  select(-c(gender, birthweight, GA)) |> 
  #na's for smoke during preg, set to 0
  mutate(smoke = as.factor(ifelse(is.na(smoke), 0, smoke))) |> 
  #otherwise, impute mean
  mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(.,na.rm = TRUE), .))) 

########
#combine
########
comb <- epi_imp |> 
  left_join(target_first, by = c("child_pid" = "child_PID")) |> 
  relocate(child_pid, Project_ID, SID, PID, mom_site, race, smoke) 

#look at outliers
mo_out <- comb |> 
  filter(Mo < 1)
sb_out <- comb |> 
  filter(Sb > 1.4)

#remove outliers
comb_small <- comb |> 
  filter(Mo >=1, Sb <= 1.4)

write_csv(comb_small, "madres_data/base_data.csv")

########
#experiment with copula
########

comb_small <- read_csv("madres_data/base_data.csv")

#log-transform target data
comb_log <- comb_small |> 
  mutate(across(10:19, log)) |> 
  #factors back to numeric
  mutate(across(where(is.factor), as.numeric))

#spearman rho
cor(comb_log[, 5:19], method = "spearman")

#create pseudo observations
u <- pobs(comb_log[, 5:19])
cor(u, method = "spearman")

#fit copulas
cfit_gaus <- fitCopula(normalCopula(dim = 15, dispstr = "un"), u)
cfit_t <- fitCopula(tCopula(dim = 15, dispstr = "un", df.fixed = FALSE), u)
cfit_t2 <- fitCopula(tCopula(dim = 15, dispstr = "un", df = 3, df.fixed = FALSE), u)
cfit_t3 <- fitCopula(tCopula(dim = 15, dispstr = "un", df.fixed = TRUE), u)

cfit_gum <- fitCopula(gumbelCopula(4, dim = 15), u)
cfit_frank <- fitCopula(frankCopula(4, dim = 15), u)
cfit_clay <- fitCopula(claytonCopula(4, dim = 15), u)
cfit_joe <- fitCopula(joeCopula(4, dim = 15), u)

#evaluate fit
aic_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, 
                          cfit_clay, cfit_joe), AIC)
names(aic_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", 
                       "cfit_clay", "cfit_joe")
sort(aic_values)

lik_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, 
                          cfit_clay, cfit_joe), logLik)
names(lik_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", 
                       "cfit_clay", "cfit_joe")
sort(lik_values)

#t copula performs best, proceed with this

########
#fit t copula and simulate data
########
# cfit_t <- fitCopula(tCopula(dim = 10, dispstr = "un"), u)
#should df.fixed = TRUE (current default), or specified beforehand? 

#get rho and degrees of freedom
rho <- coef(cfit_t)[1:105]
df <- coef(cfit_t)[106]
# rho <- coef(cfit_gaus)

#create function for simulation
simulate_data <- function(data, rho, df = 1) {
  #'data = original observed data
  #'rho = rho values from t-copula
  #'df = degrees of freedom from t-copula
  
  #simulate pseudo-observations from copula
  samp <- rCopula(nrow(data), 
                  # normalCopula(rho, dim = ncol(data), dispstr = "un"))
                  tCopula(rho, dim = ncol(data), dispstr = "un", df = df))
  #transform pseudo-observations to observed marginal distributions
  sampt <- 1:ncol(data) |> 
    purrr::map_dfc(
      \(x) {
        df <- data.frame(quantile(data[[x]], probs = samp[,x]), 
                         row.names = NULL)
        names(df) <- names(data)[x]
        return(df)
      }
    )
  return(sampt)
}

#create multiple simulated datasets
set.seed(0)
simulated <- 1:10 |> 
  purrr::map(\(x) {
    mutate(simulate_data(comb_log[,5:19], rho = rho, df = df), sim = x)
    })

comb_sim <- bind_rows(simulated)

#look at univariate distributions
comb_sim |> 
  mutate(sim = as.factor(sim)) |> 
  # mutate(across(age:Sn, scale)) |> 
  pivot_longer(cols = 1:15) |>
  ggplot(aes(x = value, group = sim)) +
  geom_density(color = "grey10", alpha = 0.0001) + 
  geom_density(
    data = comb_log |> select(5:19) |> pivot_longer(cols = 1:15),
    mapping = aes(x = value), 
    color = "darkorange", inherit.aes = FALSE
  ) +
  facet_wrap(~name, scales = "free")

#look at ggpairs plots
ggpairs(comb_sim, aes(group = sim), columns = 1:15)
ggpairs(comb_log, columns = 5:19, 
        diag = list(continuous = ggally_box_no_facet), 
        upper = list(continuous = wrap("cor", method = "spearman")))

########
#experiment with simulating response variable
########

#set covariate coefficients
#choose sigma --> R^2 around 20% (moderate signal)
#howe found Hg, Ni, Sb, Sn significant --> 1+
#use two effect sizes --> *2
  #two-way b/t Hg and Ni (marg. sig) (3 forms), --> +3
  #two-way b/t Cd and As (not marg. sig) (3 forms) --> +3
  #two-way b/t Ni and Co (highly correlated) (3 forms) --> +3
    #3 forms - mult, poly, sinusoidal 
  #three-way b/t Hg, Ni, and Tl (2 forms) --> +2
    #2 forms - mult, poly
  #two-way b/t race/eth + Ni (1 form) --> +1
    #effect modulated by level
##totals (3+3+3+2+1)*2 + 1 = 25 versions 

df <- simulated[[1]]
set.seed(1)
df2 <- df |> 
  mutate(across(age:Sn, scale)) |> 
  mutate(
    response = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      age + 0.5*bmi +
      rnorm(nrow(df), 0, 5) 
  ) |> 
  select(-sim)
# df$response <- df2$response

m1 <- lm(response ~ ., data = df2)
summary(m1)  
m2 <- lm(response ~ Hg + Sb +
           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
           age + bmi, data = df2)
summary(m2)

#fit bkmr
y <- df2$response
Z <- df2 |> 
  select(As:Sn)
X <- df2 |> 
  select(mom_site:bmi)
bk1 <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE)

#eval fit
TracePlot(fit = bk1, par = "beta")
TracePlot(fit = bk1, par = "sigsq.eps")
TracePlot(fit = bk1, par = "r", comp = 1)
#pips!
ExtractPIPs(bk1)
#plot univariate 
pred.resp.univar <- PredictorResponseUnivar(fit = bk1)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")
