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

# correlation structure
cor(target_first[, 5:14], method = "spearman")
# univariate density plots
target_first |> 
  select(5:14) |> 
  gather() |> 
  ggplot(aes(x = value)) +
  geom_density() + 
  facet_wrap(~key, scales = "free")

# log transformed
target_first |> 
  select(5:14) |> 
  # remove outliers
  filter(Mo >=1, Sb <= 1.4) |> 
  gather() |> 
  ggplot(aes(x = log(value))) +
  geom_density() + 
  facet_wrap(~key, scales = "free")

cor_mat <- cor(target_first[, 5:14], method = "spearman")
cor_mat[lower.tri(cor_mat)] <- NA
melt_cor <- reshape2::melt(cor_mat) |> 
  mutate(label = ifelse(value == 1, NA, round(value, 2)))
melt_cor |> 
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    limit = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = "Spearman's rho") +
  theme(
    panel.grid.major.x = element_line(color = "grey85",
                                      linewidth = 0.25,
                                      linetype = 2), 
    panel.border = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

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
           t1_demo_hispanic == 0 & t1_demo_race == 2 ~ 1, #non-hisp white
           t1_demo_hispanic == 0 & t1_demo_race == 4 ~ 2, #non-hisp black
           t1_demo_hispanic == 0 ~ 3, #other, non-hispanic
           t1_demo_hispanic == 1 & t1_demo_usa == 1 ~ 4, #hispanic in US
           t1_demo_hispanic == 1 & t1_demo_usa == 0 ~ 5, #hispanic NOT in US
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
  select(-c(gender, birthweight, GA, mom_site)) |> 
  #na's for smoke during preg, set to 0
  mutate(smoke = as.factor(ifelse(is.na(smoke), 0, smoke))) |> 
  #otherwise, impute mean
  mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(.,na.rm = TRUE), .))) 

#race
table(epi_imp$race, useNA = 'ifany')

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

#look at correlation between race and chemicals
comb_small |> 
  mutate(across(10:19, log)) |> 
  select(c(6, 8:19)) |> 
  pivot_longer(cols = 2:12, names_to = "key", values_to = "value") |> 
  mutate(race = as.factor(race)) |> 
  ggplot(aes(x = race, y = value, color = race)) +
  geom_boxplot() +
  facet_wrap(~key, scales = "free_y")

comb_small |> 
  mutate(across(10:19, log)) |> 
  select(c(7, 8:19)) |> 
  pivot_longer(cols = 2:12, names_to = "key", values_to = "value") |> 
  mutate(smoke = as.factor(smoke)) |> 
  ggplot(aes(x = smoke, y = value, color = smoke)) +
  geom_boxplot() +
  facet_wrap(~key, scales = "free_y")

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
cor(comb_log[, 7:19], method = "spearman")

#create pseudo observations
u <- pobs(comb_log[, 7:19])
#jitter smoke uniformly 
prop_smoke0 <- 1 - mean(comb_log$smoke)
set.seed(0)
u_smoke <- comb_log$smoke |> 
  map_dbl(\(x) {
    ifelse(x == 0, runif(1, 0, prop_smoke0), runif(1, prop_smoke0, 1))
  })
u[, 1] <- u_smoke
cor(u, method = "spearman")

#fit copulas
cfit_gaus <- fitCopula(normalCopula(dim = 13, dispstr = "un"), u)
cfit_t <- fitCopula(tCopula(dim = 13, dispstr = "un", df.fixed = FALSE), u)
# In var.mpl(copula, u) :
#   the covariance matrix of the parameter estimates is computed as if 
#   'df.fixed = TRUE' with df = 60.8950734746005
cfit_t2 <- fitCopula(tCopula(dim = 13, dispstr = "un", df = 4, df.fixed = TRUE), u)
cfit_t3 <- fitCopula(tCopula(dim = 13, dispstr = "un", df = 10, df.fixed = TRUE), u)

cfit_gum <- fitCopula(gumbelCopula(4, dim = 13), u)
cfit_frank <- fitCopula(frankCopula(4, dim = 13), u)
cfit_clay <- fitCopula(claytonCopula(4, dim = 13), u)
cfit_joe <- fitCopula(joeCopula(4, dim = 13), u)

cfit_gum2 <- fitCopula(gumbelCopula(2, dim = 13), u)
cfit_frank2 <- fitCopula(frankCopula(2, dim = 13), u)
cfit_clay2 <- fitCopula(claytonCopula(2, dim = 13), u)
cfit_joe2 <- fitCopula(joeCopula(2, dim = 13), u)

#evaluate fit
aic_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, cfit_clay, cfit_joe, 
                          cfit_gum2, cfit_frank2, cfit_clay2, cfit_joe2), AIC)
names(aic_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", "cfit_clay", "cfit_joe", 
                       "cfit_gum2", "cfit_frank2", "cfit_clay2", "cfit_joe2")
sort(aic_values)

lik_values <- sapply(list(cfit_gaus, cfit_t, cfit_t2, cfit_t3, #cfit_t4, cfit_t5,
                          cfit_gum, cfit_frank, cfit_clay, cfit_joe, 
                          cfit_gum2, cfit_frank2, cfit_clay2, cfit_joe2), logLik)
names(lik_values) <- c("cfit_gaus", "cfit_t", "cfit_t2", "cfit_t3", #"cfit_t4", "cfit_t5",
                       "cfit_gum", "cfit_frank", "cfit_clay", "cfit_joe", 
                       "cfit_gum2", "cfit_frank2", "cfit_clay2", "cfit_joe2")
sort(lik_values)

#guassian copula performs best, proceed with this
write_rds(cfit_gaus, "sim/gauscop.RDS")

########
#fit t copula and simulate data
########
cfit_gaus <- read_rds("sim/gauscop.RDS")
#should df.fixed = TRUE (current default), or specified beforehand? 

#get rho and degrees of freedom
rho <- coef(cfit_gaus)
# rho <- coef(cfit_gaus)

#create function for simulation
simulate_data <- function(data, n, rho, prop_smoke, prop_race) {
  #'data = original observed data
  #'n = sample size
  #'rho = rho values from normal copula
  #'prop_smoke = proportion smoke from observed dataset
  #'prop_race = table with race/eth values
  
  #simulate pseudo-observations from copula
  samp <- rCopula(n, 
                  normalCopula(rho, dim = ncol(data), dispstr = "un"))
  #transform pseudo-observations to observed marginal distributions
  sampt <- 1:ncol(data) |> 
    purrr::map_dfc(
      \(x) {
        if(names(data)[x] == "smoke") {
          df <- data.frame(ifelse(u_smoke < prop_smoke0, 0, 1), 
                           row.names = NULL)
        } else {
          df <- data.frame(quantile(data[[x]], probs = samp[,x]), 
                           row.names = NULL)
        }
        names(df) <- names(data)[x]
        return(df)
      }
    ) |> 
    mutate(race = sample(x = names(prop_race), prob = prop_race,
                         size = n, replace = T)) |> 
    relocate(race)
  return(sampt)
}

#create multiple simulated datasets
set.seed(0)
simulated <- 1:10 |> 
  purrr::map(\(x) {
    mutate(simulate_data(comb_log[,7:19], n = nrow(comb_log), rho = rho, 
                         prop_smoke = 1-mean(comb_log$smoke), 
                         prop_race = table(comb_log$race)), sim = x)
    })

comb_sim <- bind_rows(simulated)

#look at univariate distributions
comb_sim |> 
  mutate(sim = as.factor(sim)) |> 
  # mutate(across(age:Sn, scale)) |> 
  pivot_longer(cols = 1:13) |>
  ggplot(aes(x = value, group = sim)) +
  geom_density(color = "grey10", alpha = 0.0001) + 
  geom_density(
    data = comb_log |> select(7:19) |> pivot_longer(cols = 1:13),
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
set.seed(0)
df2 <- df |> 
  mutate(across(age:Sn, ~c(scale(.)))) |> 
  mutate(across(mom_site:smoke, ~as.factor(round(.)))) |> 
  mutate(
    r00 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r01 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      0.6*Hg*Ni + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r02 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      0.2*Hg*((Ni-1)^2) +
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r03 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      2*sin(0.25*pi*Hg*Ni) + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r11 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      Cd*As + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r12 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      0.2*Cd*((As-1)^2) +
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r13 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      2*sin(0.25*pi*Cd*As) + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r21 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      0.6*Co*Ni + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r22 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      0.2*Co*((Ni-1)^2) +
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
    r23 = 
      Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 
      2*sin(0.25*pi*Co*Ni) + 
      age + 0.5*bmi + 
      case_when(race == 1 ~ 2, 
                race == 2 ~ -0.5, 
                race == 3 ~ 1, 
                race == 4 ~ -0.25) +
      ifelse(smoke == 1, -1, 0.5) +
      rnorm(nrow(df), 0, 5), 
  ) |> 
  select(-sim)

# m1 <- lm(response ~ ., data = df2)
# summary(m1)  

#fit oracle models
m00 <- lm(r00 ~ Hg + Sb +
           I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
           age + bmi + race + smoke, data = df2)
summary(m00)$r.squared
m01 <- lm(r01 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            Hg*Ni + 
            age + bmi + race + smoke, data = df2)
summary(m01)$r.squared
summary(m01)$coefficients[,4] # p-values
m02 <- lm(r02 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(Hg*((Ni-1)^2)) + 
            age + bmi + race + smoke, data = df2)
summary(m02)$r.squared
m03 <- lm(r03 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(sin(0.25*pi*Hg*Ni)) + 
            age + bmi + race + smoke, data = df2)
summary(m03)$r.squared
m11 <- lm(r11 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            Hg*Ni + 
            age + bmi + race + smoke, data = df2)
summary(m11)$r.squared
m12 <- lm(r12 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(Hg*((Ni-1)^2)) + 
            age + bmi + race + smoke, data = df2)
summary(m12)$r.squared
m13 <- lm(r13 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(sin(0.25*pi*Hg*Ni)) + 
            age + bmi + race + smoke, data = df2)
summary(m13)$r.squared
m21 <- lm(r21 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            Hg*Ni + 
            age + bmi + race + smoke, data = df2)
summary(m21)$r.squared
m22 <- lm(r22 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(Hg*((Ni-1)^2)) + 
            age + bmi + race + smoke, data = df2)
summary(m22)$r.squared
m23 <- lm(r23 ~ Hg + Sb +
            I(1/(1+exp(-4*Ni))) + I(Sb^2) + I(1/(1+exp(-4*Sn))) +
            I(sin(0.25*pi*Hg*Ni)) + 
            age + bmi + race + smoke, data = df2)
summary(m23)$r.squared

#fit bkmr
Z <- df2 |> 
  select(As:Sn)
X <- bind_cols(
    data.frame(model.matrix(~ mom_site-1, data = df2)), 
    data.frame(model.matrix(~ race-1, data = df2)), 
    df2
  ) |> 
  select(mom_site2:mom_site5, race2:race4, smoke:bmi) |> 
  mutate(across(everything(), as.numeric))
bklist <- vector(mode='list', length=10)
set.seed(0)
for(i in 1:10) {
  y <- df2[,i+15]
  bk <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE)
  bklist[[i]] <- bk
}
write_rds(bklist, "madres_data/bk_test2.RDS")

pips <- bklist |> 
  purrr::map_df(\(x) {
    ExtractPIPs(x) |> 
      pivot_wider(names_from = variable, values_from = PIP) 
  }) |> 
  mutate(mod = names(df2)[16:25])

pips |> 
  pivot_longer(cols = As:Sn, names_to = "variable", values_to = "PIP") |> 
  ggplot(aes(x = variable, y = PIP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~mod, scales = "free_x")
ggsave("madres_data/bk_test2.png", width = 7, height = 5)

#model with no interactions
bk1 <- bklist[[1]]
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

# read back in first round of fitted datasets, simulate outcome, fit regression models
