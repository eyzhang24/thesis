library(tidyverse)
library(latex2exp)

# set up ------------------------------------------------------------------

# set theme for plots
theme_set(theme_light())
theme_update(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

# create labeller
equations1 <-  c(TeX("No inter", output = "character"), 
                 TeX("0.35Hg$*$Ni"), TeX("0.13Hg$*($Ni$-1)^2$"), 
                 TeX("0.35Cd$*$As"), TeX("0.125Cd$*($As$-1)^2$"), 
                 TeX("0.3Ni$*$Co"), TeX("0.1Ni$*($Co$-1)^2$"), 
                 TeX("0.3Hg$*$Ni$*$Tl"), TeX("0.09Hg$*($Ni$-1)^2*$Tl"), 
                 TeX("0.7Hg$*$Ni"), TeX("0.26Hg$*($Ni$-1)^2$"), 
                 TeX("0.7Cd$*$As"), TeX("0.25Cd$*($As$-1)^2$"), 
                 TeX("0.6Ni$*$Co"), TeX("0.2Ni$*($Co$-1)^2$"), 
                 TeX("0.6Hg$*$Ni$*$Tl"), TeX("0.18Hg$*($Ni$-1)^2*$Tl"))
names1 <- c("_base",
            "am1", "ap1", "bm1", "bp1", "cm1", "cp1", "dm1", "dp1",
            "am2", "ap2", "bm2", "bp2", "cm2", "cp2", "dm2", "dp2")
appender1 <- function(string) {
  return(equations1[match(string, names1)])
}
namesa <- c(1, 2, 4, 6, 8, 10, 12, 14, 16, 3, 5, 7, 9, 11, 13, 15, 17)
appendera <- function(string) {
  return(equations1[match(string, namesa)])
}

# add coloring based on true significance
get_sign <- function(case, chem) {
  case_when(
    chem %in% c("Hg", "Ni", "Sb", "Sn") ~ TRUE, 
    case %in% 6:9 & chem %in% c("Cd", "As") ~ TRUE, 
    case %in% 10:13 & chem == "Co" ~ TRUE, 
    case %in% 14:17 & chem == "Tl" ~ TRUE, 
    .default = FALSE
  )
}

# add coloring based on true significance BSR
get_sign_bsr <- function(case, chem) {
  case_when(
    case %in% 2:5 & chem %in% c(4, 5) ~ TRUE, 
    case %in% 6:9 & chem %in% c(1, 2) ~ TRUE, 
    case %in% 10:13 & chem %in% c(3, 5) ~ TRUE, 
    case %in% 14:17 & chem %in% c(4, 5, 6) ~ TRUE, 
    .default = FALSE
  )
}

# naive mlr ---------------------------------------------------------------

# p-values small
nsm_pval <- read_csv("sim/_mlr/pvalsm.csv")
nsm_pvalc <- nsm_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
nsm_pvalc_sens <- nsm_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
nsm_pvalc |> 
  ggplot(aes(x = var, y = p)) +
  geom_bar(data = nsm_pvalc_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_nsm_univ_pval.png", width = 7.5, height = 5)

# p-values large
nlg_pval <- read_csv("sim/_mlr/pvallg.csv")
nlg_pvalc <- nlg_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
nlg_pvalc_sens <- nlg_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
nlg_pvalc |> 
  ggplot(aes(x = var, y = p)) +
  geom_bar(data = nlg_pvalc_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_nlg_univ_pval.png", width = 7.5, height = 5)

# oracle mlr --------------------------------------------------------------

# p-values small
osm_pval <- read_csv("sim/_oracle/pvalsm.csv")
osm_pvalc <- osm_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sb", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
osm_pvalc_sens <- osm_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
osm_pvalc |> 
  ggplot(aes(x = var, y = p)) +
  geom_bar(data = osm_pvalc_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_osm_univ_pval.png", width = 7.5, height = 5)

# interactions

keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')

osm_pvali <- osm_pval |> 
  filter(!(var %in% keepnames)) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
osm_pvali_sens <- osm_pvali |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
osm_pvali |> 
  ggplot(aes(x = "", y = p)) +
  geom_bar(data = osm_pvali_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = NULL) +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_osm_biv_pval.png", width = 7.5, height = 5)

# p-values large
olg_pval <- read_csv("sim/_oracle/pvallg.csv")
olg_pvalc <- olg_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sb", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
olg_pvalc_sens <- olg_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
olg_pvalc |> 
  ggplot(aes(x = var, y = p)) +
  geom_bar(data = olg_pvalc_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_olg_univ_pval.png", width = 7.5, height = 5)

# interactions

keepnames <- c('(Intercept)', 'Hg', 'Sb', 'I(1/(1 + exp(-4 * Ni)))', 'Ni', 
               'I(Sb^2)', 'I(1/(1 + exp(-4 * Sn)))', 'age', 'bmi', 
               'race2', 'race3', 'race4', 'race5', 'smoke1', 'Cd', 'As', 'Co')

olg_pvali <- olg_pval |> 
  filter(!(var %in% keepnames)) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
olg_pvali_sens <- olg_pvali |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
olg_pvali |> 
  ggplot(aes(x = "", y = p)) +
  geom_bar(data = olg_pvali_sens, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = NULL) +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_olg_biv_pval.png", width = 7.5, height = 5)

# smaller size bkmr -------------------------------------------------------

# plot pips bkmr small
ksm_pips <- read_csv("sim/bkmr_sm/pips.csv")

# # boxplot of pip values
# ksm_pips |> 
#   ggplot(aes(x = variable, y = PIP)) +
#   geom_boxplot() + 
#   facet_wrap(~case)
# 
# # proportion > 0.5
# ksm_pips |> 
#   mutate(imp = PIP >= 0.5) |> 
#   group_by(case, variable) |> 
#   summarize(sensitivity = sum(imp)) |> 
#   ggplot(aes(x = variable, y = sensitivity)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~case)


ksm_pip_sig <- ksm_pips |> 
  mutate(sign = get_sign(case, variable))
ksm_pip_sen <- ksm_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
ksm_pip_sig |> 
  ggplot(aes(x = variable)) +
  geom_bar(data = ksm_pip_sen, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly significant", 
       x = "Chemical")
ggsave("index/figures/ch4_ksm_univ_pips.png", width = 7.5, height = 5)

# plot univariate relationships bkmr small


# plot bivariate relationships bkmr small
ksm_biv <- read_csv("sim/bkmr_sm/biv_expresp.csv")

# plot Hg and Ni
ksm_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(variable1 == "Hg", "Hg by Ni", "Ni by Hg"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ksm_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
ksm_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(variable1 == "Cd", "Cd by As", "As by Cd"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-4, 4) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ksm_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
ksm_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(variable1 == "Ni", "Ni by Co", "Co by Ni"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ksm_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bkmr small
ksm_triv <- read_csv("sim/bkmr_sm/triv_expresp.csv")
ksm_triv <- ksm_triv |> 
  mutate(variable1 = case_when(
    z1_name == "Hg" ~ "Hg by Ni + Tl", 
    z1_name == "Ni" ~ "Ni by Hg + Tl", 
    z1_name == "Tl" ~ "Tl by Hg + Ni"), 
  quantile = as.factor(z23_q)) 

ksm_triv |> 
  ggplot(aes(z1_val, est, color = quantile)) +
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) + 
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2+3\nquantile")
ggsave("index/figures/ch4_ksm_triv_expresp.png", width = 6, height = 4)

# int vs. rest significant visualization
ksm_ints <- read_csv("sim/bkmr_sm/ints.csv")

ksm_ints <- ksm_ints |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

summary(ksm_ints$signif)

ksm_ints |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = variable, y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3) +
  ylim(0, 0.15) +
  facet_wrap(~case, scales = "free_x", 
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = "Chemical, with all others fixed")
ggsave("index/figures/ch4_ksm_int_rest.png", width = 7.5, height = 5)

# one vs. other significant visualization
ksm_intb <- read_csv("sim/bkmr_sm/int_bivar.csv")

ksm_intb <- ksm_intb |> 
  mutate(cond = paste0(z1, "+", z2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

table(ksm_intb$signif)

ksm_intb |> 
  group_by(case) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = "", y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3) +
  ylim(0, 0.15) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_ksm_int_biv.png", width = 6, height = 4)

# one vs. two others significant visualization
ksm_intt <- read_csv("sim/bkmr_sm/int_trivar.csv")
# glimpse(ksm_intt)

ksm_intt <- ksm_intt |> 
  mutate(cond = paste0(variable, " by ", fixedat1, "+", fixedat2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

ksm_intt |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |>
  ggplot(aes(x = cond, y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3) +
  ylim(0, 0.025) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_ksm_int_triv.png", width = 6, height = 4)

# bivar and trivar together
int_comb <- bind_rows(
  select(ksm_intb, cond, trial, case, signif), 
  select(ksm_intt, cond, trial, case, signif)
)

int_comb |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(cond, sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3) +
  ylim(0, 0.15) +
  facet_wrap(~case, scales = "free_x",
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_ksm_int_bitri.png", width = 7.5, height = 5)


# larger size bkmr --------------------------------------------------------


# plot pips bkmr large
klg_pips <- read_csv("sim/bkmr_lg/pips.csv")


klg_pip_sig <- klg_pips |> 
  mutate(sign = get_sign(case, variable))
klg_pip_sen <- klg_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
klg_pip_sig |> 
  ggplot(aes(x = variable)) +
  geom_bar(data = klg_pip_sen, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly significant", 
       x = "Chemical")
ggsave("index/figures/ch4_klg_univ_pips.png", width = 7.5, height = 5)

# plot univariate relationships bkmr large


# plot bivariate relationships bkmr large
klg_biv <- read_csv("sim/bkmr_lg/biv_expresp.csv")

# plot Hg and Ni
klg_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(variable1 == "Hg", "Hg by Ni", "Ni by Hg"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_klg_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
klg_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(variable1 == "Cd", "Cd by As", "As by Cd"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-4, 4) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_klg_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
klg_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(variable1 == "Ni", "Ni by Co", "Co by Ni"), 
         quantile = as.factor(quantile)) |> 
  ggplot(aes(z1, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_klg_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bkmr large
klg_triv <- read_csv("sim/bkmr_lg/triv_expresp.csv")
klg_triv <- klg_triv |> 
  mutate(variable1 = case_when(
    z1_name == "Hg" ~ "Hg by Ni + Tl", 
    z1_name == "Ni" ~ "Ni by Hg + Tl", 
    z1_name == "Tl" ~ "Tl by Hg + Ni"), 
    quantile = as.factor(z23_q)) 

klg_triv |> 
  ggplot(aes(z1_val, est, color = quantile)) +
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) + 
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2+3\nquantile")
ggsave("index/figures/ch4_klg_triv_expresp.png", width = 6, height = 4)

# int vs. rest significant visualization
klg_ints <- read_csv("sim/bkmr_lg/ints.csv")

klg_ints <- klg_ints |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

summary(klg_ints$signif)

klg_ints |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = variable, y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3, position = "identity") +
  ylim(0, 0.9) +
  facet_wrap(~case, scales = "free_x", 
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = "Chemical, with all others fixed")
ggsave("index/figures/ch4_klg_int_rest.png", width = 7.5, height = 5)

# one vs. other significant visualization
klg_intb <- read_csv("sim/bkmr_lg/int_bivar.csv")

klg_intb <- klg_intb |> 
  mutate(cond = paste0(z1, "+", z2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

table(klg_intb$signif)

klg_intb |> 
  group_by(case) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(x = "", y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity*100), vjust = -0.5, size = 3) +
  ylim(0, 1) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_klg_int_biv.png", width = 6, height = 4)

# one vs. two others significant visualization
klg_intt <- read_csv("sim/bkmr_lg/int_trivar.csv")
# glimpse(klg_intt)

klg_intt <- klg_intt |> 
  mutate(cond = paste0(variable, " by ", fixedat1, "+", fixedat2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

klg_intt |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |>
  ggplot(aes(x = cond, y = sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity), vjust = -0.5, size = 3) +
  ylim(0, 0.1) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_klg_int_triv.png", width = 6, height = 4)

# bivar and trivar together
int_comb <- bind_rows(
  select(klg_intb, cond, trial, case, signif),
  select(klg_intt, cond, trial, case, signif)
)

int_comb |> 
  group_by(case, cond) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  ggplot(aes(cond, sensitivity)) +
  geom_bar(stat = "identity", fill = "grey85") +
  geom_text(aes(label = sensitivity*100), vjust = -0.5, size = 3) +
  ylim(0, 1) +
  facet_wrap(~case, scales = "free_x",
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  labs(y = "Sensitivity", 
       x = NULL)
ggsave("index/figures/ch4_klg_int_bitri.png", width = 7.5, height = 5)


# smaller size bsr --------------------------------------------------------

# plot pip's bsr small
ssm_pips <- read_csv("sim/bsr_sm/pips.csv")

ssm_pip_sig <- ssm_pips |> 
  mutate(sign = get_sign(case, variable))
ssm_pip_sen <- ssm_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
ssm_pip_sig |> 
  ggplot(aes(x = variable)) +
  geom_bar(data = ssm_pip_sen, 
           aes(y = sensitivity), 
           stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly significant", 
       x = "Chemical")
ggsave("index/figures/ch4_ssm_univ_pips.png", width = 7.5, height = 5)

# plot bivariate pip's small
ssm_pipb <- read_csv("sim/bsr_sm/pip_biv.csv")

ssm_pipb |> 
  # filter(case == 2) |> 
  mutate(sign = get_sign_bsr(case, Var1) & get_sign_bsr(case, Var2), 
         inter = paste0(Var1, "_", Var2)) |> 
  ggplot(aes(x = inter)) + 
  geom_pointrange(aes(y = PIP, color = sign),
                  stat = "summary",
                  fun.min = function(z) {quantile(z, 0.25)},
                  fun.max = function(z) {quantile(z, 0.75)}, 
                  fun = median, 
                  size  = 0.05) +
  scale_color_manual(values = c("#92D0E4", "darkorange2")) +
  facet_wrap(~case, 
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = c(0.9, 0.1)) +
  labs(y = "PIP value distribution", 
       color = "Truly significant", 
       x = NULL)
ggsave("index/figures/ch4_ssm_biv_pips.png", width = 7.5, height = 5)

# plot trivariate pip's small
ssm_pipt <- read_csv("sim/bsr_sm/pip_triv.csv")

ssm_pipt |> 
  ggplot(aes(x = "")) +
  geom_pointrange(aes(y = PIP),
                  stat = "summary",
                  fun.min = function(z) {quantile(z, 0.25)},
                  fun.max = function(z) {quantile(z, 0.75)}, 
                  fun = median, 
                  size  = 0.05) +
  facet_wrap(~case, 
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = c(0.9, 0.1)) +
  labs(y = "PIP value distribution", 
       color = "Truly significant", 
       x = NULL)

# plot bivariate relationships bkmr small
ssm_biv <- read_csv("sim/bsr_sm/biv_expresp.csv")

# plot Hg and Ni
ssm_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(j1 == "Hg", "Hg by Ni", "Ni by Hg"), 
         quantile = as.factor(j2quant)) |> 
  ggplot(aes(j1val, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-3, 9) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ssm_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
ssm_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(j1 == "Cd", "Cd by As", "As by Cd"), 
         quantile = as.factor(j2quant)) |> 
  ggplot(aes(j1val, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  ylim(-2, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ssm_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
ssm_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(j1 == "Ni", "Ni by Co", "Co by Ni"), 
         quantile = as.factor(j2quant)) |> 
  ggplot(aes(j1val, est, color = quantile)) + 
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) +
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  # ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")
ggsave("index/figures/ch4_ssm_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bsr small
ssm_triv <- read_csv("sim/bsr_sm/triv_expresp.csv")
ssm_triv <- ssm_triv |> 
  mutate(variable1 = case_when(
    j1 == "Hg" ~ "Hg by Ni + Tl", 
    j1 == "Ni" ~ "Ni by Hg + Tl", 
    j1 == "Tl" ~ "Tl by Hg + Ni"), 
    quantile = as.factor(j23quant)) 

## WHAT TO DO ABOUT THIS?? ##

ssm_triv |> 
  ggplot(aes(j1val, est, color = quantile)) +
  geom_line(aes(group = interaction(trial, quantile)), alpha = 0.2) + 
  ggh4x::facet_grid2(variable1~case, scales = "free", independent = "x", 
                     labeller = labeller(
                       case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "palegreen3", "darkorange"), 
                     guide = guide_legend(override.aes = list(alpha = 1), 
                                          reverse = TRUE)) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2+3\nquantile")
# ggsave("index/figures/ch4_ssm_triv_expresp.png", width = 6, height = 4)


