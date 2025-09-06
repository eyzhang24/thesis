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
                 TeX("0.475$X_{D05}*X_{194}$"), TeX("0.24$X_{D05}*(X_{194}-1)^2$"), 
                 TeX("0.5$X_{F08}*X_{F03}$"), TeX("0.17$X_{F08}*(X_{F03}-1)^2$"), 
                 TeX("0.36$X_{074}*X_{194}$"), TeX("0.25$X_{074}*(X_{094}-1)^2$"), 
                 TeX("0.3$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.125$X_{D05}*(X_{PCB}-1)^2*X_{194}$"), 
                 TeX("0.9$X_{D05}*X_{194}$"), TeX("0.48$X_{D05}*(X_{194}-1)^2$"), 
                 TeX("$X_{F08}*X_{F03}$"), TeX("0.34$X_{F08}*(X_{F03}-1)^2$"), 
                 TeX("0.6$X_{074}*X_{194}$"), TeX("0.5$X_{074}*(X_{194}-1)^2$"), 
                 TeX("0.72$X_{D05}*X_{PCB}*X_{194}$"), TeX("0.25$X_{D05}*(X_{PCB}-1)^2*X_{194}$"))
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
    chem %in% c("LBXD05LA", "LBX074LA", "LBX194LA", "LBXPCBLA", "LBXF04LA") ~ TRUE, 
    case %in% 6:9 & chem %in% c("LBXF08LA", "LBXF03LA") ~ TRUE, 
    # case %in% 10:13 & chem == "Co" ~ TRUE, 
    # case %in% 14:17 & chem == "Tl" ~ TRUE, 
    .default = FALSE
  )
}

# add coloring based on true significance BSR
get_sign_bsr <- function(case, chem) {
  case_when(
    case %in% 2:5 & chem %in% c(13, 8) ~ TRUE, 
    case %in% 6:9 & chem %in% c(18, 15) ~ TRUE, 
    case %in% 10:13 & chem %in% c(1, 8) ~ TRUE, 
    case %in% 14:17 & chem %in% c(13, 9, 8) ~ TRUE, 
    .default = FALSE
  )
}

# 1   074 
# 2   099 
# 3   138 
# 4   153 
# 5   170 
# 6   180 
# 7   187 
# 8   194 
# 9   PCB 
# 10   HXC 
# 11   118 
# 12   D03 
# 13   D05 
# 14   D07 
# 15   F03 
# 16   F04 
# 17   F05 
# 18   F08 


# base case ---------------------------------------------------------------

# smaller size
nsm_pval <- read_csv("nhanes_sim/_mlr/pvalsm.csv")
osm_pval <- read_csv("nhanes_sim/_oracle/pvalsm.csv")
ksm_pips <- read_csv("nhanes_sim/bkmr_sm/pips.csv")
ssm_pips <- read_csv("nhanes_sim/bsr_sm/pips.csv")

# p-value visualization
osm_pvalc <- osm_pval |>
  mutate(var = case_when(
    grepl("LBX194LA", var) ~ "LBX194LA*", 
    grepl("LBXPCBLA", var) ~ "LBX194LA*", 
    grepl("I\\(LBXF04LA", var) ~ "LBXF04LA^2", 
    .default = var
  )) |> 
  mutate(var = gsub("LBX", "", var), var = gsub("LA", "", var)) |> 
  filter(var %in% c('D05', '074', '194*', 'F04^2', 'F04', 'F08', 'F03', 'PCB')) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = TRUE, 
         mod = "Oracle MLR")
nsm_pvalc <- nsm_pval |> 
  filter(var %in% c('LBX074LA', 'LBX099LA', 'LBX138LA', 'LBX153LA', 'LBX170LA', 
                    'LBX180LA', 'LBX187LA', 'LBX194LA', 'LBXPCBLA', 'LBXHXCLA', 
                    'LBX118LA', 'LBXD03LA', 'LBXD05LA', 'LBXD07LA', 'LBXF03LA', 
                    'LBXF04LA', 'LBXF05LA', 'LBXF08LA')) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR") |> 
  mutate(var = substr(var, 4, 6)) 

base_mlrs <- bind_rows(nsm_pvalc, osm_pvalc) |> 
  filter(case == 1)

p1 <- base_mlrs |> 
  ggplot(aes(var, p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  labs(y = "p-value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical") +
  facet_grid(.~mod, scales = "free_x", space = "free") 

# pip visualization
base_pips <- bind_rows(
  mutate(ksm_pips, mod = "BKMR"), 
  mutate(ssm_pips, mod = "BSR")
) |> 
  filter(case == 1) |> 
  mutate(sign = get_sign(case, variable))

p2 <- base_pips |> 
  ggplot(aes(variable, PIP)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  labs(y = "PIP distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~mod, scales = "free_x") 


# larger size
nlg_pval <- read_csv("nhanes_sim/_mlr/pvallg.csv")
olg_pval <- read_csv("nhanes_sim/_oracle/pvallg.csv")
klg_pips <- read_csv("nhanes_sim/bkmr_lg/pips.csv")
slg_pips <- read_csv("nhanes_sim/bsr_lg/pips.csv")

# p-val visualization
nlg_pvalc <- nlg_pval  |> 
  filter(var %in% c('LBX074LA', 'LBX099LA', 'LBX138LA', 'LBX153LA', 'LBX170LA', 
                    'LBX180LA', 'LBX187LA', 'LBX194LA', 'LBXPCBLA', 'LBXHXCLA', 
                    'LBX118LA', 'LBXD03LA', 'LBXD05LA', 'LBXD07LA', 'LBXF03LA', 
                    'LBXF04LA', 'LBXF05LA', 'LBXF08LA')) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR") |> 
  mutate(var = substr(var, 4, 6)) 
olg_pvalc <- olg_pval |>
  mutate(var = case_when(
    grepl("LBX194LA", var) ~ "LBX194LA*", 
    grepl("LBXPCBLA", var) ~ "LBX194LA*", 
    grepl("I\\(LBXF04LA", var) ~ "LBXF04LA^2", 
    .default = var
  )) |> 
  mutate(var = gsub("LBX", "", var), var = gsub("LA", "", var)) |> 
  filter(var %in% c('D05', '074', '194*', 'F04^2', 'F04', 'F08', 'F03', 'PCB')) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = TRUE, 
         mod = "Oracle MLR")

base_mlrl <- bind_rows(nlg_pvalc, olg_pvalc) |> 
  filter(case == 1)

p3 <- base_mlrl |> 
  ggplot(aes(var, p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  labs(y = "p-value distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_grid(.~mod, scales = "free_x", space = "free") 

# pip visualization
base_pipl <- bind_rows(
  mutate(klg_pips, mod = "BKMR"), 
  mutate(slg_pips, mod = "BSR")
) |> 
  filter(case == 1) |> 
  mutate(sign = get_sign(case, variable))

p4 <- base_pipl |> 
  ggplot(aes(variable, PIP)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  labs(y = "PIP distribution", 
       color = "Truly significant", 
       x = "Chemical") +
  facet_wrap(~mod, scales = "free_x") 

# plot grid
pone <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  p4 + theme(legend.position = "none"), 
  labels = c("a", "c", "b", "d"), nrow = 2, rel_widths = c(4, 5)
)
pone

legend <- cowplot::get_legend(p1)

# stitch together
cowplot::plot_grid(pone, legend, rel_widths = c(5, 0.5))
ggsave("index/figures/ch4_basecasesig.png", width = 9, height = 5)

# create table of sensitivities
base_mlr <- bind_rows(
  mutate(base_mlrs, size = "Small"), 
  mutate(base_mlrl, size = "Large")
)

sum_base_mlr <- base_mlr |> 
  group_by(mod, size, var, sign) |> 
  summarize(sensitivity = sum(p < 0.05)/n())

base_bay <- bind_rows(
  mutate(base_pips, size = "Small"), 
  mutate(base_pipl, size = "Large")
)

sum_base_bay <- base_bay |> 
  rename(var = variable) |> 
  group_by(mod, size, var, sign) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) 

sum_base <- bind_rows(sum_base_mlr, sum_base_bay)
write_csv(sum_base, "index/data/base_sens.csv")

# naive mlr ---------------------------------------------------------------

# p-values small
nsm_pval <- read_csv("sim/_mlr/pvalsm.csv")
# nsm_pvalc <- nsm_pval |> 
#   filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
#                     "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
#   mutate(trial = case, 
#          case = ceiling(case/100), 
#          sign = get_sign(case, var), 
#          mod = "Naive MLR")
nsm_pvalc_sens <- nsm_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
nsm_pvalc |> 
  filter(case != 1) |> 
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "p-value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("nhanes_figs/ch4_nsm_univ_pval.png", width = 8.5, height = 6)

# p-values large
# nlg_pval <- read_csv("sim/_mlr/pvallg.csv")
# nlg_pvalc <- nlg_pval |> 
#   filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
#                     "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
#   mutate(trial = case, 
#          case = ceiling(case/100), 
#          sign = get_sign(case, var), 
#          mod = "Naive MLR")
nlg_pvalc_sens <- nlg_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
nlg_pvalc |> 
  filter(case != 1) |> 
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "p-value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("nhanes_figs/ch4_nlg_univ_pval.png", width = 8.5, height = 6)

# save sensitivity as table
naive_all <- bind_rows(
  mutate(nsm_pvalc, size = "Small", mod = "Naive MLR"),
  mutate(nlg_pvalc, size = "Large", mod = "Naive MLR")
) |> 
  group_by(mod, case, size, var, sign) |> 
  summarize(sensitivity = sum(p < 0.05)/n())
write_csv(naive_all, "nhanes_figs/naive_sens.csv")

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
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = TRUE)
osm_pvalc_sens <- osm_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = sum(imp)/n())

# line plot of p-values
# osm_pvalc |> 
#   filter(case != 1) |> 
#   ggplot(aes(x = var, y = p)) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
#   geom_pointrange(aes(color = sign), 
#                   stat = "summary",
#                   fun.min = function(z) {quantile(z,0.25)},
#                   fun.max = function(z) {quantile(z,0.75)},
#                   fun = median, 
#                   size = 0.2) +
#   scale_color_manual(values = c("darkorange2")) +
#   # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   labs(y = "p-value distribution", 
#        color = "Truly\nsignificant", 
#        x = "Chemical") +
#   facet_wrap(~case, 
#              labeller = labeller(
#                case = as_labeller(appendera, default = label_parsed))) 
# ggsave("index/figures/ch4_osm_univ_pval.png", width = 7.5, height = 5)

# interactions

keepnames <- c(
  '(Intercept)', 'LBXD05LA', 'LBX074LA', 'I(1/(1 + exp(-4 * LBX194LA)))', 
  'I(1/(1 + exp(-4 * LBXPCBLA)))', 'I(LBXF04LA^2)', 'LBXF04LA', 'BMXBMI', 
  'LBXCOT', 'LBXWBCSI', 'LBXLYPCT', 'LBXMOPCT', 'LBXNEPCT', 'LBXEOPCT', 
  'LBXBAPCT', 'RIDAGEYR', 'RIAGENDR1', 'RIDRETH12', 'RIDRETH13', 'RIDRETH14', 
  'RIDRETH15', 'DMDEDUC22', 'DMDEDUC23', 'DMDEDUC24', 'DMDEDUC25', 'LBX194LA', 
  'LBXF08LA', 'LBXF03LA')

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
# osm_pvali |> 
#   ggplot(aes(x = "", y = p)) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
#   geom_pointrange(stat = "summary",
#                   fun.min = function(z) {quantile(z,0.25)},
#                   fun.max = function(z) {quantile(z,0.75)},
#                   fun = median, 
#                   size = 0.2, color = "darkorange") +
#   labs(y = "p-value distribution", 
#        color = "Truly significant", 
#        x = NULL) +
#   facet_wrap(~case, 
#              labeller = labeller(
#                case = as_labeller(appendera, default = label_parsed))) 
# ggsave("index/figures/ch4_osm_biv_pval.png", width = 7.5, height = 5)

# put interactions and univ together
osm_comb <- bind_rows(
  mutate(osm_pvali, variable = "Int"), 
  filter(mutate(osm_pvalc[osm_pvalc$case != 1, ], variable = var), !(var %in% c("F08", "F03")))
)

osm_comb |> 
  mutate(variable = factor(variable, levels = c(unique(osm_pvalc$var), "Int")), 
         interaction = (variable == "Int")) |> 
  ggplot(aes(x = variable, y = p, color = interaction)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2"))+
  labs(y = "p-value distribution", 
       color = "Interaction", 
       x = "Term") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("nhanes_figs/ch4_osm_pval.png", width = 8.5, height = 6)

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
# olg_pvalc |> 
#   ggplot(aes(x = var, y = p)) +
#   geom_bar(data = olg_pvalc_sens, 
#            aes(y = sensitivity), 
#            stat = "identity", fill = "grey85") +
#   geom_pointrange(stat = "summary",
#                   fun.min = function(z) {quantile(z,0.25)},
#                   fun.max = function(z) {quantile(z,0.75)},
#                   fun = median, 
#                   size = 0.2) +
#   scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#         legend.position = c(0.95, 0.05)) +
#   labs(y = "p-value distribution", 
#        color = "Truly significant", 
#        x = "Chemical") +
#   facet_wrap(~case, 
#              labeller = labeller(
#                case = as_labeller(appendera, default = label_parsed))) 
# ggsave("index/figures/ch4_olg_univ_pval.png", width = 7.5, height = 5)

# interactions

keepnames <- c(
  '(Intercept)', 'LBXD05LA', 'LBX074LA', 'I(1/(1 + exp(-4 * LBX194LA)))', 
  'I(1/(1 + exp(-4 * LBXPCBLA)))', 'I(LBXF04LA^2)', 'LBXF04LA', 'BMXBMI', 
  'LBXCOT', 'LBXWBCSI', 'LBXLYPCT', 'LBXMOPCT', 'LBXNEPCT', 'LBXEOPCT', 
  'LBXBAPCT', 'RIDAGEYR', 'RIAGENDR1', 'RIDRETH12', 'RIDRETH13', 'RIDRETH14', 
  'RIDRETH15', 'DMDEDUC22', 'DMDEDUC23', 'DMDEDUC24', 'DMDEDUC25', 'LBX194LA', 
  'LBXF08LA', 'LBXF03LA')

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
# olg_pvali |> 
#   ggplot(aes(x = "", y = p)) +
#   geom_bar(data = olg_pvali_sens, 
#            aes(y = sensitivity), 
#            stat = "identity", fill = "grey85") +
#   geom_pointrange(stat = "summary",
#                   fun.min = function(z) {quantile(z,0.25)},
#                   fun.max = function(z) {quantile(z,0.75)},
#                   fun = median, 
#                   size = 0.2) +
#   scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#         legend.position = c(0.95, 0.05)) +
#   labs(y = "p-value distribution", 
#        color = "Truly significant", 
#        x = NULL) +
#   facet_wrap(~case, 
#              labeller = labeller(
#                case = as_labeller(appendera, default = label_parsed))) 
# ggsave("index/figures/ch4_olg_biv_pval.png", width = 7.5, height = 5)

# put together
olg_comb <- bind_rows(
  mutate(olg_pvali, variable = "Int"), 
  filter(mutate(olg_pvalc[olg_pvalc$case != 1, ], variable = var), !(var %in% c("F08", "F03")))
)

olg_comb |> 
  mutate(variable = factor(variable, levels = c(unique(olg_pvalc$var), "Int")), 
         interaction = (variable == "Int")) |> 
  ggplot(aes(x = variable, y = p, color = interaction)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2"))+
  labs(y = "p-value distribution", 
       color = "Interaction", 
       x = "Term") +
  facet_wrap(~case, 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("nhanes_figs/ch4_olg_pval.png", width = 8.5, height = 6)

# save sensitivity as table
oracle_all <- bind_rows(
  mutate(osm_comb, size = "Small", mod = "Oracle MLR"),
  mutate(olg_comb, size = "Large", mod = "Oracle MLR")
) |> 
  group_by(mod, case, size, var, variable) |> 
  summarize(sensitivity = sum(p < 0.05)/n())
write_csv(oracle_all, "nhanes_figs/oracle_sens.csv")

# smaller size bkmr -------------------------------------------------------

# plot pips bkmr small
ksm_pips <- read_csv("nhanes_sim/bkmr_sm/pips.csv")

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
  filter(case != 1) |> 
  mutate(variable = substr(variable, 4, 6)) |> 
  ggplot(aes(x = variable)) +
  # geom_bar(data = ksm_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  #       legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")
ggsave("nhanes_figs/ch4_ksm_univ_pips.png", width = 8.5, height = 6)

# plot bivariate relationships bkmr small
ksm_biv <- read_csv("nhanes_sim/bkmr_sm/biv_expresp.csv")

# plot Hg and Ni
ksm_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(variable1 == "LBXD05LA", "D05 by 194", "194 by D05"), 
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
ggsave("nhanes_figs/ch4_ksm_biv_expresp_1.png", width = 7.5, height = 5)

# plot Cd and As
ksm_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(variable1 == "LBXF03LA", "F03 by F08", "F08 by F03"), 
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
ggsave("nhanes_figs/ch4_ksm_biv_expresp_2.png", width = 7.5, height = 5)

# plot Ni and Co
ksm_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(variable1 == "LBX074LA", "074 by 194", "194 by 074"), 
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
ggsave("nhanes_figs/ch4_ksm_biv_expresp_3.png", width = 7.5, height = 5)

# plot trivariate relationships bkmr small
ksm_triv <- read_csv("nhanes_sim/bkmr_sm/triv_expresp.csv")
ksm_triv <- ksm_triv |> 
  mutate(variable1 = case_when(
    z1_name == "LBXD05LA" ~ "D05 by PCB + 194", 
    z1_name == "LBXPCBLA" ~ "PCB by D05 + 194", 
    z1_name == "LBX194LA" ~ "194 by D05 + PCB"), 
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
ggsave("nhanes_figs/ch4_ksm_triv_expresp.png", width = 7.5, height = 5)

# int vs. rest significant visualization
ksm_ints <- read_csv("nhanes_sim/bkmr_sm/ints.csv")

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
ggsave("nhanes_figs/ch4_ksm_int_rest.png", width = 7.5, height = 5)

# one vs. other significant visualization
ksm_intb <- read_csv("nhanes_sim/bkmr_sm/int_bivar.csv")

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
ggsave("nhanes_figs/ch4_ksm_int_biv.png", width = 6, height = 4)

# one vs. two others significant visualization
ksm_intt <- read_csv("nhanes_sim/bkmr_sm/int_trivar.csv")
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
ggsave("nhanes_figs/ch4_ksm_int_triv.png", width = 6, height = 4)

# bivar and trivar together
int_combs <- bind_rows(
  select(ksm_intb, cond, trial, case, signif), 
  select(ksm_intt, cond, trial, case, signif)
)

int_combs |> 
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
ggsave("nhanes_figs/ch4_ksm_int_bitri.png", width = 7.5, height = 5)

# save sensitivity as table

# fdr's for bivariate and trivariate
ksm_bifdrs <- read_csv("nhanes_sim/bkmr_sm/int_bivar_full.csv")
ksm_trifdrs <- read_rds("nhanes_sim/bkmr_sm/tri_fdrs.RDS")


# larger size bkmr --------------------------------------------------------


# plot pips bkmr large
klg_pips <- read_csv("nhanes_sim/bkmr_lg/pips.csv")

klg_pip_sig <- klg_pips |> 
  mutate(sign = get_sign(case, variable))
klg_pip_sen <- klg_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
klg_pip_sig |> 
  filter(case != 1) |> 
  mutate(variable = substr(variable, 4, 6)) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = klg_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        # legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")
ggsave("nhanes_figs/ch4_klg_univ_pips.png", width = 8.5, height = 6)

# plot univariate relationships bkmr large


# plot bivariate relationships bkmr large
klg_biv <- read_csv("nhanes_sim/bkmr_lg/biv_expresp.csv")

# plot Hg and Ni
klg_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(variable1 == "LBXD05LA", "D05 by 194", "194 by D05"), 
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
ggsave("nhanes_figs/ch4_klg_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
klg_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(variable1 == "LBXF03LA", "F03 by F08", "F08 by F03"), 
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
ggsave("nhanes_figs/ch4_klg_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
klg_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(variable1 == "LBX074LA", "074 by 194", "194 by 074"), 
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
ggsave("nhanes_figs/ch4_klg_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bkmr large
klg_triv <- read_csv("nhanes_sim/bkmr_lg/triv_expresp.csv")
klg_triv <- klg_triv |> 
  mutate(variable1 = case_when(
    z1_name == "LBXD05LA" ~ "D05 by PCB + 194", 
    z1_name == "LBXPCBLA" ~ "PCB by D05 + 194", 
    z1_name == "LBX194LA" ~ "194 by D05 + PCB"), 
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
ggsave("nhanes_figs/ch4_klg_triv_expresp.png", width = 6, height = 4)

# int vs. rest significant visualization
klg_ints <- read_csv("nhanes_sim/bkmr_lg/ints.csv")

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
ggsave("nhanes_figs/ch4_klg_int_rest.png", width = 7.5, height = 5)

# one vs. other significant visualization
klg_intb <- read_csv("nhanes_sim/bkmr_lg/int_bivar.csv")

klg_intb <- klg_intb |> 
  mutate(across(c(z1, z2), ~substr(., 4, 6))) |> 
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
ggsave("nhanes_figs/ch4_klg_int_biv.png", width = 6, height = 4)

# one vs. two others significant visualization
klg_intt <- read_csv("nhanes_sim/bkmr_lg/int_trivar.csv")
# glimpse(klg_intt)

klg_intt <- klg_intt |> 
  mutate(across(c(variable, fixedat1, fixedat2), ~substr(., 4, 6))) |> 
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
ggsave("nhanes_figs/ch4_klg_int_triv.png", width = 6, height = 4)

# bivar and trivar together
int_combl <- bind_rows(
  select(klg_intb, cond, trial, case, signif),
  select(klg_intt, cond, trial, case, signif)
)

int_combl |> 
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
ggsave("nhanes_figs/ch4_klg_int_bitri.png", width = 7.5, height = 5)




# bkmr combine sensitivity ------------------------------------------------

# univariate pips
bkmr_pips <- bind_rows(
  mutate(ksm_pips, size = "Small"), 
  mutate(klg_pips, size = "Large")
) |> 
  group_by(size, case, variable) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) |> 
  mutate(sign = get_sign(case, variable))
write_csv(bkmr_pips, "nhanes_figs/bkmr_pip_sens.csv")

# interactions
bkmr_comb <- bind_rows(
  mutate(int_combs, size = "Small"), 
  mutate(int_combl, size = "Large")
)

bkmr_sens <- bkmr_comb |> 
  group_by(size, case, cond) |> 
  summarize(sensitivity = sum(signif)/n())
write_csv(bkmr_sens, "nhanes_figs/bkmr_int_sens.csv")

# smaller size bsr --------------------------------------------------------

# plot pip's bsr small
ssm_pips <- read_csv("nhanes_sim/bsr_sm/pips.csv")

ssm_pip_sig <- ssm_pips |> 
  mutate(sign = get_sign(case, variable))
ssm_pip_sen <- ssm_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
ssm_pip_sig |> 
  filter(case != 1) |> 
  mutate(variable = substr(variable, 4, 6)) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = ssm_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        # legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")
ggsave("nhanes_figs/ch4_ssm_univ_pips.png", width = 7.5, height = 5)

# plot bivariate pip's small
ssm_pipb <- read_csv("nhanes_sim/bsr_sm/pip_biv.csv")

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
ggsave("nhanes_figs/ch4_ssm_biv_pips.png", width = 7.5, height = 5)

# plot trivariate pip's small
ssm_pipt <- read_csv("nhanes_sim/bsr_sm/pip_triv.csv")

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
ssm_biv <- read_csv("nhanes_sim/bsr_sm/biv_expresp.csv")

# plot Hg and Ni
ssm_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(j1 == "LBX194LA", "194 by D05", "D05 by 194"), 
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
ggsave("nhanes_figs/ch4_ssm_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
ssm_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(j1 == "LBXF08LA", "F08 by F03", "F03 by F08"), 
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
ggsave("nhanes_figs/ch4_ssm_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
ssm_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(j1 == "LBX074LA", "074 by 194", "194 by 074"), 
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
ggsave("nhanes_figs/ch4_ssm_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bsr small
ssm_triv <- read_csv("nhanes_sim/bsr_sm/triv_expresp.csv")
ssm_triv <- ssm_triv |> 
  mutate(variable1 = case_when(
    j1 == "LBXD05LA" ~ "D05 by PCB + 194", 
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


# larger size bsr ---------------------------------------------------------

# plot pip's bsr large
slg_pips <- read_csv("nhanes_sim/bsr_lg/pips.csv")

slg_pip_sig <- slg_pips |> 
  mutate(sign = get_sign(case, variable))
slg_pip_sen <- slg_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = sum(imp)/n())

# point and lineplot
slg_pip_sig |> 
  filter(case != 1) |> 
  mutate(variable = substr(variable, 4, 6)) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = slg_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # legend.position = c(0.95, 0.05)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")
ggsave("nhanes_Figs/ch4_slg_univ_pips.png", width = 7.5, height = 5)

# plot bivariate pip's large
slg_pipb <- read_csv("nhanes_sim/bsr_lg/pip_biv.csv")

slg_pipb |> 
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
ggsave("nhanes_figs/ch4_slg_biv_pips.png", width = 7.5, height = 5)

# plot trivariate pip's large
slg_pipt <- read_csv("nhanes_sim/bsr_lg/pip_triv.csv")

slg_pipt |> 
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

# plot bivariate relationships bkmr large
slg_biv <- read_csv("nhanes_sim/bsr_lg/biv_expresp.csv")

# plot Hg and Ni
slg_biv |> 
  filter(case %in% 2:5) |> 
  mutate(variable1 = ifelse(j1 == "LBX194LA", "194 by D05", "D05 by 194"), 
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
ggsave("nhanes_Figs/ch4_slg_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
slg_biv |> 
  filter(case %in% 6:9) |> 
  mutate(variable1 = ifelse(j1 == "LBXF08LA", "F08 by F03", "F03 by F08"), 
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
ggsave("nhanes_figs/ch4_slg_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
slg_biv |> 
  filter(case %in% 10:13) |> 
  mutate(variable1 = ifelse(j1 == "LBX074LA", "074 by 194", "194 by 074"), 
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
ggsave("nhanes_figs/ch4_slg_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bsr large
slg_triv <- read_csv("nhanes_sim/bsr_lg/triv_expresp.csv")
slg_triv <- slg_triv |> 
  mutate(variable1 = case_when(
    j1 == "Hg" ~ "Hg by Ni + Tl", 
    j1 == "Ni" ~ "Ni by Hg + Tl", 
    j1 == "Tl" ~ "Tl by Hg + Ni"), 
    quantile = as.factor(j23quant)) 

## WHAT TO DO ABOUT THIS?? ##

slg_triv |> 
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
# ggsave("index/figures/ch4_slg_triv_expresp.png", width = 6, height = 4)


# bsr combine sens --------------------------------------------------------

# univariate pips
bsr_pips <- bind_rows(
  mutate(ssm_pips, size = "Small"), 
  mutate(slg_pips, size = "Large")
) |> 
  group_by(size, case, variable) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) |> 
  mutate(sign = get_sign(case, variable))
write_csv(bsr_pips, "nhanes_figs/bsr_pip_sens.csv")

# interactions
# cnames <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")
cnames <- c('LBX074LA', 'LBX099LA', 'LBX138LA', 'LBX153LA', 'LBX170LA', 
            'LBX180LA', 'LBX187LA', 'LBX194LA', 'LBXPCBLA', 'LBXHXCLA', 'LBX118LA', 
            'LBXD03LA', 'LBXD05LA', 'LBXD07LA', 'LBXF03LA', 'LBXF04LA', 'LBXF05LA', 
            'LBXF08LA')

bsr_comb <- bind_rows(
  mutate(ssm_pipb, size = "Small"), 
  mutate(slg_pipb, size = "Large")
) |> 
  mutate(sign = get_sign_bsr(case, Var1) & get_sign_bsr(case, Var2), 
         v1 = cnames[Var1], v2 = cnames[Var2], 
         v1 = substr(v1, 4, 6), v2 = substr(v2, 4, 6), 
         inter = paste0(v1, "-", v2)) |> 
  mutate(inter2 = ifelse(sign, inter, "none")) |> 
  group_by(size, case, inter2, sign) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n())
write_csv(bsr_comb, "nhanes_figs/bsr_int_sens.csv")
