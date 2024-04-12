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


# base case ---------------------------------------------------------------

# smaller size
nsm_pval <- read_csv("sim/_mlr/pvalsm.csv")
osm_pval <- read_csv("sim/_oracle/pvalsm.csv")
ksm_pips <- read_csv("sim/bkmr_sm/pips.csv")
ssm_pips <- read_csv("sim/bsr_sm/pips.csv")

# p-value visualization
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
         sign = TRUE, 
         mod = "Oracle MLR")
nsm_pvalc <- nsm_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR")

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
nlg_pval <- read_csv("sim/_mlr/pvallg.csv")
olg_pval <- read_csv("sim/_oracle/pvallg.csv")
klg_pips <- read_csv("sim/bkmr_lg/pips.csv")
slg_pips <- read_csv("sim/bsr_lg/pips.csv")

# p-val visualization
nlg_pvalc <- nlg_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR")
olg_pvalc <- olg_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sn*", "Sb^2")) |> 
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
nsm_pvalc <- nsm_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR")
nsm_pvalc_sens <- nsm_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(p, 0.75))

# line plot of p-values
nsm_pvalc |> 
  filter(case != 1) |> 
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_label(data = filter(nsm_pvalc_sens, case != 1), 
             mapping = aes(x = var, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
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
ggsave("index/figures/ch4_nsm_univ_pval.png", width = 7.5, height = 5)

# p-values large
nlg_pval <- read_csv("sim/_mlr/pvallg.csv")
nlg_pvalc <- nlg_pval |> 
  filter(var %in% c("As", "Cd", "Co", "Hg", "Ni", 
                    "Tl", "Pb", "Mo", "Sb", "Sn")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var), 
         mod = "Naive MLR")
nlg_pvalc_sens <- nlg_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(p, 0.75))

# line plot of p-values
nlg_pvalc |> 
  filter(case != 1) |> 
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_label(data = filter(nlg_pvalc_sens, case != 1), 
           mapping = aes(x = var, y = upper+0.1, label = sensitivity), 
           size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
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
ggsave("index/figures/ch4_nlg_univ_pval.png", width = 7.5, height = 5)

# save sensitivity as table
naive_all <- bind_rows(
  mutate(nsm_pvalc, size = "Small", mod = "Naive MLR"),
  mutate(nlg_pvalc, size = "Large", mod = "Naive MLR")
) |> 
  group_by(mod, case, size, var, sign) |> 
  summarize(sensitivity = sum(p < 0.05)/n())
write_csv(naive_all, "index/data/naive_sens.csv")

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
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(p, 0.75))

# line plot of p-values
osm_pvalc |>
  filter(case != 1) |>
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_label(data = filter(osm_pvalc_sens, case != 1), 
             mapping = aes(x = var, y = upper+0.12, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
  geom_pointrange(aes(color = sign),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median,
                  size = 0.2, color = "darkorange2") +
  # scale_color_manual(values = c("darkorange2")) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "p-value distribution",
       color = "Truly\nsignificant",
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

# density plot of p-values
osm_pvali |>
  ggplot(aes(x = p)) +
  geom_density(color = "darkorange2", fill = "darkorange", alpha = 0.1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey30") +
  labs(x = "p-value distribution") +
  facet_wrap(~case, scales = "free_y", 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed)))
ggsave("index/figures/ch4_osm_biv_pval.png", width = 7.5, height = 5)

# put interactions and univ together
osm_comb <- bind_rows(
  mutate(osm_pvali, variable = "Int"), 
  mutate(osm_pvalc[osm_pvalc$case != 1, ], variable = var)
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
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_osm_pval.png", width = 7.5, height = 5)

# p-values large
olg_pval <- read_csv("sim/_oracle/pvallg.csv")
olg_pvalc <- olg_pval |> 
  mutate(var = case_when(
    grepl("Ni", var) ~ "Ni*", 
    grepl("Sn", var) ~ "Sn*", 
    grepl("I\\(Sb", var) ~ "Sb^2", 
    .default = var
  )) |> 
  filter(var %in% c("Hg", "Ni*", "Tl", "Pb", "Mo", "Sn*", "Sb^2")) |> 
  mutate(trial = case, 
         case = ceiling(case/100), 
         sign = get_sign(case, var))
olg_pvalc_sens <- olg_pvalc |> 
  mutate(imp = p < 0.05) |> 
  group_by(case, var) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(p, 0.75))

# line plot of p-values
olg_pvalc |>
  filter(case != 1) |> 
  ggplot(aes(x = var, y = p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  # geom_bar(data = olg_pvalc_sens,
  #          aes(y = sensitivity),
  #          stat = "identity", fill = "grey85") +
  geom_label(data = filter(olg_pvalc_sens, case != 1), 
             mapping = aes(x = var, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
  geom_pointrange(stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median,
                  size = 0.2, color = "darkorange2") +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Sensitivity")) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  #       legend.position = c(0.95, 0.05)) +
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

# density plot of p-values
olg_pvali |>
  ggplot(aes(x = p)) +
  geom_density(color = "darkorange2", fill = "darkorange", alpha = 0.1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey30") +
  labs(x = "p-value distribution") +
  facet_wrap(~case, scales = "free_y", 
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed)))
ggsave("index/figures/ch4_olg_biv_pval.png", width = 7.5, height = 5)

# put together
olg_comb <- bind_rows(
  mutate(olg_pvali, variable = "Int"), 
  mutate(olg_pvalc[olg_pvalc$case != 1, ], variable = var)
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
               case = as_labeller(appendera, default = label_parsed))) 
ggsave("index/figures/ch4_olg_pval.png", width = 7.5, height = 5)

# save sensitivity as table
oracle_all <- bind_rows(
  mutate(osm_comb, size = "Small", mod = "Oracle MLR"),
  mutate(olg_comb, size = "Large", mod = "Oracle MLR")
) |> 
  group_by(mod, case, size, var, variable) |> 
  summarize(sensitivity = sum(p < 0.05)/n())
write_csv(oracle_all, "index/data/oracle_sens.csv")

# smaller size bkmr -------------------------------------------------------

# plot pips bkmr small
ksm_pips <- read_csv("sim/bkmr_sm/pips.csv")

ksm_pip_sig <- ksm_pips |> 
  mutate(sign = get_sign(case, variable))
ksm_pip_sen <- ksm_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(PIP, 0.75))

# point and lineplot
ksm_pip_sig |> 
  filter(case != 1) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_label(data = filter(ksm_pip_sen, case != 1), 
             mapping = aes(x = variable, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0, "lines")) +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_wrap(~case,
             labeller = as_labeller(appendera, 
                                    default = label_parsed)) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")
ggsave("index/figures/ch4_ksm_univ_pips.png", width = 7.5, height = 5)

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
       color = "Chem 2+3\nquantile") +
  theme(strip.text.x = element_text(size = 7))
ggsave("index/figures/ch4_ksm_triv_expresp.png", width = 6, height = 4)

# int vs. rest significant visualization
ksm_ints <- read_csv("sim/bkmr_sm/ints.csv")

ksm_ints <- ksm_ints |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

summary(ksm_ints$signif)

# one vs. other significant visualization
ksm_intb <- read_csv("sim/bkmr_sm/int_bivar.csv")

ksm_intb <- ksm_intb |> 
  mutate(cond = paste0(z1, "+", z2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

# one vs. two others significant visualization
ksm_intt <- read_csv("sim/bkmr_sm/int_trivar.csv")
# glimpse(ksm_intt)

ksm_intt <- ksm_intt |> 
  mutate(cond = paste0(variable, " by ", fixedat1, "+", fixedat2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))


# bivar and trivar together
int_combs <- bind_rows(
  select(ksm_intb, cond, trial, case, signif), 
  select(ksm_intt, cond, trial, case, signif)
)

# save sensitivity as table


# larger size bkmr --------------------------------------------------------


# plot pips bkmr large
klg_pips <- read_csv("sim/bkmr_lg/pips.csv")

klg_pip_sig <- klg_pips |> 
  mutate(sign = get_sign(case, variable))
klg_pip_sen <- klg_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2), 
            upper = quantile(PIP, 0.75))

# point and lineplot
klg_pip_sig |> 
  filter(case != 1) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = klg_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_label(data = filter(klg_pip_sen, case != 1), 
             mapping = aes(x = variable, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(.05, "lines")) +
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
       color = "Chem 2+3\nquantile") +
  theme(strip.text.x = element_text(size = 7))
ggsave("index/figures/ch4_klg_triv_expresp.png", width = 6, height = 4)

# int vs. rest significant visualization
klg_ints <- read_csv("sim/bkmr_lg/ints.csv")

klg_ints <- klg_ints |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

summary(klg_ints$signif)

# one vs. other significant visualization
klg_intb <- read_csv("sim/bkmr_lg/int_bivar.csv")

klg_intb <- klg_intb |> 
  mutate(cond = paste0(z1, "+", z2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

table(klg_intb$signif)


# one vs. two others significant visualization
klg_intt <- read_csv("sim/bkmr_lg/int_trivar.csv")
# glimpse(klg_intt)

klg_intt <- klg_intt |> 
  mutate(cond = paste0(variable, " by ", fixedat1, "+", fixedat2)) |> 
  rowwise() |> 
  mutate(signif = ifelse(
    between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
  ))

# bivar and trivar together
int_combl <- bind_rows(
  select(klg_intb, cond, trial, case, signif),
  select(klg_intt, cond, trial, case, signif)
)


# bkmr combine sensitivity ------------------------------------------------

# univariate pips
bkmr_pips <- bind_rows(
  mutate(ksm_pips, size = "Small"), 
  mutate(klg_pips, size = "Large")
) |> 
  group_by(size, case, variable) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) |> 
  mutate(sign = get_sign(case, variable))
write_csv(bkmr_pips, "index/data/bkmr_pip_sens.csv")

# interactions
bkmr_comb <- bind_rows(
  mutate(int_combs, size = "Small"), 
  mutate(int_combl, size = "Large")
)

bkmr_sens <- bkmr_comb |> 
  group_by(size, case, cond) |> 
  summarize(sensitivity = sum(signif)/n())
write_csv(bkmr_sens, "index/data/bkmr_int_sens.csv")

# smaller size bsr --------------------------------------------------------

# plot pip's bsr small
ssm_pips <- read_csv("sim/bsr_sm/pips.csv")

ssm_pip_sig <- ssm_pips |> 
  mutate(sign = get_sign(case, variable))
ssm_pip_sen <- ssm_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2, 
                                      max = 2), 
            upper = quantile(PIP, 0.75))


# point and lineplot
ssm_pip_sig |> 
  filter(case != 1) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = ssm_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_label(data = filter(ssm_pip_sen, case != 1), 
             mapping = aes(x = variable, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
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
       color = "Chem 2+3\nquantile") +
  theme(strip.text.x = element_text(size = 7))
ggsave("index/figures/ch4_ssm_triv_expresp.png", width = 6, height = 4)


# larger size bsr ---------------------------------------------------------

# plot pip's bsr large
slg_pips <- read_csv("sim/bsr_lg/pips.csv")

slg_pip_sig <- slg_pips |> 
  mutate(sign = get_sign(case, variable))
slg_pip_sen <- slg_pips |> 
  mutate(imp = PIP >= 0.5) |> 
  group_by(case, variable) |> 
  summarize(sensitivity = weights::rd(sum(imp)/n(), 2, max = 0), 
            upper = quantile(PIP, 0.75))

# point and lineplot
slg_pip_sig |> 
  filter(case != 1) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  # geom_bar(data = slg_pip_sen, 
  #          aes(y = sensitivity), 
  #          stat = "identity", fill = "grey85") +
  geom_label(data = filter(slg_pip_sen, case != 1), 
             mapping = aes(x = variable, y = upper+0.1, label = sensitivity), 
             size = 2, label.size = NA, label.padding = unit(0.05, "lines")) +
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
ggsave("index/figures/ch4_slg_univ_pips.png", width = 7.5, height = 5)

# plot bivariate pip's large
slg_pipb <- read_csv("sim/bsr_lg/pip_biv.csv")

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
ggsave("index/figures/ch4_slg_biv_pips.png", width = 7.5, height = 5)

# plot trivariate pip's large
slg_pipt <- read_csv("sim/bsr_lg/pip_triv.csv")

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
slg_biv <- read_csv("sim/bsr_lg/biv_expresp.csv")

# plot Hg and Ni
slg_biv |> 
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
ggsave("index/figures/ch4_slg_biv_expresp_1.png", width = 6, height = 4)

# plot Cd and As
slg_biv |> 
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
ggsave("index/figures/ch4_slg_biv_expresp_2.png", width = 6, height = 4)

# plot Ni and Co
slg_biv |> 
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
ggsave("index/figures/ch4_slg_biv_expresp_3.png", width = 6, height = 4)

# plot trivariate relationships bsr large
slg_triv <- read_csv("sim/bsr_lg/triv_expresp.csv")
slg_triv <- slg_triv |> 
  mutate(variable1 = case_when(
    j1 == "Hg" ~ "Hg by Ni + Tl", 
    j1 == "Ni" ~ "Ni by Hg + Tl", 
    j1 == "Tl" ~ "Tl by Hg + Ni"), 
    quantile = as.factor(j23quant)) 

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
       color = "Chem 2+3\nquantile") +
  theme(strip.text = element_text(size = 7)) 
ggsave("index/figures/ch4_slg_triv_expresp.png", width = 6, height = 4)


# bsr combine sens --------------------------------------------------------

# univariate pips
bsr_pips <- bind_rows(
  mutate(ssm_pips, size = "Small"), 
  mutate(slg_pips, size = "Large")
) |> 
  group_by(size, case, variable) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) |> 
  mutate(sign = get_sign(case, variable))
write_csv(bsr_pips, "index/data/bsr_pip_sens.csv")

# interactions
cnames <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")

bsr_comb <- bind_rows(
  mutate(ssm_pipb, size = "Small"), 
  mutate(slg_pipb, size = "Large")
) |> 
  mutate(sign = get_sign_bsr(case, Var1) & get_sign_bsr(case, Var2), 
         v1 = cnames[Var1], v2 = cnames[Var2], 
         inter = paste0(v1, "-", v2)) |> 
  mutate(inter2 = ifelse(sign, inter, "none")) |> 
  group_by(size, case, inter2, sign) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n())
write_csv(bsr_comb, "index/data/bsr_int_sens.csv")


# bivar FDR ---------------------------------------------------------------

bsr_comb <- read_csv("index/data/bsr_int_sens.csv")

ksmfdr<- read_csv("sim/bkmr_sm/int_bivarfull.csv")
klgfdr <- read_csv("sim/bkmr_lg/int_bivarfull.csv")

bkmr_comb <- bind_rows(
  mutate(ksmfdr, size = "Small"), 
  mutate(klgfdr, size = "Large")
) |> 
  rowwise() |> 
  mutate(inter1 = paste0(z1, "-", z2), 
         sign = ifelse(
           (inter1 == "Hg-Ni" & case %in% 2:5) | 
             (inter1 == "As-Cd" & case %in% 6:9) | 
             (inter1 == "Co-Ni" & case %in% 10:13), TRUE, FALSE
         ), 
         signif = ifelse(
           between(0, est - 1.96*sd, est + 1.96*sd), FALSE, TRUE
         ), 
         inter2 = ifelse(sign, inter1, "none")) |> 
  group_by(size, case, inter2, sign) |> 
  summarize(sensitivity = sum(signif)/n())

fdr_comb <- bind_rows(
  mutate(bkmr_comb, mod = "BKMR"), 
  mutate(bsr_comb, mod = "BSR")
) |> 
  relocate(mod)

write_csv(fdr_comb, "index/data/comb_int_fdr.csv")

# HgNi only ---------------------------------------------------------------

# univariate significance
hgni_univ <- bind_rows(
  nsm_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Small", mod = "Naive"), 
  nlg_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Large", mod = "Naive"), 
  osm_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Small", mod = "Oracle"), 
  olg_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Large", mod = "Oracle"), 
  ksm_pips |> mutate(met = PIP, size = "Small", mod = "BKMR"), 
  klg_pips |> mutate(met = PIP, size = "Large", mod = "BKMR"), 
  ssm_pips |> mutate(met = PIP, size = "Small", mod = "BSR"), 
  slg_pips |> mutate(met = PIP, size = "Large", mod = "BSR"), 
) |> 
  filter(case %in% 2:5) |> 
  mutate(variable = gsub("[^[:alpha:]]", "", variable), 
         sign = get_sign(case, variable),  
         modsize = paste0(mod, "\n", size))
hgni_univ$modsize <- with(hgni_univ, 
                          factor(modsize, levels = unique(modsize)))
# levels(hgni_univ$modsize)

# mlr plot
hgni1 <- hgni_univ |> 
  filter(mod %in% c("Naive", "Oracle")) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(y = p, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_grid(modsize~case,
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "p-value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")

# bayesian plot
hgni2 <- hgni_univ |> 
  filter(mod %in% c("BKMR", "BSR")) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_grid(modsize~case,
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")

# plot grid
hgnione <- cowplot::plot_grid(
  hgni1 + theme(legend.position = "none"),
  hgni2 + theme(legend.position = "none"),
  labels = "auto", nrow = 2
)
hgnione

hgnilegend <- cowplot::get_legend(hgni1)

# stitch together
cowplot::plot_grid(hgnione, hgnilegend, rel_widths = c(5, 0.6))
ggsave("index/figures/ch4_hgni_univ.png", width = 7.5, height = 7)

# interaction significance
hgni_univ <- bind_rows(
  nsm_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Small", mod = "Naive"), 
  nlg_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Large", mod = "Naive"), 
  osm_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Small", mod = "Oracle"), 
  olg_pvalc |> rename(variable = var) |> 
    mutate(met = p, size = "Large", mod = "Oracle"), 
  ksm_pips |> mutate(met = PIP, size = "Small", mod = "BKMR"), 
  klg_pips |> mutate(met = PIP, size = "Large", mod = "BKMR"), 
  ssm_pips |> mutate(met = PIP, size = "Small", mod = "BSR"), 
  slg_pips |> mutate(met = PIP, size = "Large", mod = "BSR"), 
) |> 
  filter(case %in% 2:5) |> 
  mutate(variable = gsub("[^[:alpha:]]", "", variable), 
         sign = get_sign(case, variable),  
         modsize = paste0(mod, "\n", size))
hgni_univ$modsize <- with(hgni_univ, 
                          factor(modsize, levels = unique(modsize)))
# levels(hgni_univ$modsize)

# mlr plot
hgni1 <- hgni_univ |> 
  filter(mod %in% c("Naive", "Oracle")) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(y = p, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_grid(modsize~case,
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "p-value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")

# bayesian plot
hgni2 <- hgni_univ |> 
  filter(mod %in% c("BKMR", "BSR")) |> 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(y = PIP, color = sign), 
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  size = 0.2) +
  facet_grid(modsize~case,
             labeller = labeller(
               case = as_labeller(appendera, default = label_parsed))) +
  scale_color_manual(values = c("deepskyblue3", "darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "PIP value distribution", 
       color = "Truly\nsignificant", 
       x = "Chemical")

# plot grid
hgnione <- cowplot::plot_grid(
  hgni1 + theme(legend.position = "none"),
  hgni2 + theme(legend.position = "none"),
  labels = "auto", nrow = 2
)
hgnione

hgnilegend <- cowplot::get_legend(hgni1)

# stitch together
cowplot::plot_grid(hgnione, hgnilegend, rel_widths = c(5, 0.6))
ggsave("index/figures/ch4_hgni_univ.png", width = 7.5, height = 7)

# bivariate surfaces

# bkmr small
hgni3 <- ksm_biv |> 
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
  theme(strip.text = element_text(size = 8)) + 
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")

# bkmr large
hgni4 <- klg_biv |> 
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
  theme(strip.text = element_text(size = 8)) + 
  ylim(-6, 6) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")

# bsr small
hgni5 <- ssm_biv |> 
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
  theme(strip.text = element_text(size = 8)) + 
  ylim(-3, 9) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")

# bsr large
hgni6 <- slg_biv |> 
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
  theme(strip.text = element_text(size = 8)) + 
  ylim(-3, 9) +
  labs(x = "Chem 1 value",  
       y = "Estimated response", 
       color = "Chem 2\nquantile")

hgnitwo <- cowplot::plot_grid(
  hgni3 + theme(legend.position = "none"),
  hgni4 + theme(legend.position = "none"),
  hgni5 + theme(legend.position = "none"),
  hgni6 + theme(legend.position = "none"),
  labels = "auto", nrow = 2
)
hgnitwo

hgnilegend2 <- cowplot::get_legend(hgni3 + theme(legend.position = "bottom"))

# stitch together
cowplot::plot_grid(hgnitwo, hgnilegend2, nrow = 2, rel_heights = c(5, 0.4))
# ggsave("temp.png", width = 9, height = 7)
ggsave("index/figures/ch4_hgni_biv.png", width = 9, height = 6)


# run times ---------------------------------------------------------------

## ALSO HAVE TO ADD IN RACE BY ETHNICITY!!! 

t_ksm <- read_rds("sim/bkmr_sm/times.RDS")
t_klg <- read_rds("sim/bkmr_lg/times.RDS")

t_ssm <- read_rds("sim/bsr_sm/times.RDS") |> 
  pivot_longer(cols = c(time_selection, time), 
               names_to = "mod", values_to = "time") |> 
  mutate(mod = ifelse(mod == "time", "BSR mod", "BSR df"))
t_slg <- read_rds("sim/bsr_lg/times.RDS") |> 
  pivot_longer(cols = c(time_selection, time), 
               names_to = "mod", values_to = "time") |> 
  mutate(mod = ifelse(mod == "time", "BSR mod", "BSR df"))

t_nsm1 <- read_rds("sim/_mlr/mlr_mods_sm_times.RDS") |> unlist()
t_nsm <- data.frame(
  case_full = names(t_nsm1), 
  time = t_nsm1
) |> 
  mutate(case = cumsum(case_full != lag(case_full, default = first(case_full))) + 1, 
         time = as.difftime(time, units = "secs"), 
         case = as.character(case)) 

t_nlg1 <- read_rds("sim/_mlr/mlr_mods_lg_times.RDS") |> unlist()
t_nlg <- data.frame(
  case_full = names(t_nlg1), 
  time = t_nlg1
) |> 
  mutate(case = cumsum(case_full != lag(case_full, default = first(case_full))) + 1, 
         time = as.difftime(time, units = "secs"), 
         case = as.character(case)) 

t_osm1 <- read_rds("sim/_oracle/oracle_mods_sm_times.RDS") |> unlist()
t_osm <- data.frame(
  case_full = names(t_osm1), 
  time = t_osm1
) |> 
  mutate(case = cumsum(case_full != lag(case_full, default = first(case_full))) + 1, 
         time = as.difftime(time, units = "secs"), 
         case = as.character(case)) 

t_olg1 <- read_rds("sim/_oracle/oracle_mods_lg_times.RDS") |> unlist()
t_olg <- data.frame(
  case_full = names(t_olg1), 
  time = t_olg1
) |> 
  mutate(case = cumsum(case_full != lag(case_full, default = first(case_full))) + 1, 
         time = as.difftime(time, units = "secs"), 
         case = as.character(case)) 

times <- bind_rows(
  mutate(t_nsm, mod = "Naive", size = "Small"), 
  mutate(t_nlg, mod = "Naive", size = "Large"), 
  mutate(t_osm, mod = "Oracle", size = "Small"), 
  mutate(t_olg, mod = "Oracle", size = "Large"), 
  mutate(t_ksm, mod = "BKMR", size = "Small"), 
  mutate(t_klg, mod = "BKMR", size = "Large"), 
  mutate(t_ssm, size = "Small"), 
  mutate(t_slg, size = "Large"), 
)

# function to format difftime objects with appropriate units
format_difftime <- function(difftime_obj) {
  value <- as.numeric(difftime_obj)
  # determine appropriate units
  if (abs(value) < 60) {
    duration <- as.numeric(difftime_obj, units = "secs")
    result <- sprintf("%.2g %s", duration, "s")
  } else if (abs(value) < 3600) {
    duration <- as.numeric(difftime_obj, units = "mins")
    result <- sprintf("%.2f %s", duration, "m")
  } else if (abs(value) < 86400) {
    duration <- as.numeric(difftime_obj, units = "hours")
    result <- sprintf("%.2f %s", duration, "h")
  } else {
    duration <- as.numeric(difftime_obj, units = "days")
    result <- sprintf("%.2f %s", duration, "d")
  }
  return(result)
}

time_table <- times |> 
  mutate(mod = factor(mod, levels = c("Naive", "Oracle", "BKMR", "BSR df", "BSR mod")), 
         size = factor(size, levels = c("Small", "Large")), 
         case = as.numeric(case), 
         scenario = case_when(
           case == 1 ~ "Base", 
           case %in% c(2, 6, 10, 14) ~ "Mult. Lower", 
           case %in% c(3, 7, 11, 15) ~ "Mult. Higher", 
           case %in% c(4, 8, 12, 16) ~ "Poly. Lower", 
           case %in% c(5, 9, 13, 17) ~ "Poly. Higher"
         ), 
         scenario = factor(scenario, levels = c("Base", "Mult. Lower", "Mult. Higher", 
                                                "Poly. Lower", "Poly. Higher", "Three-way"))) |> 
  group_by(mod, size, scenario) |> 
  summarize(mean_time = mean(time)) |> 
  ungroup() |> 
  rowwise() |> 
  mutate(mean_time = format_difftime(mean_time))
  # mutate(mean_time = format_largest_unit(format_difftime(mean_time)))

time_table_wide <- time_table |> 
  pivot_wider(names_from = scenario, values_from = mean_time)

write_csv(time_table_wide, "sim/tables/time1.csv")
write_csv(time_table_wide, "index/data/time1.csv")


# three-way sens ----------------------------------------------------------

oracle_senst <- bind_rows(
  filter(mutate(osm_comb, size = "Small"), case %in% 14:17), 
  filter(mutate(olg_comb, size = "Large"), case %in% 14:17)
) |> 
  group_by(case, size) |> 
  summarize(sensitivity = sum(p < 0.05)/n()) |> 
  mutate(mod = "Oracle") |> 
  relocate(mod)
bkmr_senst <- bkmr_comb |> 
  filter(case %in% 14:17) |> 
  pivot_wider(names_from = cond, values_from = signif) |> 
  mutate(signif = unname(pick(4) | pick(5) | pick(6))) |> 
  rename(signif = 7) |> 
  group_by(case, size) |> 
  summarize(sensitivity = sum(signif)/n()) |> 
  mutate(mod = "BKMR") |> 
  relocate(mod)
bsr_senst <- bind_rows(
  mutate(ssm_pipt, size = "Small"), 
  mutate(slg_pipt, size = "Large")
) |> 
  group_by(case, size) |> 
  summarize(sensitivity = sum(PIP >= 0.5)/n()) |> 
  mutate(mod = "BSR") |> 
  relocate(mod)

all_senst <- bind_rows(oracle_senst, bkmr_senst, bsr_senst) |> 
  mutate(inter_type = ifelse(case %in% c(14, 15), "Multiplicative", "Polynomial"), 
         effect_size = ifelse(case %in% c(14, 16), "Lower", "Higher"))
write_csv(all_senst, "index/data/triv_sens.csv")
write_csv(all_senst, "sim/tables/triv_sens.csv")








