library(tidyverse)
library(ggforce)

# this code creates visualizations for the final presentation

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

# natural gas emissions ---------------------------------------------------

# read in data 
df <- read_csv("misc_code/cleaned_02-24.csv")

# get natural gas monthly amounts
ng <- df |> 
  select(month, year, cogen = 4, boilers = 5) |> 
  filter(year != 2024) |> 
  rowwise() |> 
  mutate(total_cf = (cogen + boilers)*100, 
         date = as.Date(paste(year, str_pad(month, 2, pad = "0"), "01", sep = "-")))

# make line plot by month
ng |> 
  ggplot(aes(date, total_cf)) +
  geom_line()

# make bar plot by year
ng |> 
  filter(year != 2015) |> 
  group_by(year) |> 
  summarize(total_cf = sum(total_cf)) |> 
  mutate(year = factor(year)) |> 
  ggplot(aes(year, total_cf/1000000)) +
  geom_bar(stat = "identity", fill = "deepskyblue3", alpha = 0.75) +
  # geom_text(aes(label = format(round(total_cf), big.mark = ",")), 
  #           angle = 90, hjust = 1.25) +
  labs(y = "Natural gas combustion\n(Millions of cubic feet)", 
       x = "Calendar year") 

ggsave("misc_code/naturalgas-byyear.png", width = 6, height = 4)


# natural gas components --------------------------------------------------

comps <- readxl::read_xlsx("misc_code/naturalgas_emissions.xlsx", sheet = 1)
metals <- readxl::read_xlsx("misc_code/naturalgas_emissions.xlsx", sheet = 2)

all_poll <- bind_rows(
  mutate(comps, type = "Compounds"), 
  mutate(metals[,2:4], type = "Metals")
)

ng2023 <- ng |> 
  group_by(year) |> 
  summarize(total_cf = sum(total_cf)) |> 
  filter(year == 2023) |> 
  select(-year) |> 
  as.numeric()

all_poll2023 <- all_poll |> 
  janitor::clean_names() |> 
  mutate(emission_factor_lb_106_scf = 
           as.numeric(str_replace_all(emission_factor_lb_106_scf, 
                                      "[,<]", "")), 
         emissions_lb = emission_factor_lb_106_scf * ng2023 / 106)

metals2023 <- all_poll2023 |> 
  filter(type == "Metals")

radius <- sqrt(metals2023$emissions_lb/pi)
set.seed(7)
packing_result <- packcircles::circleLayout(radius)[[1]]

pmetals2023 <- bind_cols(metals2023, packing_result)
# pmetals2023$abbrev <- c("As", "Ba", "Be", "Cd", "Cr",
#                         "Co", "Cu", "Pb", "Mn", "Hg",
#                         "Mo", "Ni", "Se", "V", "Zn")

pmetals2023 |> 
  mutate(pollutant = ifelse(pollutant == "Zinc", 
                            paste0(pollutant, "\n(", 
                                   format(round(emissions_lb), big.mark = ","), 
                                   " lbs)"), 
                            pollutant)) |> 
  ggplot(aes(x = x, y = y, label = pollutant)) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = emission_factor_rating), 
              color = NA) + 
  geom_text(aes(size = radius), show.legend = F) +
  coord_fixed() +
  theme_void() +
  scale_size_continuous(range = c(3, 7)) +
  scale_fill_manual(values = c("#FFFF66", "#FF9966", "#FF6666")) +
  labs(fill = "Grade (A-E)")
ggsave("misc_code/emission_amts.png", width = 6, height = 6)


# slab and spike prior ----------------------------------------------------

ggplot(data.frame(x = c(2, 2)), aes(x = x)) +
  stat_function(fun = ~1/.x, n=1000) + #, xlim = c(0.02, 2)
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 25), 
               linewidth = 4, alpha = 0.5, color = "deepskyblue3") + 
  ylim(0, 50) +
  labs(x = "x",
       y = "f(x)") 
ggsave("misc_code/slabandspike.png", width = 4, height = 2.5)



# smaller hg-ni -----------------------------------------------------------

# bkmr small
ksm_biv <- read_csv("sim/bkmr_sm/biv_expresp.csv") |> 
  mutate(mod = "BKMR", size = "Small") |> 
  filter(case == 3)
# bkmr large
klg_biv <- read_csv("sim/bkmr_lg/biv_expresp.csv") |> 
  mutate(mod = "BKMR", size = "Large") |> 
  filter(case == 3)
# bsr small
ssm_biv <- read_csv("sim/bsr_sm/biv_expresp.csv") |> 
  select(variable1 = j1, variable2 = j2, z1 = j1val, quantile = j2quant, est, trial, case) |> 
  filter(case == 3) |> 
  mutate(mod = "BSR", size = "Small")
# bsr large
slg_biv <- read_csv("sim/bsr_lg/biv_expresp.csv") |> 
  select(variable1 = j1, variable2 = j2, z1 = j1val, quantile = j2quant, est, trial, case) |> 
  filter(case == 3) |> 
  mutate(mod = "BSR", size = "Large")

# combine
hgni7 <- bind_rows(ksm_biv, klg_biv, ssm_biv, slg_biv)

# labels 
label <- data.frame(
  mod = c("BKMR", "BKMR", "BSR", "BSR"), 
  size = factor(c("Small", "Large", "Small", "Large"), 
                levels = c("Small", "Large")), 
  z1 = rep(-1.5, 4), 
  est = rep(7.5, 4), 
  val = c(0.11, 0.88, 0.10, 0.77)
)

# bivariate surfaces
hgni7 |>
  filter(variable1 == "Hg") |>
  mutate(quantile = as.factor(quantile), 
         size = factor(size, levels = c("Small", "Large"))) |>
  ggplot(aes(z1, est)) +
  geom_line(aes(color = quantile,group = interaction(trial, quantile)), 
            alpha = 0.2) +
  geom_text(data = label, mapping = aes(label = val)) +
  ggh4x::facet_grid2(size ~ mod) +
  scale_color_manual(
    values = c("deepskyblue3", "palegreen3", "darkorange"),
    guide = guide_legend(override.aes = list(alpha = 1),
                         reverse = TRUE)
  ) +
  xlim(-2, 2.5) +
  ylim(-5, 9) +
  labs(x = "Hg value",
       y = "Estimated response",
       color = "Ni quantile")

ggsave("misc_code/hgni7.png", width = 6, height = 5)

# re relationship -------------------------------------------------------

getsign_re <- function(chem) {
  return(chem %in% c("Hg", "Ni", "Sb", "Sn"))
}

which_race <- function(race) {
  case_when(
    race == 1 ~ "Non-Hisp. white", 
    race == 2 ~ "Non-Hisp. black", 
    race == 3 ~ "Non-Hisp. other", 
    race == 4 ~ "Hispanic born in US", 
    race == 5 ~ "Hispanic born outside US",
    race == 6 ~ "Collapsed non-Hisp."
  )
}

size_race <- function(race) {
  case_when(
    race == 1 ~ "16", 
    race == 2 ~ "27", 
    race == 3 ~ "13", 
    race == 4 ~ "87", 
    race == 5 ~ "109",
    race == 6 ~ "56"
  )
}

#### BKMR #####

# small
ksmre_expresp <- read_csv("sim/re/ksmre_expresp.csv")
# large
klgre_expresp <- read_csv("sim/re/klgre_expresp.csv")
# combine
bkmr_re_expresp <- bind_rows(
  mutate(ksmre_expresp, size = "Small"), 
  mutate(klgre_expresp, size = "Large")
) |> 
  filter(race != 6) 

# second two cases
kre34 <- bkmr_re_expresp |>
  filter(case == 4) |>
  mutate(
    which_race = which_race(race),
    size_race = size_race(race),
    size = factor(size, levels = c("Small", "Large"))
  ) |>
  mutate(race2 = factor(race,
                        levels = c(1, 2, 3, 4, 5),
                        labels = c("Non-Hisp. white", "Non-Hisp. other", "Non-Hisp. black", 
                                   "Hispanic born in US", "Hispanic born outside US")),
         case = factor(ifelse(case == 3, "Lower", "Higher"),
                       levels = c("Lower", "Higher")))
#### BSR #####

# small
ssmre_expresp <- read_csv("sim/re/ssmre_expresp.csv")
# large
slgre_expresp <- read_csv("sim/re/slgre_expresp.csv")
# combine
bsr_re_expresp <- bind_rows(
  mutate(ssmre_expresp, size = "Small"), 
  mutate(slgre_expresp, size = "Large")
) |> 
  filter(race != 6) 

# second two cases
sre34 <- bsr_re_expresp |> 
  filter(case == 3) |> 
  mutate(which_race = which_race(race), 
         size_race = size_race(race), 
         size = factor(size, levels = c("Small", "Large")))  |>
  mutate(
    race2 = factor(
      race,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Non-Hisp. white", "Non-Hisp. other", "Non-Hisp. black", 
                 "Hispanic born in US", "Hispanic born outside US")
    ), 
    case = factor(ifelse(case == 3, "Lower", "Higher"), 
                  levels = c("Lower", "Higher"))
  ) |> 
  rename(z1 = j1val) |> 
  select(-lower, -upper)

#### all together ####

# bind bsr and bkmr together
all34 <- bind_rows(
  mutate(kre34, mod = "BKMR"), 
  mutate(sre34, mod = "BSR")
)

label2 <- data.frame(
  mod = c("BKMR", "BKMR", "BSR", "BSR"), 
  size = factor(c("Small", "Large", "Small", "Large"), 
                levels = c("Small", "Large")), 
  z1 = rep(-1.5, 4), 
  est = rep(8, 4), 
  val = c("0.02", "0.21", "", "")
)

all34 |> 
  ggplot(aes(z1, est)) +
  geom_line(aes(color = race2, group = interaction(trial, race2)), alpha = 0.15) +
  geom_text(data = label2, aes(label = val)) +
  ylim(-6, 10) + 
  xlim(-2, 2) +
  facet_grid(size ~ mod) + #breaks = c(1, 2, 3, 4, 5)
  scale_color_manual(values = rev(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"))) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Hg", y = "Estimate", color = "Race")

ggsave("misc_code/hgrace.png", width = 6, height = 4)


# list of chems -----------------------------------------------------------

chemlist <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")

chempairs <- t(combn(chemlist, 2))
pairlist <- paste(chempairs[, 1], chempairs[, 2], sep = "-")

chemtrios <- t(combn(chemlist, 3))
triolist <- paste(chemtrios[, 1], chemtrios[, 2], chemtrios[, 3], sep = "-")

sink("misc_code/chemlist.txt")
cat(chemlist, sep = ", ")
cat(pairlist, sep = ", ")
cat(triolist, sep = ", ")
sink()


# toy code ----------------------------------------------------------------

# generate data from distribution
set.seed(0) # reproducibility
x <- seq(0, 25, length.out = 51)
Y <- exp(x/10) + 2*sin(x/2) + rnorm(51, mean = 0, sd = 0.5)
df <- data.frame(x, Y)

for(mean.val in seq(0, 25, by = 1.25)) {
  # get normal distribution of weights around query points
  df$Weight <- dnorm(df$x, mean = mean.val, sd = 1)
  
  # plot points colored by their weights
  p1 <- ggplot(df, aes(x, Y)) +
    geom_point(aes(color = Weight)) +
    geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                  linetype = "dashed", color = "darkorange") + 
    geom_vline(xintercept = mean.val, linetype = "dotted") +
    theme(legend.position = "none") +
    scale_color_gradient(low = "lightblue", high = "black")
  
  # plot a curve of weights
  normcurv <- data.frame(x = seq(0, 25, length.out = 251)) 
  normcurv$Weight <- dnorm(normcurv$x, mean = mean.val, sd = 1)
  p2 <- ggplot(normcurv, aes(x, Weight, color = Weight)) +
    geom_line() +
    scale_y_continuous(breaks = c(0, 0.2, 0.4)) +
    scale_color_gradient(low = "lightblue", high = "black") +
    theme(legend.position = "none", 
          plot.margin = unit(c(5.5, 5.5, 5.5, 10), "points")) 
  
  # stitch plots together
  q2 <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.7, 0.3))
  
  ggsave(paste0("misc_code/kernel/plot", 
                str_pad(sprintf("%.2f", mean.val), 5, pad = "0"), ".png"), 
         plot = q2, width = 5, height = 4)
}



# save plot

