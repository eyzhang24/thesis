library(tidyverse)
library(NLinteraction)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

#extract names of files
main_folder <- "sim/bsr_df_sm"
subfolders <- list.dirs(main_folder, full.names = TRUE, recursive = TRUE)
mod_subfolders <- subfolders[grepl("mods", subfolders)]
mod_labels <- gsub("\\D", "", mod_subfolders)
mod_labels <- ifelse(mod_labels == "", 1, as.numeric(mod_labels))

mod_paths <- mod_subfolders |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = mod_labels)

bsr1 <- read_rds("sim/bsr_df_sm/_rslurm_bsr_sm/mods/bsr_sm__base_1.RDS")
bsr1_1 <- bsr1[[1]]

waic <- names(mod_paths) |> 
  purrr::map_df(\(x) {
    mod_paths[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- data.frame(
          df = c(1, 2, 3, 4, 5), 
          waic = c(bsr[[1]]$waic, 
                   bsr[[2]]$waic, 
                   bsr[[3]]$waic, 
                   bsr[[4]]$waic, 
                   bsr[[5]]$waic), 
          trial = rep(substr(y, nchar(y) - 4, nchar(y) - 4))
        )
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })

write_csv(waic, "sim/tables/test_waic.csv")

##look at waic values
waic <- read_csv("sim/tables/test_waic.csv")

min_waic <- waic |> 
  group_by(case, trial) |> 
  filter(waic == min(waic)) |> 
  arrange(as.numeric(case))

waic |> 
  filter(trial <= 5) |> 
  mutate(df = as.factor(df)) |> 
  ggplot(aes(x = "", y = waic, color = df)) +
  geom_point() +
  facet_wrap(~interaction(str_pad(case, 2, pad = "0"), trial), 
             scales = "free_y", ncol = 17) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  )

ggsave("sim/figs/test_waic.png", height = 7.5, width = 10)

#second trial
#extract names of files
main_folder2 <- "sim/bsr_df_sm2"
subfolders2 <- list.dirs(main_folder2, full.names = TRUE, recursive = TRUE)
mod_subfolders2 <- subfolders2[grepl("mods", subfolders2)]
mod_labels2 <- gsub("\\D", "", substr(mod_subfolders2, 16, nchar(mod_subfolders2)))
mod_labels2 <- ifelse(mod_labels2 == "", 1, as.numeric(mod_labels2))

mod_paths2 <- mod_subfolders2 |> 
  purrr::map(\(x) {
    list.files(x, full.names = TRUE)
  }) |> 
  setNames(nm = mod_labels2)

waic2 <- names(mod_paths2) |> 
  purrr::map_df(\(x) {
    mod_paths2[[x]] |> 
      purrr::map_df(\(y) {
        bsr <- read_rds(y)
        result <- data.frame(
          df = c(1, 2, 3, 4, 5), 
          waic = c(bsr[[1]]$waic, 
                   bsr[[2]]$waic, 
                   bsr[[3]]$waic, 
                   bsr[[4]]$waic, 
                   bsr[[5]]$waic), 
          trial = rep(substr(y, nchar(y) - 4, nchar(y) - 4))
        )
        rm(bsr)
        return(result)
      }) |> 
      mutate(case = x)
  })

write_csv(waic2, "sim/tables/test_waic2.csv")

####compare them
waic <- read_csv("sim/tables/test_waic.csv")
waic2 <- read_csv("sim/tables/test_waic2.csv")
waic_comb <- waic |> 
  mutate(iter = 1) |> 
  bind_rows(mutate(waic2, iter = 2))

waic_comb |> 
  filter(trial <= 5) |> 
  mutate(df = as.factor(df), 
         iter = factor(ifelse(iter == 1, "F", "P"), levels = c("P", "F"))) |> 
  ggplot(aes(x = iter, y = waic, color = df)) +
  geom_point() +
  geom_line(aes(group = df)) + 
  ggh4x::facet_grid2(paste0("Trial ", trial) ~ case, 
                     scales = "free_y", independent = "y") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        plot.caption = element_text(hjust = 0, size = 7)) +
  labs(y = "WAIC", x = "# MCMC iterations (F = 50,000, P = 5,000)", 
       color = "Degrees\nfreedom",
       caption = paste0(
         "Scenarios are labelled in the top strip as follows:\n", 
         "1 = base case; 2 = HgNi mult. small; 3 = HgNi mult. large; ", 
         "4 = HgNi poly. small; 5 = HgNi poly. large; ", 
         "6 = CdAs mult. small; 7 = CdAs mult. large;", 
         "8 = CdAs poly. small; 9 = CdAs poly. large;\n", 
         "10 = NiCo mult. small; 11 = NiCo mult. large; ", 
         "12 = NiCo poly. small; 13 = NiCo poly. large; ", 
         "14 = HgNiTl mult. small; 15 = HgNiTl mult. large; ", 
         "16 = HgNiTl poly. small; 17 = HgNiTl poly. large"))

ggsave("index/figures/test_waic2.png", height = 6, width = 9)

waic_min <- waic_comb |> 
  mutate(iter = ifelse(iter == 1, "full", "partial")) |> 
  arrange(iter, case, trial) |> 
  filter(trial <= 5) |> 
  group_by(iter, trial, case) |> 
  filter(waic == min(waic)) |> 
  ungroup() |> 
  pivot_wider(id_cols = c(trial, case), 
              names_from = iter, values_from = df) |> 
  mutate(equal = (full == partial))
mean(waic_min$equal)  
