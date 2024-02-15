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

  