# move bsr
library(tidyverse)

# one <- read_csv("sim/bkmr_sm/int_bivar_fullhgni.csv")
# two <- read_csv("sim/bkmr_sm/int_bivar_fullrest.csv")
# thee <- bind_rows(one, two)
# write_csv(thee, "sim/bkmr_sm/int_bivarfull.csv")

# check progress on bkmr
list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_bkmr <- list_files[grepl("_rslurm_klgre(\\d+)", list_files)]

bkmr_subf <- list.dirs(list_bkmr, full.names = TRUE, recursive = TRUE)
bkmr_mod <- bkmr_subf[grepl("mods", bkmr_subf)]
for (f in bkmr_mod) {
  message(f)
  x <- list.files(f[1])
  message("   ",  x[1])
  message("   ",  x[length(x)])
  
}
# bkmr_time <- bkmr_subf[grepl("times", bkmr_subf)]

# extracting file names

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_bsri <- list_files[grepl("_rslurm_ssmre(\\d+)_", list_files)]

bsri_subf <- list.dirs(list_bsri, full.names = TRUE, recursive = TRUE)
bsri_mod <- bsri_subf[grepl("mods", bsri_subf)]
bsri_time <- bsri_subf[grepl("times", bsri_subf)]

# copy mods to correct folder and delete original folder
f <- bsri_mod[[1]]

for (f in bsri_mod) {
  file.copy(from = list.files(f, full.names = TRUE),
            to   = paste0(substr(f, 1, 15), "/mods"))
}

# bsri_mod2 <- bsri_mod[-100]
# 
# for (f in bsri_mod2) {
#   unlink(substr(f, 1, 19), recursive = TRUE)
# }

# copy times to correct folder and delete original folder
## I MESSED THIS UP
## MIGHT HAVE TO GO BACK AND RE-RUN IT!!

for (f in bsri_time) {
  file.copy(from = list.files(f, full.names = TRUE),
            to   = paste0(substr(f, 1, 15), "/times"))
}

for (f in bsri_mod) {
  unlink(substr(f, 1, 19), recursive = TRUE)
}

## LARGE

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_bsri <- list_files[grepl("_rslurm_slgre(\\d+)_", list_files)]

bsri_subf <- list.dirs(list_bsri, full.names = TRUE, recursive = TRUE)
bsri_mod <- bsri_subf[grepl("mods", bsri_subf)]
bsri_time <- bsri_subf[grepl("times", bsri_subf)]

# copy mods to correct folder and delete original folder
f <- bsri_mod[[1]]

for (f in bsri_mod) {
  file.copy(from = list.files(f, full.names = TRUE),
            to   = paste0(substr(f, 1, 15), "/mods"))
}

# bsri_mod2 <- bsri_mod[-100]
# 
# for (f in bsri_mod2) {
#   unlink(substr(f, 1, 19), recursive = TRUE)
# }

# copy times to correct folder and delete original folder

for (f in bsri_time) {
  file.copy(from = list.files(f, full.names = TRUE),
            to   = paste0(substr(f, 1, 15), "/times"))
}

# delete files
for (f in bsri_mod) {
  unlink(substr(f, 1, 19), recursive = TRUE)
}


