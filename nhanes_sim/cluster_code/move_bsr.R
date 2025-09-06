# move bsr
library(tidyverse)

# # check progress on bkmr
# list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
# list_bkmr <- list_files[grepl("_rslurm_klgre(\\d+)", list_files)]
# 
# bkmr_subf <- list.dirs(list_bkmr, full.names = TRUE, recursive = TRUE)
# bkmr_mod <- bkmr_subf[grepl("mods", bkmr_subf)]
# for (f in bkmr_mod) {
#   message(f)
#   x <- list.files(f[1])
#   message("   ",  x[1])
#   message("   ",  x[length(x)])
#   
# }
# bkmr_time <- bkmr_subf[grepl("times", bkmr_subf)]

# extracting file names

list_files <- list.dirs(".", full.names = FALSE, recursive = FALSE)
list_bsri <- list_files[grepl("_rslurm_ssmre(\\d+)_", list_files)]

bsri_subf <- list.dirs(list_bsri, full.names = TRUE, recursive = TRUE)
bsri_mod <- bsri_subf[grepl("mods", bsri_subf)]
bsri_time <- bsri_subf[grepl("times", bsri_subf)]

# check number of mods w/ two files
result_df <- data.frame(
  folder_name = character(),
  file_names = I(list()),  # Store file names as a list column
  num_files = integer(),
  stringsAsFactors = FALSE
)

# Iterate over each folder name in the vector
for (f in bsri_mod) {
  # Get all file names in the current folder
  files <- list.files(f)
  
  # Append the information to the data frame
  result_df <- rbind(
    result_df,
    data.frame(
      folder_name = f,
      file_names = I(list(files)),  # Store as a list
      num_files = length(files),
      stringsAsFactors = FALSE
    )
  )
}

print("stop")

# first collapse any mod folders with multiple rds's
for (f in bsri_mod) {
  files <- list.files(f, full.names = T) 
  
  if (length(files) == 2) {
    print(f)
    # Read in the .RDS files
    list1 <- readRDS(files[1])
    list2 <- readRDS(files[2])
    
    # Check if both lists are of size 6
    if (length(list1) != 6 || length(list2) != 6) {
      stop(paste("Error: Lists in folder", f, "are not of size 6"))
    }
    
    # Initialize a combined list
    combined_list <- vector("list", 6)
    
    # Iterate through the indices of the lists
    for (i in seq_along(combined_list)) {
      # Get elements at the current index
      elem1 <- list1[[i]]
      elem2 <- list2[[i]]
      
      # Apply the combination logic
      if (is.null(elem1) && is.null(elem2)) {
        combined_list[[i]] <- NULL
      } else if (is.null(elem1)) {
        combined_list[[i]] <- elem2
      } else if (is.null(elem2)) {
        combined_list[[i]] <- elem1
      } else {
        # Stop if both elements are not NULL
        stop(paste("Conflict in folder", f, "at index", i))
      }
    }
    
    saveRDS(combined_list, 
            paste0(sub("(df\\d).*", "\\1", files[1]), ".RDS"))
  }
}

# copy mods to correct folder and delete original folder
# f <- bsri_mod[[1]]

for (f in bsri_mod[289:400]) {
  file_list <- list.files(f, full.names = TRUE)
  # remove anything with pone or ptwo in the name
  filtered <- file_list[!grepl("pone|ptwo", file_list)]
  if (length(filtered) != 1) {
    stop(paste("Folder", f, "does not contain exactly one file after filtering."))
  }
  file.copy(from = filtered,
            to   = paste0(substr(f, 1, 15), "/mods"))
}

# 108, 136, 224, 288 needs to be fixed, first four are null

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

list_files <- list.dirs(".", full.names = F, recursive = FALSE)
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
            to   = paste0(substr(f, 6, 20), "/times"))
}

# move files into temporary folder
for (f in list_bsri) {
  file.copy(from = f, 
            to = "temp", 
            recursive = T)
}

# delete old copy
for (f in list_bsri) {
  unlink(f, recursive = TRUE)
}

## FIX THE TIME THING

list_files <- list.dirs("temp", full.names = T, recursive = FALSE)
list_bsri <- list_files[grepl("_rslurm_slgre(\\d+)_", list_files)]

bsri_time <- bsri_subf[grepl("times", bsri_subf)]

bsri_timepath <- bsri_time |> 
  purrr::map(\(x) {
    temp <- list.files(x, full.names = F)
    temp <- temp[which(!grepl("df", temp))]
    return(temp)
  }) |> 
  unlist()

for(i in 1:400) {
  from.path <- paste0(list_bsri[i], '/results_0.RDS')
  to.path <- paste0(list_bsri[i], "/", bsri_timepath[i])
  file.rename(from.path, to.path)
}

for(i in 1:400) {
  from.path <- paste0(list_bsri[i], "/", bsri_timepath[i])
  to.path <- paste0(substr(list_bsri[i], 6, 20), "/times")
  file.copy(from.path, to.path, recursive = T)
}



