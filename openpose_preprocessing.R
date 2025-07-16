# Script inspired by Fujiwara et al. (2022), shared at https://osf.io/2t6a5/ 
# This script was written by I. S. Plank

library(zoo)
library(tidyverse)

# clean workspace
rm(list = ls())

# set the task
task = "MG"

# set path
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = "/Users/vilya/Documents/MSYNC/data/preprocessedOP"
} else {
  dt.path = "/media/emba/emba-2/MSYNC/data/preprocessedOP"
}

# set paths and input files
files = list.files(path = dt.path, pattern = sprintf(".*%s.*.csv", task))

# settings for the videos and durations
fps      = 120
skip     = 10
duration = 600 # only applies for CT, but for MG we use all of it anyway

# initialise the data frames
df.ref = data.frame()
df     = data.frame()

# total number of files
total = length(files)
x = 0

for (f in files){
  
  if (task == "CT") {
    max.key = 8
  } else {
    max.key = 14
  }
  
  x = x + 1
  print(sprintf("%s: %i of %i", Sys.time(), x, total))
  
  # load the dyad's data
  tmp = read_csv(file.path(dt.path, f), show_col_types = F) %>%
    rename_with(~ "frame", .cols = where(is.character)) %>%
    pivot_longer(cols = where(is.numeric)) %>%
    mutate(
      frame = as.numeric(frame),
      dyad  = gsub(sprintf("(.+)_%s.*", task), "\\1", f),
      side  = substr(name, 1, 1), 
      key   = as.numeric(substr(name, 2, nchar(name)-1)),
      name  = substr(name, nchar(name), nchar(name))
    ) %>% 
    # focus on the relevant keypoints
    filter(key <= max.key) %>%
    pivot_wider(id_cols = c(dyad, side, key, frame)) %>%
    # set all frames to NA where the confidence is below 2/3
    mutate(
      x = if_else(c < 2/3, NA, x),
      y = if_else(c < 2/3, NA, y)
    )
  
  # add phase information if MG
  if (task == "MG") {
    tmp$phase = gsub(".*_SC_p(.+)_cut.*", "\\1", f)
  } else {
    tmp$phase = "CT"
  }
  
  if (sum(is.na(tmp$x))/nrow(tmp) >= 1/3) {
    warning(sprintf("Dyad %s has %.1f%% missing data.", f,
                    100*sum(is.na(tmp$x))/nrow(tmp)))
  }
  
  tmp = tmp %>% 
    arrange(side, key, frame) %>%
    group_by(side, key) %>%
    mutate(
      # interpolate the missing values with interpolation
      x.ip = na.approx(x, na.rm = F),
      y.ip = na.approx(y, na.rm = F),
      # compute a rolling mean
      x.ma = rollmean(x.ip, k = fps+1, fill = NA),
      y.ma = rollmean(y.ip, k = fps+1, fill = NA),
      # convert frames to numbers
      frame = as.numeric(frame),
      # convert keys to relevant names
      key = case_match(key, 
                       0  ~ "head",
                       1  ~ "neck",
                       2  ~ "shoulderL",
                       3  ~ "ellbowL",
                       4  ~ "handL", 
                       5  ~ "shoulderR", 
                       6  ~ "ellbowR", 
                       7  ~ "handR", 
                       8  ~ "ref",
                       9  ~ "hipL",
                       10 ~ "kneeL",
                       11 ~ "footL",
                       12 ~ "hipR",
                       13 ~ "kneeR",
                       14 ~ "footR")
    ) %>%
    filter(frame >= skip*fps & frame < (duration*fps + skip*fps)) %>%
    # add the movement quantity
    group_by(side, key) %>%
    arrange(side, key, frame) %>%
    mutate(
      # difference in pixels
      x.diff = lag(x.ma) - x.ma,
      y.diff = lag(y.ma) - y.ma,
      #  Euclidean distance in pixels
      dist   = sqrt(x.diff**2 + y.diff**2)
    ) %>% ungroup()
  
  # use keypoint 8 as reference for postural mirroring
  tmp.ref = merge(
    tmp %>% ungroup() %>%
      filter(key == "ref") %>%
      mutate(
        x.ref = x.ma, 
        y.ref = y.ma
      ) %>%
      select(side, phase, frame, x.ref, y.ref), 
    tmp) %>%
    mutate(
      x = x.ma - x.ref,
      y = y.ma - y.ref
    ) %>% filter(key != "ref") %>%
    select(dyad, phase, side, key, frame, x, y) %>%
    # if on the right side, then flip it
    mutate(
      x = if_else(side == "R", x * (-1), x)
    ) %>%
    pivot_longer(cols = c(x, y), names_to = "axis") %>%
    pivot_wider(names_from = side, values_from = value) %>%
    arrange(key, axis, frame)
  
  df     = rbind(df, tmp)
  df.ref = rbind(df.ref, tmp.ref)
  
}

# max data points
if (task == "CT") {
  max.dt = duration*fps
} else {
  max.dt = 3*60*fps
}

# check how many dyads with too few valid data points
df.ref %>% group_by(dyad, phase, key, axis) %>%
  summarise(
    valid = min(sum(!is.na(L))/max.dt, sum(!is.na(R))/max.dt)
  ) %>% filter(valid <= 2/3) %>% arrange(valid)

# need to exclude MSYNC_15 for CT
if (task == "CT") {
  df     = df %>% filter(dyad != "MSYNC_15")
  df.ref = df.ref %>% filter(dyad != "MSYNC_15")
}  

# save the data
saveRDS(df.ref, file = file.path(dt.path, sprintf("MSYNC_OP_%s_ref.rds", task)))
saveRDS(df, file = file.path(dt.path, sprintf("MSYNC_OP_%s.rds", task)))
