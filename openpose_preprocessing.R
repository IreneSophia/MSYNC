# Script inspired by Fujiwara et al. (2022), shared at https://osf.io/2t6a5/ 
# This script was written by I. S. Plank

library(zoo)
library(tidyverse)

# clean workspace
rm(list = ls())

# set paths and input files
dt.path = '/Users/vilya/Documents/MSYNC/data/preprocessedOP'
files = list.files(path = dt.path, pattern = "*.csv")

# settings for the videos and durations
fps      = 120
skip     = 10
duration = 600

# initialise the data frames
df.ref = data.frame()
df     = data.frame()

for (f in files){
  
  # load the dyad's data
  tmp = read_csv(file.path(dt.path, f)) %>%
    pivot_longer(cols = where(is.numeric)) %>%
    mutate(
      dyad  = gsub("(.+)_CT.*", "\\1", f),
      side  = substr(name, 1, 1), 
      key   = as.numeric(substr(name, 2, nchar(name)-1)),
      name  = substr(name, nchar(name), nchar(name))
    ) %>% 
    # focus on the relevant keypoints
    filter(key <= 8) %>%
    pivot_wider(id_cols = c(dyad, side, key, frame)) %>%
    # set all frames to NA where the confidence is below 2/3
    mutate(
      x = if_else(c < 2/3, NA, x),
      y = if_else(c < 2/3, NA, y)
    )
  
  if (sum(is.na(tmp$x))/nrow(tmp) >= 1/3) {
    warning(sprintf("Dyad has %.1f%% missing data.", 
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
                       0 ~ "head",
                       1 ~ "neck",
                       2 ~ "shoulderL",
                       3 ~ "ellbowL",
                       4 ~ "handL", 
                       5 ~ "shoulderR", 
                       6 ~ "ellbowR", 
                       7 ~ "handR", 
                       8 ~ "ref")
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
      select(side, frame, x.ref, y.ref), 
    tmp) %>%
    mutate(
      x = x.ma - x.ref,
      y = y.ma - y.ref
    ) %>% filter(key != "ref") %>%
    select(dyad, side, key, frame, x, y) %>%
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


# save the data
saveRDS(df.ref, file = file.path(dt.path, "MSYNC_OP_ref.rds"))
saveRDS(df, file = file.path(dt.path, "MSYNC_OP.rds"))
