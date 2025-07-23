# (C) Irene Sophia Plank
# 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = c("/Users/vilya/Documents/MSYNC/data/MSYNC_CT_ET_FACE-MAPPER", 
              "/Users/vilya/Documents/MSYNC/data")
} else {
  dt.path = c("/media/emba/emba-2/MSYNC/data/MSYNC_CT_ET_FACE-MAPPER", 
              "/media/emba/emba-2/MSYNC/data")
}

# select time window length in seconds
tiwo = 600

# Read in data ------------------------------------------------------------

# read in the sections as a dataframe
df.sec = read_csv(file.path(dt.path[2], "MSYNC_CT_events.csv")) %>%
  pivot_wider(names_from = name, values_from = `timestamp [ns]`) %>%
  mutate(
    # shift the end trigger to keep specific tiwo
    end = clap + 1000000000 * tiwo,
    # if the recording is too short, then use the full recording
    end = if_else(end > recording.end, recording.end, end)
  ) %>% select(-recording.end)

# read in the data, merge with sections and remove everything outside
df = merge(read_csv(file.path(dt.path[1], "fixations_on_face.csv")), 
           df.sec) %>%
  # now we finally filter out the fixations outside the section
  filter(
    `start timestamp [ns]` >= clap & `end timestamp [ns]` <= end
  ) %>%
  # now we compute the duration of the fixations in ms
  mutate(
    dur.ms = (`end timestamp [ns]` - `start timestamp [ns]`) / 1000000
  )

# Looking at the face of the other one ------------------------------------

# now we can aggregate the durations and compute dwell times for the individual
df.agg = df %>%
  mutate(dyad = substr(subID, 1, 8)) %>%
  group_by(subID, dyad) %>%
  mutate(
    total.dur.ms = sum(dur.ms),
    ROI = if_else(`fixation on face`, 'face', 'other')
  ) %>%
  group_by(subID, dyad, ROI, total.dur.ms) %>%
  summarise(
    dur.ms = sum(dur.ms)
  ) %>%
  mutate(
    dwell = dur.ms / total.dur.ms
  )

# save data frame
write_csv(df.agg, file.path(dt.path[2], "MSYNC_ET_CT_indi.csv"))


# Shared face contact -----------------------------------------------------

# when were two people looking at each other at the same time
df.face = df %>%
  filter(`fixation on face`) %>%
  mutate(
    start.aligned = `start timestamp [ns]` - clap,
    end.aligned = `end timestamp [ns]` - clap,
    start.ms = round(start.aligned/1000000),
    end.ms   = round(end.aligned/1000000),
    dyad     = substr(subID, 1, 8),
    position = substr(subID, nchar(subID), nchar(subID))
  )

# initialise dataframe
df.dyad = data.frame()

# loop through dyad
for (d in unique(df.face$dyad)) {
  # get the positions of the associated interaction partners
  pos = unique(df.face %>% filter(dyad == d) %>% select(position))$position
  if (length(pos) != 2) {
    df.dyad = rbind(
      df.dyad, 
      data.frame(dyad = d, shared.ms = NA, 
                 comment = 'Not two interaction partners')
    )
    next
  } else {
    comment = ''
  }
  # initialise empty list for the vectors
  t = list()
  # loop through the positions
  for (po in pos) {
    # focus on this one person
    df.sel = df.face %>% filter(position == po & dyad == d)
    if (nrow(df.sel) == 0) {
      comment = paste0(
        comment,  
        sprintf("No relevant fixations for %s; ", po))
      next
    }
    # initialise vector
    x = c()
    # loop through fixations
    for (i in 1:nrow(df.sel)) {
      x = c(x, df.sel$start.ms[i]:df.sel$end.ms[i])
    }
    # add this info to the list
    t[[po]] = x
  }
  # calculate overlap
  if (length(t) != 2) {
    shared.ms = NA
  } else {
    shared.ms = length(intersect(t$L, t$R))
  }
  df.dyad = rbind(
    df.dyad, 
    data.frame(dyad = d, shared.ms, comment)
  )
  comment = ''
}

# add mean of total fixation duration for each dyad
df.dyad = merge(df.dyad, 
                df.agg %>% group_by(dyad) %>% 
                  summarise(total.dur.ms.mean = mean(total.dur.ms))) %>%
  mutate(
    shared.perc = shared.ms/total.dur.ms.mean
  )

write_csv(df.dyad, file.path(dt.path[2], "MSYNC_ET_CT_dyad.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "ET_CT.Rdata"))

