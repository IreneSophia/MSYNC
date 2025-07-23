# (C) Irene Sophia Plank
# 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = c("/Users/vilya/Documents/MSYNC/data/MSYNC_MG_ET_FACE-MAPPER", 
              "/Users/vilya/Documents/MSYNC/data")
} else {
  dt.path = c("/media/emba/emba-2/MSYNC/data/MSYNC_MG_ET_FACE-MAPPER", 
              "/media/emba/emba-2/MSYNC/data")
}

# Read in data ------------------------------------------------------------

# read in the sections as a dataframe
df.sec = read_csv(file.path(dt.path[2], "MSYNC_MG_events.csv")) %>%
  pivot_wider(names_from = c(phase, type), values_from = `timestamp [ns]`)

# read in the data, merge with sections and remove everything outside
df = merge(read_csv(file.path(dt.path[1], "fixations_on_face.csv")),
           df.sec) %>%
  # now we figure out to which phase which fixation belonged
  mutate(
    phase = case_when(
      `start timestamp [ns]` >= `1_clap` & `end timestamp [ns]` <= `1_end` ~ "p1",
      `start timestamp [ns]` >= `2_clap` & `end timestamp [ns]` <= `2_end` ~ "p2",
      `start timestamp [ns]` >= `3_clap` & `end timestamp [ns]` <= `3_end` ~ "p3"
    )
  ) %>%
  # filter out everything not happening during one of the phases
  filter(!is.na(phase)) %>%
  # now we compute the duration of the fixations in ms
  mutate(
    dur.ms = (`end timestamp [ns]` - `start timestamp [ns]`) / 1000000
  )

# Looking at the face of the other one ------------------------------------

# now we can aggregate the durations and compute dwell times for the individual
df.agg = df %>%
  mutate(dyad = substr(subID, 1, 8)) %>%
  group_by(subID, dyad, phase) %>%
  mutate(
    total.dur.ms = sum(dur.ms),
    ROI = if_else(`fixation on face`, 'face', 'other')
  ) %>%
  group_by(subID, dyad, phase, ROI, total.dur.ms) %>%
  summarise(
    dur.ms = sum(dur.ms)
  ) %>%
  mutate(
    dwell = dur.ms / total.dur.ms
  )

# save data frame
write_csv(df.agg, file.path(dt.path[2], "MSYNC_ET_MG_indi.csv"))

# Shared face contact -----------------------------------------------------

# when were two people looking at each other at the same time
df.face = df %>%
  filter(`fixation on face`) %>%
  mutate(
    start.aligned = case_when(
      phase == "p1" ~ `start timestamp [ns]` - `1_clap`,
      phase == "p2" ~ `start timestamp [ns]` - `2_clap`,
      phase == "p3" ~ `start timestamp [ns]` - `3_clap`),
    end.aligned = case_when(
      phase == "p1" ~ `end timestamp [ns]` - `1_clap`,
      phase == "p2" ~ `end timestamp [ns]` - `2_clap`,
      phase == "p3" ~ `end timestamp [ns]` - `3_clap`),
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
    comment = "Not two interaction partners"
    df.dyad = rbind(
      df.dyad, 
      data.frame(dyad = d, phase = 'p1', shared.ms = NA, comment), 
      data.frame(dyad = d, phase = 'p2', shared.ms = NA, comment), 
      data.frame(dyad = d, phase = 'p3', shared.ms = NA, comment)
    )
    next
  } else {
    comment = ''
  }
  # loop through phases
  for (ph in unique(df.face %>% filter(dyad == d) %>% select(phase))$phase) {
    # initialise empty list for the vectors
    t = list()
    # loop through the positions
    for (po in pos) {
      # focus on this one person
      df.sel = df.face %>% filter(position == po & dyad == d & phase == ph)
      if (nrow(df.sel) == 0) {
        comment = paste0(
          comment,  
          sprintf("%s: No relevant fixations for %s; ", ph, po))
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
      data.frame(dyad = d, phase = ph, shared.ms, comment)
    )
    comment = ''
  }
}

# add mean of total fixation duration for each dyad
df.dyad = merge(df.dyad, 
                df.agg %>% group_by(dyad, phase) %>% 
                  summarise(total.dur.ms.mean = mean(total.dur.ms))) %>%
  mutate(
    shared.perc = shared.ms/total.dur.ms.mean
  )

write_csv(df.dyad, file.path(dt.path[2], "MSYNC_ET_MG_dyad.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "ET_MG.Rdata"))

