# (C) Irene Sophia Plank
# 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
dt.path = c("/media/emba/emba-2/MSYNC/data/MSYNC_CT_ET_FACE-MAPPER", 
            "/media/emba/emba-2/MSYNC/data")

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
  group_by(subID) %>%
  mutate(
    total.dur.ms = sum(dur.ms),
    ROI = if_else(`fixation on face`, 'face', 'other')
  ) %>%
  group_by(subID, ROI, total.dur.ms) %>%
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
    start.ms = round(`start timestamp [ns]`/1000000),
    end.ms   = round(`end timestamp [ns]`/1000000),
    dyad     = substr(subID, 1, 8),
    position = substr(subID, nchar(subID), nchar(subID))
  )

# initialise dataframe
df.dyad = df.face %>%
  select(dyad) %>% 
  distinct() %>%
  mutate(
    shared.ms = NA,
    Comment = ""
  )

# loop through dyad
for (d in df.dyad$dyad) {
  # initialise list of vectors
  t = list()
  # get the positions of the associated interaction partners
  pos = unique(df.face %>% filter(dyad == d) %>% select(position))$position
  if (length(pos) != 2) {
    df.dyad[df.dyad$dyad == d,'Comment'] = "Not two interaction partners"
    next
  }
  # loop through the positions
  for (po in pos) {
    # focus on this one person
    df.sel = df.face %>% filter(position == po & dyad == d)
    if (nrow(df.sel) == 0) {
      df.dyad[df.dyad$dyad == d,'Comment'] = paste0(
        df.dyad[df.dyad$dyad == d,'Comment'], 
        sprintf("No relevant fixatios for %s; ", po))
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
  shared.ms = length(intersect(t$L, t$R))
  df.dyad[df.dyad$dyad == d, 'shared.ms'] = shared.ms
}

write_csv(df.dyad, file.path(dt.path[2], "MSYNC_ET_CT_dyad.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "ET_CT.Rdata"))

