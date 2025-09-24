# Summarise the triggers for the eye-tracking data.
# (C) Irene Sophia Plank
# 
# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
dt.path = c("/media/emba/emba-2/MSYNC/data")

# summarise the triggers for the CT task
list.files(path = file.path(dt.path, "ET_events"), 
                pattern = "*_CT_ET*", full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(., show_col_types = F), 
         .id = "fln") %>%
  mutate(
    dyad     = gsub(".*ET_events/(.+)_CT_ET.*", "\\1", fln),
    position = gsub(".*_CT_ET_(.+)-events.csv", "\\1", fln),
    subID    = paste0(dyad, "_", position)
  ) %>%
  # check if everyone has a clap
  group_by(subID) %>%
  mutate(
    clap = sum(name == "clap")
  ) %>%
  # filter out the people without a clap and the recording start
  filter(clap == 1 & name != "recording.begin") %>%
  select(subID, `recording id`, name, `timestamp [ns]`) %>%
  write_csv(., file.path(dt.path, "MSYNC_CT_events.csv"))


# summarise the triggers for the MG task
list.files(path = file.path(dt.path, "ET_events"), 
           pattern = "*_MG_ET*", full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(., show_col_types = F), 
         .id = "fln") %>%
  mutate(
    dyad     = gsub(".*ET_events/(.+)_MG_ET.*", "\\1", fln),
    position = gsub(".*_MG_ET_(.+)-events.csv", "\\1", fln),
    subID    = paste0(dyad, "_", position),
    phase    = as.numeric(gsub(".*_p(.+)", "\\1", name)),
    type     = if_else(substr(name, 1, 3) == "End", "end", "clap")
  ) %>%
  # filter out everything not associated with a phase
  filter(!is.na(phase)) %>%
  # check if everyone has a clap
  group_by(subID, phase) %>%
  mutate(
    start = sum(type == "clap"),
    end   = sum(type == "end"),
    # check if end is after the clap
    diff = `timestamp [ns]` - lag(`timestamp [ns]`)
  ) %>%
  # filter out the people where one or both are missing or end is before clap
  filter(start == 1 & end == 1 & (diff > 0 | is.na(diff))) %>%
  select(subID, `recording id`, phase, type, `timestamp [ns]`) %>%
  write_csv(., file.path(dt.path, "MSYNC_MG_events.csv"))
