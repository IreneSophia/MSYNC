# (C) Irene Sophia Plank
# 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
dt.path = c("/media/emba/emba-2/MSYNC/data/MSYNC_MG_ET_FACE-MAPPER", 
            "/media/emba/emba-2/MSYNC/data")

# Read in data ------------------------------------------------------------

# read in the sections as a dataframe
df.sec = read_csv(file.path(dt.path[1], "sections.csv")) %>%
  select(`section id`, 
         `section start time [ns]`, `section end time [ns]`) %>%
  pivot_longer(cols = where(is.numeric), names_to = "trigger")

##### [!ADAPT]: here are only the recording start and end, where phases?

# read in the data, merge with sections and remove everything outside
df = read_csv(file.path(dt.path[1], "fixations_on_face.csv"))

# merge and find start and end point
df.edge = merge(
  df %>%
    select(`section id`, 
           `start timestamp [ns]`) %>%
    pivot_longer(cols = where(is.numeric)) %>%
    merge(., df.sec, all = T) %>%
    fill(trigger, .direction = "down") %>%
    # focus on all the data within the section
    filter(trigger == "section start time [ns]") %>%
    drop_na() %>% group_by(`section id`) %>% 
    mutate(
      row = row_number()
    ) %>% filter(row == min(row)) %>%
    rename("start" = "value") %>%
    select(`section id`, start),
  df %>%
    select(`section id`, 
           `end timestamp [ns]`) %>%
    pivot_longer(cols = where(is.numeric)) %>%
    merge(., df.sec, all = T) %>%
    fill(trigger, .direction = "down") %>%
    # focus on all the data within the section
    filter(trigger == "section start time [ns]") %>%
    drop_na() %>% group_by(`section id`) %>% 
    mutate(
      row = row_number()
    ) %>% filter(row == max(row)) %>%
    rename("end" = "value") %>%
    select(`section id`, end)
  )

# merge the starting and end points with the data
df = merge(df, df.edge) %>%
  # now we finally filter out the fixations outside the section
  filter(
    `start timestamp [ns]` >= start & `end timestamp [ns]` <= end
  ) %>%
  # now we compute the duration of the fixations in ms
  mutate(
    dur.ms = (`end timestamp [ns]` - `start timestamp [ns]`) / 1000000
  )

# now we can aggregate the durations and compute dwell times
df.agg = df %>%
  group_by(`section id`) %>%
  mutate(
    total.dur.ms = sum(dur.ms),
    ROI = if_else(`fixation on face`, 'face', 'other')
  ) %>%
  group_by(`section id`, ROI, total.dur.ms) %>%
  summarise(
    dur.ms = sum(dur.ms)
  ) %>%
  mutate(
    dwell = dur.ms / total.dur.ms
  )

# save data frame
write_csv(df.agg, file.path(dt.path[2], "MSYNC_ET_CT.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "ET_CT.Rdata"))

