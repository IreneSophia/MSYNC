# Script to extract some subject-specific information, including their
# dyad and order of tasks
# (C) Irene Sophia Plank
# 
# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# set path to MEA files
dt.path = c("/media/emba/emba-2/MSYNC/data")

# goal: data frame with the PID, dyad, subID (in MSYNC), order of the tasks

# get the order, PID, dyad and position of each participant
df = merge(
  read_csv(file.path(dt.path, "MSYNC_ParticipantLog.csv")) %>%
    rename("order" = `Order of Tasks`, "dyad" = "Dyad Number") %>%
    select(dyad, PID, order), 
  read_csv(file.path(dt.path, "df_centraXX_final.csv")) %>%
    rename("PID" = "internalStudyMemberID") %>%
    select(PID, dyad, role)) %>%
  mutate(
    dyad     = sprintf("MSYNC_%02d", dyad),
    position = if_else(role == "Leader1", "L", "R"), 
    subID    = paste0(dyad, "_", position)
  )

write_csv(df, file.path(dt.path, "MSYNC_subID_info.csv"))
