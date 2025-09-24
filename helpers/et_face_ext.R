# Extract features from the mobile eye-tracking data of the MSYNC project:
# Face attention, shared face attention, eye contact initiation
# (C) Irene Sophia Plank
# 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)

# select the task
task = "MG"

# set path to MEA files
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = "/Users/vilya/Documents/MSYNC/data/"
} else {
  dt.path = "/media/emba/emba-2/MSYNC/data"
}

# sampling rate
fps = 200

# Read in data ------------------------------------------------------------

# differs between the tasks
if (task == "CT") {
  
  # select time window length in seconds
  tiwo = 600
  # read in the sections as a dataframe
  df.sec = read_csv(file.path(dt.path, sprintf("MSYNC_%s_events.csv", task))) %>%
    pivot_wider(names_from = name, values_from = `timestamp [ns]`) %>%
    mutate(
      # shift the end trigger to keep specific tiwo
      end = clap + 1000000000 * tiwo,
      # if the recording is too short, then use the full recording
      end = if_else(end > recording.end, recording.end, end)
    ) %>% select(-recording.end)
  
  # read in the data, merge with sections and remove everything outside
  df = merge(read_csv(file.path(dt.path, sprintf("MSYNC_%s_ET_FACE-MAPPER", task), "gaze_on_face.csv")), 
             df.sec) %>%
    # now we finally filter out the fixations outside the section
    filter(
      `timestamp [ns]` >= clap & `timestamp [ns]` <= end
    ) %>%
    group_by(subID) %>% arrange(`timestamp [ns]`) %>%
    mutate(
      dyad     = substr(subID, 1, 8), 
      t.adjust = `timestamp [ns]` - min(`timestamp [ns]`),
      t        = row_number(),
      dur      = max(t.adjust)/1000000000
    ) 
  
} else if (task == "MG") {
  
  # normal length of phases
  tiwo = 3*60
  
  # read in the sections as a dataframe
  df.sec = read_csv(file.path(dt.path, sprintf("MSYNC_%s_events.csv", task))) %>%
    pivot_wider(names_from = c(phase, type), values_from = `timestamp [ns]`)
  
  # read in the data, merge with sections and remove everything outside
  df = merge(read_csv(file.path(dt.path, sprintf("MSYNC_%s_ET_FACE-MAPPER", task), "gaze_on_face.csv")),
             df.sec) %>%
    # now we figure out to which phase which fixation belonged
    mutate(
      phase = case_when(
        `timestamp [ns]` >= `1_clap` & `timestamp [ns]` <= `1_end` ~ "p1",
        `timestamp [ns]` >= `2_clap` & `timestamp [ns]` <= `2_end` ~ "p2",
        `timestamp [ns]` >= `3_clap` & `timestamp [ns]` <= `3_end` ~ "p3"
      )
    ) %>%
    # filter out everything not happening during one of the phases
    filter(!is.na(phase)) %>%
    # arrange it
    group_by(subID, phase) %>% arrange(`timestamp [ns]`) %>%
    mutate(
      dyad     = paste0(substr(subID, 1, 8), "_", phase), 
      t.adjust = `timestamp [ns]` - min(`timestamp [ns]`),
      t        = row_number(),
      dur      = max(t.adjust)/1000000000
    )
  
}

# check the length: only "MSYNC_03_R" in CT too short
df = df %>%
  # only keep participants where at least 2/3 were tracked
  filter(dur > 2*tiwo/3) %>% ungroup()

# get a list of all dyads where both are included
df.inc = df %>% select(subID, dyad) %>% distinct() %>%
  group_by(dyad) %>% count() %>% filter(n == 2)

# Face attention: individual and shared -----------------------------------

df.wide = df %>% ungroup() %>%
  select(subID, dyad, t, `gaze on face`) %>%
  mutate(
    side  = substr(subID, nchar(subID), nchar(subID))
  ) %>% select(-subID) %>% 
  pivot_wider(names_from = side, values_from = `gaze on face`) %>%
  mutate(
    B = L & R
  )

# aggregate for the individual
df.indi = df.wide %>% 
  group_by(dyad) %>%
  summarise(
    L = mean(L, na.rm = T)*100,
    R = mean(R, na.rm = T)*100
  ) %>%
  pivot_longer(cols = c(L, R), names_to = "side", values_to = "dwell.face") %>%
  mutate(
    subID = paste0(substr(dyad, 1, 8), "_", side), side = NULL
  ) %>% ungroup()

# aggregate for the dyad > shared face attention
df.dyad = df.wide %>%
  # exclude dyads with only one participant
  filter(dyad %in% df.inc$dyad) %>%
  group_by(dyad) %>%
  summarise(shared.face = mean(B, na.rm = T)*100)

# Initiations -------------------------------------------------------------

# initialise dataframe
df.ini = data.frame()

# focus on the rows where something changed
df.change = df.wide %>% 
  # exclude dyads with only one participant
  filter(dyad %in% df.inc$dyad) %>% 
  fill(c(L, R), .direction = "down") %>%
  select(-B) %>% drop_na() %>%
  filter(L != lag(L) | R != lag(R))

# loop through the dyads
for (d in unique(df.change$dyad)) {
  print(sprintf('%s : %s', Sys.time(), d))
  currentStatus = startRow = initiator = F
  df.sel = df.change %>% filter(dyad == d)
  # loop through the changed rows
  for (i in 1:nrow(df.sel)) {
    # is our status that we have found one? If no...
    if (currentStatus == F) {
      # check whether we have found one now!
      if (df.sel$L[i] == T & df.sel$R[i] == F) {
        startRow = df.sel$t[i]
        initiator = "L"
        currentStatus = T
      } else if (df.sel$L[i] == F & df.sel$R[i] == T) {
        startRow = df.sel$t[i]
        initiator = "R"
        currentStatus = T
      }
      # if we had already found one, check if it is now successful
    } else if (df.sel$L[i] == T & df.sel$R[i] == T){
      df.ini = rbind(df.ini, 
                     data.frame(dyad = d, startRow, endRow = df.sel$t[i], initiator, cat = "SI"))
      currentStatus = startRow = initiator = F # start new
      # if we had already found one, check if it failed
    } else if (df.sel$L[i] == F & df.sel$R[i] == F) {
      df.ini = rbind(df.ini, 
                     data.frame(dyad = d, startRow, endRow = df.sel$t[i], initiator, cat = "FI"))
      currentStatus = startRow = initiator = F # start new
    }
  }
}

# aggregate the information and add to the individual dataframe
df.indi = merge(df.indi, 
                df.ini %>%
                  group_by(dyad, initiator, cat) %>%
                  count() %>%
                  pivot_wider(names_from = cat, values_from = n) %>%
                  ungroup() %>%
                  mutate(
                    subID = paste0(substr(dyad, 1, 8), "_", initiator), 
                    initiator = NULL,
                    TI = SI + FI,
                    SI.ratio = SI / TI,
                    SI.FI    = SI / FI
                  ), all = T)

# aggregate the information and add to the dyad dataframe
df.dyad = merge(df.dyad, 
                df.ini %>%
                  group_by(dyad, cat) %>%
                  count() %>%
                  pivot_wider(names_from = cat, values_from = n) %>%
                  ungroup() %>%
                  mutate(
                    TI = SI + FI,
                    SI.ratio = SI / TI,
                    SI.FI    = SI / FI
                  ), all = T)

# if MG, then separate dyad and phase again
if (task == "MG") {
  df.indi = df.indi %>%
    mutate(
      phase = substr(dyad, nchar(dyad)-1, nchar(dyad)),
      dyad  = substr(dyad, 1, 8)
    )
  df.dyad = df.dyad %>%
    mutate(
      phase = substr(dyad, nchar(dyad)-1, nchar(dyad)),
      dyad  = substr(dyad, 1, 8)
    )
}

# Save workspace ----------------------------------------------------------

# save the csvs
write_csv(df.dyad, file.path(dt.path, sprintf("MSYNC_ET_%s_dyad.csv", task)))

# save the csvs
write_csv(df.indi, file.path(dt.path, sprintf("MSYNC_ET_%s_indi.csv", task)))

# save workspace
save.image(file = file.path(dt.path, sprintf("MSYNC_%s_ET_FACE-MAPPER", task), 
                            sprintf("ET_%s.Rdata", task)))

