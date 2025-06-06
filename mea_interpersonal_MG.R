# (C) Irene Sophia Plank
# 
# This script takes the output of Motion Energy analysis and calculates INTER-
# personal synchrony. It is an adaptation of a script written by Jana Koehler, 
# published in https://github.com/jckoe/MLASS-study. 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)
library(rMEA)

# set path to MEA files
dt.path = c("/media/emba/emba-2/MSYNC/data/mea", 
            "/media/emba/emba-2/MSYNC/data")

# set frame rate
fps = 120

# Initialize function to create a fake MEA object out of two vectors
# Input: 
#     * s1, s2: numeric vectors containing the values to be correlated
#     * sampRate: sampling rate per second
#     * s1Name, s2Name: name for the values to be correlated, default is "s1Name" and "s2Name"
# Output:
#     * fake MEA object that pretends to be a MEA object
#
fakeMEA = function(s1, s2, sampRate, ROI, id, s) {
  mea = structure(list(all_01_01 = structure(list(MEA = structure(list(
    L = s1, R = s2), row.names = c(NA, -length(s1)), class = "data.frame"), 
    ccf = NULL, ccfRes = NULL), id = id, session = sprintf("%02d", s), group = ROI, sampRate = sampRate, 
    filter = "raw", ccf = "", s1Name = "L", s2Name = "R", uid = sprintf("%s_%s_%02d", ROI, id, i), 
    class = c("MEA","list"))), class = "MEAlist", nId = 1L, n = 1L, groups = ROI, sampRate = sampRate, 
    filter = "raw", s1Name = "L", s2Name = "R", ccf = "")
  return(mea)
}

# Read in data ------------------------------------------------------------

# read in the data as a dataframe
df.mea = list.files(path = dt.path[1], pattern = "*_MG_SC_p*", full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_delim(., show_col_types = F, delim = " ", 
                     col_names = c("L", "R", "lamp", "flicker1", "flicker2")), 
         .id = "fln") %>%
  group_by(fln) %>%
  mutate(
    frame = row_number(),
    dyad  = gsub(".*mea/(.+)_MG_SC.*", "\\1", fln),
    phase = as.numeric(gsub(".*SC_p(.+)$", "\\1", gsub(".txt", "", fln)))
  ) %>% ungroup() %>% select(-fln) %>%
  # exclude participants with "no" in name
  filter(!is.na(phase)) %>%
  # exclude the frames where there is light flicker aka "paranormal activity"
  mutate(
    L = if_else(flicker1 > 0 | flicker2 > 0, NA, L),
    R = if_else(flicker1 > 0 | flicker2 > 0, NA, R)
  ) %>%
  group_by(dyad, phase) %>%
  arrange(dyad, phase, frame) %>%
  mutate(
    # check how much each dyad lost to the flicker
    total_flick = mean(flicker1 > 0 | flicker2 > 0),
    # use linear interpolation to replace the missing values
    L.ip = approx(frame, L, frame)$y,
    R.ip = approx(frame, R, frame)$y,
    # check whether one of both people moved
    move = if_else(L.ip > 0 | R.ip > 0, 1, 0)
    )

# check how much data is still lost despite interpolation
df.miss = df.mea %>% 
  group_by(dyad) %>%
  summarise(
    missing = round(mean(is.na(R.ip) | is.na(L.ip)),3)
  )

# initialise mea list
mea = c()

# loop through the dyads
for (d in unique(df.mea$dyad)) {
  for (i in 1:3) {
    # extract relevant data
    df.sel = df.mea %>%
      filter(dyad == d & phase == i)
    if (nrow(df.sel) > 0) {
      id  =  gsub("MSYNC_", "", d)
      ROI = 'all'
      # create fakeMEA object for this dyad
      mea.sel = fakeMEA(df.sel$L.ip, df.sel$R.ip, fps, ROI, id, i)
      names(mea.sel) = sprintf("%s_%s_%02d", ROI, id, i)
      # add it to the mea list
      mea = c(mea, mea.sel)
    }
  }
}

# Preprocessing -----------------------------------------------------------

# scaling
mea.scaled = MEAscale(mea)

# Time series synchronisation ---------------------------------------------

# compute windowed lagged cross correlation
mea.ccf = MEAccf(mea.scaled,
                 lagSec = 2,
                 winSec = 30, 
                 incSec = 1, 
                 r2Z = T,
                 ABS = T)

# visual inspection
pdf(file = paste(dt.path[1], "heatmaps_MG.pdf", sep = "/"))  
for (i in 1:length(mea.ccf)){
  MEAheatmap(mea.ccf[[i]], legendSteps = 20, rescale = T) 
}
dev.off()

# convert from mea list to list
ls.ccf = getCCF(mea.ccf, type = "fullMatrix")

# peak picking
for (i in 1:length(ls.ccf)){
  # append maximum of positive lag (L movement happening before L movement)
  ls.ccf[[i]]$R_peak = apply(ls.ccf[[i]][,(floor(ncol(ls.ccf[[i]])/2)+2):ncol(ls.ccf[[i]])], 1, max, na.rm = T)
  ls.ccf[[i]]$R_mean = apply(ls.ccf[[i]][,(floor(ncol(ls.ccf[[i]])/2)+2):ncol(ls.ccf[[i]])], 1, mean, na.rm = T)
  # append maximum of negative lag (R movement happening before R movement)
  ls.ccf[[i]]$L_peak = apply(ls.ccf[[i]][,1:floor(ncol(ls.ccf[[i]])/2)], 1, max, na.rm = T) 
  ls.ccf[[i]]$L_mean = apply(ls.ccf[[i]][,1:floor(ncol(ls.ccf[[i]])/2)], 1, mean, na.rm = T) 
  # append maximum of both lags
  ls.ccf[[i]]$B_peak = apply(ls.ccf[[i]], 1, max, na.rm = T) 
  ls.ccf[[i]]$B_mean = apply(ls.ccf[[i]], 1, mean, na.rm = T) 
  # keep only relevant columns
  ls.ccf[[i]] = ls.ccf[[i]][,c("L_peak", "L_mean","R_peak", "R_mean", "B_peak", "B_mean")]
  # transpose all df in list
  ls.ccf[[i]] = as.data.frame(t(ls.ccf[[i]]))
  # set rownames as first column
  data.table::setDT(ls.ccf[[i]], keep.rownames = TRUE)
  colnames(ls.ccf[[i]])[1] = "feature"
}

# create one overall dataframe in the format ID-peaks
df.ccf = bind_rows(ls.ccf, .id = "ID") %>% 
  separate(ID, c("ROI", "dyad", "phase")) %>%
  separate(feature, c("position", "aggregation")) %>%
  pivot_longer(cols = starts_with("w"), names_to = "window", values_to = "MEA.sync") %>%
  mutate(
    MEA.sync = if_else(MEA.sync != -Inf, MEA.sync, NA),
    dyad = paste0("MSYNC_", dyad),
    phase = as.numeric(phase)
  ) %>% 
  group_by(dyad, position, phase, aggregation) %>%
  summarise(
    MEA.sync = mean(MEA.sync, na.rm = T)
  )

# aggregate total movement and merge with ccf
df.mov = df.mea %>%
  group_by(dyad, phase) %>%
  summarise(
    L = sum(L.ip > 0, na.rm = T)/ n(),
    R = sum(R.ip > 0, na.rm = T) / n(),
    B = sum(move, na.rm = T) / n()
  ) %>%
  pivot_longer(cols = c(L, R, B), names_to = "position", values_to = "MEA.mov")

# merge together
df = merge(df.ccf, df.mov)

# save data frame
write_csv(df, file.path(dt.path[2], "MSYNC_mea_MG.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "MEA_MG.Rdata"))

