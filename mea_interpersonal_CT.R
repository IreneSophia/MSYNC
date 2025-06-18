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

# length of time window of interest (seconds)
tiwo = 600

# seconds to ignore at the start
skip = 10

# Initialize function to create a fake MEA object out of two vectors
# Input: 
#     * s1, s2: numeric vectors containing the values to be correlated
#     * sampRate: sampling rate per second
#     * s1Name, s2Name: name for the values to be correlated, default is "s1Name" and "s2Name"
# Output:
#     * fake MEA object that pretends to be a MEA object
#
fakeMEA = function(s1, s2, sampRate, ROI, id) {
  mea = structure(list(all_01_01 = structure(list(MEA = structure(list(
    L = s1, R = s2), row.names = c(NA, -length(s1)), class = "data.frame"), 
    ccf = NULL, ccfRes = NULL), id = id, session = "01", group = ROI, sampRate = sampRate, 
    filter = "raw", ccf = "", s1Name = "L", s2Name = "R", uid = paste0(ROI, "_", id, "_01"), 
    class = c("MEA","list"))), class = "MEAlist", nId = 1L, n = 1L, groups = ROI, sampRate = sampRate, 
    filter = "raw", s1Name = "L", s2Name = "R", ccf = "")
  return(mea)
}

# Read in data ------------------------------------------------------------

# read in the data as a dataframe
df.mea = list.files(path = dt.path[1], pattern = "*CT_SC*", full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_delim(., show_col_types = F, delim = " ", 
                     col_names = c("L_head", "R_head", "L_body", "R_body", "flicker")), 
         .id = "fln") %>%
  group_by(fln) %>%
  mutate(
    frame = row_number(),
    dyad  = gsub(".*mea/(.+)_CT_SC.*", "\\1", fln)
  ) %>% ungroup() %>% select(-fln) %>%
  # exclude the frames where there is light flicker aka "paranormal activity"
  mutate(
    L_head = if_else(flicker > 0, NA, L_head),
    L_body = if_else(flicker > 0, NA, L_body),
    R_head = if_else(flicker > 0, NA, R_head),
    R_body = if_else(flicker > 0, NA, R_body),
    L      = L_head + L_body,
    R      = R_head + R_body
  ) %>%
  # only keep seconds 10-610
  filter(frame > fps*skip & frame <= fps*(tiwo+skip)) %>%
  group_by(dyad) %>%
  arrange(dyad, frame) %>%
  mutate(
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
  # extract relevant data
  df.sel = df.mea %>%
    filter(dyad == d)
  id  =  gsub("MSYNC_", "", d)
  ROI = 'all'
  # create fakeMEA object for this dyad
  mea.sel = fakeMEA(df.sel$L.ip, df.sel$R.ip, fps, ROI, id)
  names(mea.sel) = paste0(ROI, "_", id, "_01")
  # add it to the mea list
  mea = c(mea, mea.sel)
}

# Preprocessing -----------------------------------------------------------

# scaling
mea.scaled = MEAscale(mea)

# Time series synchronisation ---------------------------------------------

# compute windowed lagged cross correlation
mea.ccf = MEAccf(mea.scaled,
                 lagSec = 5,
                 winSec = 30, 
                 incSec = 15, 
                 r2Z = T,
                 ABS = T)

# visual inspection
pdf(file = paste(dt.path[1], "heatmaps_CT.pdf", sep = "/"))  
for (i in 1:length(mea.ccf)){
  MEAheatmap(mea.ccf[[i]], legendSteps = 20, rescale = T) 
}
dev.off()

# convert from mea list to list
ls.ccf = getCCF(mea.ccf, type = "fullMatrix")

# create a dataframe in which to put the information
df.ccf = data.frame()

# peak picking
for (i in 1:length(ls.ccf)){
  # drop rows with NAs
  all_lags = ls.ccf[[i]] %>% drop_na()
  idx.lag0 = which(colnames(all_lags) == "lag0")
  # extract information on positive lag (L movement happening before L movement)
  R_peak = apply(all_lags, 1, max, na.rm = T)
  R_mean = apply(all_lags, 1, mean, na.rm = T)
  R_plag = apply(all_lags, 1, which.max) + idx.lag0
  # extract information on negative lag (R movement happening before R movement)
  L_peak = apply(all_lags, 1, max, na.rm = T) 
  L_mean = apply(all_lags, 1, mean, na.rm = T) 
  L_plag = apply(all_lags, 1, which.max) 
  # extract info of both lags
  B_mean = apply(all_lags, 1, mean, na.rm = T) 
  B_peak = apply(all_lags, 1, max, na.rm = T) 
  B_plag = apply(all_lags, 1, which.max) 
  # extract lag0 synchrony
  B_zero = all_lags$lag0
  # add the information to the dataframe
  df.ccf = rbind(df.ccf, 
                 data.frame(R_peak, R_mean, R_plag, 
                            L_peak, L_mean, L_plag, 
                            B_peak, B_mean, B_plag, B_zero) %>% 
                   mutate(ID = names(ls.ccf)[i])
  )
}

# create one overall dataframe in the format ID-peaks
df.ccf = df.ccf %>% 
  pivot_longer(cols = where(is.numeric), names_to = "feature", 
               values_to = "MEA.sync") %>%
  separate("feature", sep = "_", into = c("position", "measure")) %>%
  separate("ID", sep = "_", into = c("ROI", "dyad", "phase")) %>%
  mutate(
    MEA.sync = if_else(MEA.sync != -Inf, MEA.sync, NA),
    dyad = paste0("MSYNC_", dyad),
    phase = as.numeric(phase)
  ) 

df.ccf.agg = df.ccf %>% 
  group_by(dyad, position, phase, measure) %>%
  summarise(
    MEA.sync = mean(MEA.sync, na.rm = T)
  )

# aggregate total movement and merge with ccf
df.mov = df.mea %>%
  group_by(dyad) %>%
  summarise(
    L = sum(L.ip > 0, na.rm = T)/ n(),
    R = sum(R.ip > 0, na.rm = T) / n(),
    B = sum(move, na.rm = T) / n()
  ) %>%
  pivot_longer(cols = c(L, R, B), names_to = "position", values_to = "MEA.mov")

# merge together
df = merge(df.ccf.agg, df.mov)

# save data frame
write_csv(df, file.path(dt.path[2], "MSYNC_mea_CT.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "MEA_CT.Rdata"))

