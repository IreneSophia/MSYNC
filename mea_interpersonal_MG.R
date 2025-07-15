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
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = c("/Users/vilya/Documents/MSYNC/data/preprocessedMEA", 
              "/Users/vilya/Documents/MSYNC/data")
} else {
  dt.path = c("/media/emba/emba-2/MSYNC/data/preprocessedMEA", 
              "/media/emba/emba-2/MSYNC/data")
}

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
    dyad  = gsub(".*MEA/(.+)_MG_SC.*", "\\1", fln),
    phase = as.numeric(gsub(".*SC_p(.+)$", "\\1", gsub(".txt", "", fln))),
    flicker = flicker1 + flicker2
  ) %>% ungroup() %>% select(-fln) %>%
  # filter out the first second of data
  filter(frame > fps) %>%
  # exclude the frames where there is light flicker aka "paranormal activity"
  mutate(
    L      = if_else(flicker > 0, NA, L),
    R      = if_else(flicker > 0, NA, R)
  ) %>%
  group_by(dyad, phase) %>%
  arrange(dyad, phase, frame) %>%
  mutate(
    # use linear interpolation to replace the missing values
    L.ip = approx(frame, L, frame)$y,
    R.ip = approx(frame, R, frame)$y
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

# create a dataframe in which to put the information
df.ccf = data.frame()

# peak picking
for (i in 1:length(ls.ccf)){
  # drop rows with NAs
  all_lags = ls.ccf[[i]] %>% drop_na()
  idx.lag0 = which(colnames(all_lags) == "lag0")
  # extract information on positive lag
  R_peak = apply(all_lags[,(idx.lag0+1):ncol(all_lags)], 1, max, na.rm = T)
  R_mean = apply(all_lags[,(idx.lag0+1):ncol(all_lags)], 1, mean, na.rm = T)
  R_plag = abs(idx.lag0 - apply(all_lags[,(idx.lag0+1):ncol(all_lags)], 1, which.max))/fps
  # extract information on negative lag
  L_peak = apply(all_lags[,1:(idx.lag0-1)], 1, max, na.rm = T) 
  L_mean = apply(all_lags[,1:(idx.lag0-1)], 1, mean, na.rm = T) 
  L_plag = abs(idx.lag0 - apply(all_lags[,1:(idx.lag0-1)], 1, which.max))/fps
  # extract info of both lags
  B_mean = apply(all_lags, 1, mean, na.rm = T) 
  B_peak = apply(all_lags, 1, max, na.rm = T) 
  B_plag = abs(idx.lag0 - apply(all_lags, 1, which.max)) /fps
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
    L_total.mv = sum(L.ip > 0, na.rm = T)/ n(),
    R_total.mv = sum(R.ip > 0, na.rm = T) / n(),
    B_total.mv = sum(L.ip > 0 | R.ip > 0, na.rm = T)/ n()
  ) %>%
  pivot_longer(cols = where(is.numeric)) %>%
  separate(name, into = c("position", "AOI"), sep = "_") %>%
  pivot_wider(names_from = AOI, values_from = value)

# merge together
df = merge(df.ccf.agg, df.mov)

# save data frame
write_csv(df, file.path(dt.path[2], "MSYNC_mea_MG.csv"))

# Save workspace ----------------------------------------------------------

# save workspace
save.image(file = file.path(dt.path[1], "MEA_MG.Rdata"))

