# (C) Irene Sophia Plank
# 
# This script takes the output of Motion Energy analysis and calculates INTER-
# personal synchrony. 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)
library(rMEA)

# set path to OpenPose files
if (Sys.getenv("LOGNAME") == "vilya") {
  dt.path = c("/Users/vilya/Documents/MSYNC/data/preprocessedOP", 
              "/Users/vilya/Documents/MSYNC/data")
} else {
  dt.path = c("/media/emba/emba-2/MSYNC/data/preprocessedOP", 
              "/media/emba/emba-2/MSYNC/data")
}

# set the task 
task = "MG"

# set frame rate
fps = 120

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

df.ref = readRDS(file.path(dt.path[1], sprintf("MSYNC_OP_%s_ref.rds", task)))

# initialise "mea" list
mea = c()

for (d in unique(df.ref$dyad)) {
  
  # select the relevant data
  df.sel = df.ref %>%
    filter(dyad == d)

  # extract relevant info
  id  =  gsub("MSYNC_(.+)", "\\1", d)
  
  # loop through the keypoints
  for (k in unique(df.sel$key)) {
    # extract the relevant datapoints: mirrored for left/right
    if (grepl("_R", k)) {
      L = df.sel[df.sel$key == k,]$L
      R = df.sel[df.sel$key == gsub("R", "L", k),]$R
    } else if (grepl("_L", k)) {
      L = df.sel[df.sel$key == k,]$L
      R = df.sel[df.sel$key == gsub("L", "R", k),]$R
    } else {
      L = df.sel[df.sel$key == k,]$L
      R = df.sel[df.sel$key == k,]$R
    }
    # create fakeMEA object for this dyad and keypoint
    mea.sel = fakeMEA(L, R, fps, k, id)
    names(mea.sel) = paste0(k, "_", id, "_01")
    # add it to the mea list
    mea = c(mea, mea.sel)
  }
}

# Preprocessing -----------------------------------------------------------

# scaling
mea.scaled = MEAscale(mea)

# Time series synchronisation ---------------------------------------------

# compute windowed lagged cross correlation
if (task == "CT") {
  mea.ccf = MEAccf(mea.scaled,
                   lagSec = 5,
                   winSec = 30, 
                   incSec = 15, 
                   r2Z = T,
                   ABS = T)
} else {
  mea.ccf = MEAccf(mea.scaled,
                   lagSec = 2,
                   winSec = 30, 
                   incSec = 1, 
                   r2Z = T,
                   ABS = T)
}

# visual inspection
pdf(file = file.path(dt.path[1], sprintf("heatmaps_%s.pdf", task)))  
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

# create one overall dataframe
df.ccf = df.ccf %>% 
  pivot_longer(cols = where(is.numeric), names_to = "feature", 
               values_to = "IPSmov") %>%
  separate("feature", sep = "_", into = c("position", "measure")) %>%
  separate("ID", sep = "_", into = c("key", "dyad", "phase")) %>%
  mutate(
    IPSmov = if_else(IPSmov != -Inf, IPSmov, NA),
    dyad = paste0("MSYNC_", dyad)
  ) 

if (task == "CT") {
  df.ccf$phase = "CT"
}

# aggregate the synchrony values
df.ccf.agg = df.ccf %>% 
  group_by(dyad, key, position, phase, measure) %>%
  summarise(
    IPSmov = mean(IPSmov, na.rm = T)
  )

# load information on movement in general
df.mov = readRDS(file.path(dt.path[1], sprintf("MSYNC_OP_%s.rds", task))) %>%
  rename("position" = "side") %>%
  select(dyad, position, key, frame, dist)

# summarise the movement quantity
df.mov.agg = rbind(
  df.mov %>% 
    mutate(
      key = case_when(
        grepl("L", key) & position == "R" ~ gsub("L", "R", key),
        grepl("R", key) & position == "R" ~ gsub("R", "L", key),
        T ~ key
      ),
      position = "B"
    ),
  df.mov) %>%
  filter(key != "ref") %>%
  group_by(dyad, position, key) %>%
  summarise(
    QNTmov = sum(dist, na.rm = T)
  )

# merge the dataframes
df = merge(
  df.mov.agg,
  df.ccf.agg)

# save data frame
write_csv(df, file.path(dt.path[2], sprintf("MSYNC_OP_%s.csv", task)))

# Save workspace ----------------------------------------------------------

# clean workspace
rm(list = setdiff(ls(), c("mea", "mea.ccf", "dt.path", "fps", "task")))

# save workspace
save.image(file = file.path(dt.path[1], sprintf("MSYNC_OP_%s.Rdata", task)))

