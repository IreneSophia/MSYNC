# (C) Irene Sophia Plank
# 
# This script takes the output of OpenFace and calculates interpersonal 
# synchrony. 

# clean workspace
rm(list = ls())

# load libraries
library(tidyverse)
library(rMEA)
library(data.table)    # setDT
library(moments)       # kurtosis, skewness

# set path to OpenFace files
dt.path = c("/Users/vilya/Documents/MSYNC/data/preprocessedOF", 
            "/Users/vilya/Documents/MSYNC/data")

# timing info
skip     = 10
fps      = 120
duration = 600

# set options
options(datatable.fread.datatable = F)

# Initialize function to create a fake MEA object out of two vectors
# Input: 
#     * s1, s2: numeric vectors containing the values to be correlated
#     * sampRate: sampling rate per second
#     * s1Name, s2Name: name for the values to be correlated, default is "s1Name" and "s2Name"
# Output:
#     * fake MEA object that pretends to be a MEA object
#
fakeMEA = function(s1, s2, sampRate, s1Name = "s1Name", s2Name = "s2Name") {
  mea = structure(list(all_01_01 = structure(list(MEA = structure(list(
    s1Name = s1, s2Name = s2), row.names = c(NA, -length(s1)), class = "data.frame"), 
    ccf = NULL, ccfRes = NULL), id = "01", session = "01", group = "all", sampRate = sampRate, 
    filter = "raw", ccf = "", s1Name = s1Name, s2Name = s2Name, uid = "all_01_01", 
    class = c("MEA","list"))), class = "MEAlist", nId = 1L, n = 1L, groups = "all", sampRate = sampRate, 
    filter = "raw", s1Name = s1Name, s2Name = s2Name, ccf = "")
  return(mea)
}

# Read in data ------------------------------------------------------------

if (!file.exists(file.path(dt.path[1], "MSYNC_OF_CT.rds"))) {
  # lists all files
  ls.fls = list.files(path = dt.path[1],
                      pattern = ".*aligned.csv")
  ls.IDs = gsub("(.+)_CT(.+)_a.*", "\\1\\2", ls.fls)
  
  # reads in the data
  ls      = lapply(file.path(dt.path[1], ls.fls), fread, 
                   header = F, # preserves columns for frame number
                   skip = fps*skip + 1, # skips start + header
                   nrows = duration*fps 
  ) 
  
  # read in header separately 
  header  = fread(file.path(dt.path[1], ls.fls[1]), 
                  header = F, nrows = 1, stringsAsFactors = F)
  for (i in 1:length(ls)){
    colnames(ls[[i]]) = unlist(header)
  }
  
  # add names to list elements
  names(ls) = ls.IDs
  
  # Inspect data ------------------------------------------------------------
  
  # check for NAs 
  for (i in 1:length(ls)){
    if (anyNA(ls[[i]]) == T){ 
      warning(sprintf('%s has NAs', names(ls)[i]))
    }
  }
  
  # check if every dataframe has correct number of lines
  for (i in 1:length(ls)){
    if (nrow(ls[[i]]) != fps*duration){ 
      warning(sprintf('%s is incomplete: %i', names(ls)[i], nrow(ls[[i]])))
    }
  }
  
  # MSYNC_24_FL is incomplete: 47007 
  
  # filter data based on mean confidence and successfully tracked frames 
  outlier = c()
  for (i in 1:length(ls)) {
    if (mean(ls[[i]]$confidence) < 0.75 # mean confidence of tracked frames lower than 75%?
        || # OR
        sum(ls[[i]]$success) < fps*duration*0.9) # less than 90% of goal frames successfully tracked?
    { 
      warning(sprintf('face tracking not reliable enough for %s', names(ls)[i]))
      outlier = c(outlier, gsub("(.+)_F.*", "\\1", names(ls[i]))) # add dyad to outliers
    }
  }
  
  # face tracking not reliable enough for MSYNC_24_FL, this dyad is excluded
  ls[which(grepl(outlier, names(ls)))] = NULL
  
  # convert to dataframe
  df = bind_rows(ls, .id = "ID") %>%
    separate(col = ID, into = c("dyad1", "dyad2", "speaker"), remove = F) %>%
    mutate(dyad = paste(dyad1, dyad2, sep = "_"),
           ID = paste0(dyad, "_", speaker)) %>%
    select(-dyad1, -dyad2) %>%
    relocate(ID, dyad, speaker)
  
  saveRDS(df, file.path(dt.path[1], "MSYNC_OF_CT.rds"))
} else {
  df = readRDS(file.path(dt.path[1], "MSYNC_OF_CT.rds"))
}

# clean workspace
rm(list = setdiff(ls(), c("df", "fakeMEA", "dt.path", "fps")))

# list of AUs
ls.AUs = str_subset(names(df), "AU.*_r|pose_R*")

# Time series synchronisation ---------------------------------------------

if (!file.exists(file.path(dt.path[1], "MSYNC_OF_CT_sync.rds"))) {
  
  # Steps: 
  # 1) create fake MEA object using the fakeMEA function
  # 2) calculate ccf according to rMEA
  
  # create list to be filled with fakeMEA objects of AUs
  ls.fakeAU = c()
  
  # loop through all dyads
  for (i in unique(df$dyad)){ 
    
    # initialise heatmaps
    dir.create(file.path(dt.path[1], 'pics'), showWarnings = FALSE)
    pdf(file.path(dt.path[1], 'pics', paste0(i, ".pdf")))
    
    ## START
    
    # grab only relevant portions of df
    df.sel = df %>%
      filter(dyad == i) %>%
      arrange(ID, dyad, speaker, frame)
    
    # check if data frame is present
    if (nrow(df.sel) > 0) { 
      
      print(i)
      
      # loop over AUs
      for (j in ls.AUs){ 
        
        # prepare fake MEA components
        s1 = df.sel[df.sel$speaker == "FL", j] # AU for left L participant
        s2 = df.sel[df.sel$speaker == "FR", j] # AU for right R participant
        
        # create fake MEA object
        mea = fakeMEA(s1, s2, fps) 
        
        # time lagged windowed cross-correlations
        mea = MEAccf(mea, lagSec = 2, winSec = 7, incSec = 4, r2Z = T, ABS = T) 
        names(mea) = paste(i, "CT", j, sep = "_")
        
        # add object to fakeMEA list
        ls.fakeAU = c(ls.fakeAU, mea)
        
        # configure heatmap
        par(col.main='white')                  # set plot title to white
        heatmap = MEAheatmap(mea[[1]])
        par(col.main='black')                  # set plot title back to black
        title(main = paste("CT", j, sep = "_")) # alternative title
        
      }
    }
    
    dev.off()
    # show progress
    print(paste(i, "done"))
  }
  
  saveRDS(ls.fakeAU, file = file.path(dt.path[1], "MSYNC_OF_CT_sync.rds"))
  
} else {
  
  ls.fakeAU = readRDS(file.path(dt.path[1], "MSYNC_OF_CT_sync.rds"))
  
}

# clean workspace
rm(list = setdiff(ls(), c("df", "fakeMEA", "ls.fakeAU", "fps",
                          "dt.path", "ls.AUs")))

# convert from mea list to list
ls.ccf = getCCF(ls.fakeAU, type = "fullMatrix")
names(ls.ccf) = names(ls.fakeAU)

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


# create one overall dataframe in the long format
df.ccf = df.ccf %>% 
  pivot_longer(cols = where(is.numeric), names_to = "feature", 
               values_to = "OF.sync") %>%
  separate("feature", sep = "_", into = c("position", "measure")) %>%
  mutate(
    OF.sync = if_else(OF.sync != -Inf, OF.sync, NA),
    dyad  = gsub("(.+)_CT.*", "\\1", ID),
    phase = "CT",
    key   = gsub(".*_CT_(.+)", "\\1", ID)
  ) 

# aggregate the ccf values
df.ccf.agg = df.ccf %>% 
  group_by(dyad, position, phase, key, measure) %>%
  summarise(
    OF.sync = mean(OF.sync, na.rm = T)
  )

saveRDS(df.ccf.agg, file.path(dt.path[2], "MSYNC_AU_sync_CT.rds"))

# Facial expressiveness ---------------------------------------------------

# full expressiveness: mean of all AUs and frames
df.exp = df %>%
  select(ID, dyad, speaker, frame, matches("AU.*r")) %>%
  pivot_longer(names_to = "input", values_to = "exp", cols = matches("AU.*r")) %>%
  group_by(ID, dyad, speaker) %>%
  summarise(
    exp = mean(exp, na.rm = T)
  )

# save facial expressiveness
saveRDS(df.exp, file.path(dt.path[2], "MSYNC_AU_intensity_CT.rds"))

# Save workspace ----------------------------------------------------------

# clean workspace
rm(list = setdiff(ls(), c("df", "ls.fakeAU", "ls.AUs", "dt.path", "df.ccf")))

# save workspace
save.image(file = file.path(dt.path[1], "OpenFace_CT.RData"))

