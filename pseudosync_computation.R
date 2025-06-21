# (C) Irene Sophia Plank
# 
# This script takes lists of MEA and fakeMEA objects and calculates pseudo-
# synchrony using the pseudosync_seg function. 

library(tidyverse)
library(rMEA)

# clean workspace
rm(list = ls())

# load pseudosync_seg function
source("pseudosync_seg.R")

# set paths
dt.path  = "/Users/vilya/Documents/MSYNC/data/mea"

# set number of permutations
n = 50

# Means -------------------------------------------------------------------

# create a dataframe
df.psync = data.frame()

# start with the conversation task 
attach(file.path(dt.path, "MEA_CT.Rdata"), name = 'CT')
df.mea = pseudosync_seg(mea.scaled, 
                        sampRate = 120, 
                        lagSec = 5,
                        winSec = 30, 
                        incSec = 15, 
                        n = n)
df.mea$task = "CT"
df.psync = rbind(df.psync, df.mea)
# detach the CT data
detach('CT')

# at the mirror game  
attach(file.path(dt.path, "MEA_MG.Rdata"), name = 'MG')
df.mea = pseudosync_seg(mea.scaled, 
                        sampRate = 120, 
                        lagSec = 2,
                        winSec = 30, 
                        incSec = 1, 
                        n = n)
df.mea$task = "MG"
df.psync = rbind(df.psync, df.mea)
# detach the MG data
detach('MG')

# save the pseudosynchrony values
write_csv(df.psync, file.path(dt.path, "df_psync.csv"))

# Peaks -------------------------------------------------------------------


# create a dataframe
df.psync = data.frame()

# start with the conversation task 
attach(file.path(dt.path, "MEA_CT.Rdata"), name = 'CT')
df.mea = pseudosync_seg(mea.scaled, 
                        sampRate = 120, 
                        lagSec = 5,
                        winSec = 30, 
                        incSec = 15, 
                        n = n, peak = T)
df.mea$task = "CT"
df.psync = rbind(df.psync, df.mea)
# detach the CT data
detach('CT')

# at the mirror game  
attach(file.path(dt.path, "MEA_MG.Rdata"), name = 'MG')
df.mea = pseudosync_seg(mea.scaled, 
                        sampRate = 120, 
                        lagSec = 2,
                        winSec = 30, 
                        incSec = 1, 
                        n = n, peak = T)
df.mea$task = "MG"
df.psync = rbind(df.psync, df.mea)
# detach the MG data
detach('MG')

# save the pseudosynchrony values
write_csv(df.psync, file.path(dt.path, "df_psync_peaks.csv"))
